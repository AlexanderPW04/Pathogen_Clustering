"""Hosted (or local) end-to-end check of the website, driven like a browser.

Each scenario uploads files to /review, reads the review page's own form
(roles, columns, and settings exactly as a visitor would see them), applies a
few documented choices, presses Run (/run), and downloads every result file.
The same script is run against a local server (baseline) and against a staging
URL, so the two can be compared file by file.

Usage:
    python validation/hosted_smoke_test.py --base http://127.0.0.1:5055 --output <folder> --data <folder>
        [--only name ...] [--skip-network]

`--data` must contain the files not stored in Git: mmc2.xlsx (Fontenele et al.
2021 Supplementary Table 2), and sequences.fasta / sequences.csv (the ten-record
NCBI Virus H5N1 example). Missing files skip their scenario; nothing is faked.

The script never sends credentials and never enables paid AI. It sends one
request at a time, except the deliberate two-request busy-server check.
"""

import argparse
import hashlib
import html
import io
import json
import re
import sys
import threading
import time
from html.parser import HTMLParser
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
TIMEOUT = 300


class FormReader(HTMLParser):
    """Fields of every <form> on a page, as a browser would submit them by default."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.forms, self.current, self.select, self.textarea = [], None, None, None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "form":
            self.current = {"action": attrs.get("action", ""), "fields": [], "selects": {}}
            self.forms.append(self.current)
        elif self.current is None:
            return
        elif tag == "input" and attrs.get("name"):
            kind = attrs.get("type", "text").lower()
            if kind in {"checkbox", "radio"} and "checked" not in attrs:
                return
            if kind in {"file", "submit", "button"}:
                return
            self.current["fields"].append((attrs["name"], attrs.get("value", "on" if kind == "checkbox" else "")))
        elif tag == "select" and attrs.get("name"):
            self.select = attrs["name"]
            self.current["selects"][self.select] = None
        elif tag == "option" and self.select:
            if self.current["selects"][self.select] is None or "selected" in attrs:
                if "selected" in attrs or self.current["selects"][self.select] is None:
                    self.current["selects"][self.select] = (attrs.get("value", ""), "selected" in attrs)
        elif tag == "textarea" and attrs.get("name"):
            self.textarea = attrs["name"]
            self.current["fields"].append((self.textarea, ""))

    def handle_endtag(self, tag):
        if tag == "select":
            self.select = None
        elif tag == "textarea":
            self.textarea = None
        elif tag == "form":
            self.current = None

    def handle_data(self, data):
        if self.textarea and self.current:
            name, value = self.current["fields"][-1]
            self.current["fields"][-1] = (name, value + data)


def form_fields(page, action="/review"):
    reader = FormReader()
    reader.feed(page)
    for form in reader.forms:
        if form["action"] == action:
            fields = dict(form["fields"])
            for name, chosen in form["selects"].items():
                if chosen is not None:
                    fields[name] = chosen[0]
            return fields
    raise ValueError(f"No form with action {action}")


def words(page):
    text = re.sub(r"<[^>]+>", " ", re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", page, flags=re.S))
    return re.sub(r"\s+", " ", html.unescape(text))


class Site:
    def __init__(self, base, output):
        self.base = base.rstrip("/")
        self.output = Path(output)
        self.output.mkdir(parents=True, exist_ok=True)
        self.http = requests.Session()
        self.http.headers["User-Agent"] = "hosted-validation-script"
        self.record = {}

    def get(self, path, **kwargs):
        return self.http.get(self.base + path, timeout=TIMEOUT, **kwargs)

    def post(self, path, **kwargs):
        return self.http.post(self.base + path, timeout=TIMEOUT, **kwargs)

    def upload(self, paths, extra=None):
        files = [("files", (Path(p).name, Path(p).read_bytes())) for p in paths]
        start = time.perf_counter()
        response = self.post("/review", files=files, data=extra or {})
        return response, time.perf_counter() - start

    def save(self, folder, name, content):
        path = self.output / folder / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content if isinstance(content, bytes) else content.encode("utf-8"))
        return path

    def download_all(self, folder, page):
        links = dict((name, link) for link, name in re.findall(r'href="(/results/[^"]+/([^"/]+))"', page))
        files = {}
        for name, link in sorted(links.items()):
            start = time.perf_counter()
            response = self.get(link)
            seconds = time.perf_counter() - start
            self.save(folder + "/downloads", name, response.content)
            files[name] = {"status": response.status_code, "bytes": len(response.content), "seconds": round(seconds, 1),
                           "md5": hashlib.md5(response.content).hexdigest()}
        return files

    def scenario(self, name, paths, choices=None, extra=None, expect_ready=True, download=True):
        """Upload -> review -> (choices) -> run -> downloads; returns the run page text."""

        entry = {"files": [Path(p).name for p in paths]}
        self.record[name] = entry
        response, seconds = self.upload(paths, extra)
        entry.update(review_status=response.status_code, review_seconds=round(seconds, 1))
        page = response.text
        self.save(name, "review.html", page)
        token = re.search(r'name="token" value="([^"]+)"', page)
        if not token:
            entry["review_error"] = words(page)[:600]
            return None
        fields = form_fields(page)
        fields.update(choices or {})
        if choices:
            start = time.perf_counter()
            page = self.post("/review", data=fields).text
            entry["update_seconds"] = round(time.perf_counter() - start, 1)
            self.save(name, "review_updated.html", page)
            fields = form_fields(page)
            fields.update(choices)
        text = words(page)
        entry["ready"] = "disabled>Run analysis" not in page
        entry["review_problems"] = re.findall(r"Required information to fix before running(.{0,600})", text)[:1]
        entry["readiness_badges"] = re.findall(r"(\d+ (?:available|partially available|unavailable|experimental))",
                                               text)
        profile = self.post("/review/profile", data=fields)
        if profile.status_code == 200:
            self.save(name, "mapping_profile.json", profile.content)
            entry["mapping_profile_bytes"] = len(profile.content)
        if not expect_ready or not entry["ready"]:
            return text
        start = time.perf_counter()
        run = self.post("/run", data=fields)
        entry.update(run_status=run.status_code, run_seconds=round(time.perf_counter() - start, 1))
        self.save(name, "results.html", run.text)
        entry["results_html_bytes"] = len(run.content)
        text = words(run.text)
        entry["errors"] = re.findall(r"(The analysis could not be completed[^.]*\.[^.]*\.|could not be completed\.[^.]*\.)",
                                     text)[:3]
        major = re.search(r"Major findings(.*?)(?:Files, connection details|Important limitations)", text)
        entry["major_findings"] = major.group(1).strip()[:1500] if major else None
        if download:
            entry["downloads"] = self.download_all(name, run.text)
        return text


def check_basics(site):
    for path in ["/health", "/"]:
        start = time.perf_counter()
        response = site.get(path)
        site.record[f"GET {path}"] = {"status": response.status_code, "seconds": round(time.perf_counter() - start, 2),
                                      "bytes": len(response.content),
                                      "text": response.text[:20] if path == "/health" else None}
    response = site.post("/system-status")
    site.save("system_status", "page.html", response.text)
    text = words(response.text)
    status = re.search(r"Service What it does here Status(.{0,1500})", text)
    site.record["system-status"] = {"status": response.status_code, "table": status.group(1) if status else None}


def run_scenarios(site, data, only, skip_network):
    integrated = EXAMPLES / "integrated"
    annotation = EXAMPLES / "annotation"
    wanted = (lambda name: not only or name in only)
    if wanted("sprint2_synthetic"):
        site.scenario("sprint2_synthetic", [EXAMPLES / "synthetic_sequences.fasta", EXAMPLES / "synthetic_metadata.csv"])
    ncbi = [data / "sequences.fasta", data / "sequences.csv"]
    if wanted("sprint2_ncbi10") and all(p.exists() for p in ncbi):
        site.scenario("sprint2_ncbi10", ncbi)
    if wanted("integrated_example"):
        site.scenario("integrated_example", [integrated / name for name in [
            "synthetic_individual_sequences.fasta", "synthetic_individual_metadata.csv", "synthetic_reference.fasta",
            "synthetic_pooled_snvs.csv", "synthetic_pooled_base_counts.csv", "synthetic_variant_signatures.csv"]])
    mmc2 = data / "mmc2.xlsx"
    if wanted("fontenele_mmc2") and mmc2.exists():
        site.scenario("fontenele_mmc2", [mmc2, annotation / "NC_045512.2.fasta", annotation / "NC_045512.2.gb"],
                      choices={"role_0": "snv", "role_1": "reference", "role_2": "annotation"})
    if wanted("fontenele_only") and mmc2.exists():  # missing optional files: no FASTA, no annotation
        site.scenario("fontenele_only", [mmc2], download=False)
    if wanted("formula_safety"):
        table = ("Sample_ID,Reference_Genome,Position,Ref,Alt,Alt_Frequency,Read_Depth,Alt_Depth,Collection_Date\n"
                 '"=HYPERLINK(""https://example.invalid/"",""x"")",SYN.1,2,C,T,0.1,100,10,2026-01-05\n'
                 "S2,SYN.1,7,G,A,0.25,120,30,2026-01-12\nS3,SYN.1,7,G,A,0.25,160,40,2026-01-19\n")
        path = site.save("formula_safety", "formula_snvs.csv", table)
        reference = site.save("formula_safety", "syn_reference.fasta", ">SYN.1 synthetic\n" + "ACGTACGTAC" * 12 + "\n")
        site.scenario("formula_safety", [path, reference], choices={"role_0": "snv", "role_1": "reference"})
        check_formulas(site)
    if wanted("unsupported"):
        check_unsupported(site)
    if wanted("busy") and mmc2.exists():
        check_busy(site, mmc2)
    if wanted("literature") and not skip_network and all(p.exists() for p in ncbi):
        check_literature(site, ncbi)


def check_formulas(site):
    import openpyxl

    folder = site.output / "formula_safety" / "downloads"
    result = {"xlsx_formula_cells": None, "csv_cells_starting_with_formula": []}
    if (folder / "analysis_results.xlsx").exists():
        book = openpyxl.load_workbook(folder / "analysis_results.xlsx")
        result["xlsx_formula_cells"] = sum(1 for sheet in book for row in sheet.iter_rows() for cell in row
                                           if cell.data_type == "f")
    for path in folder.glob("*.csv"):
        text = path.read_text(encoding="utf-8")
        if re.search(r'(^|,)"?=HYPERLINK', text, flags=re.M):
            result["csv_cells_starting_with_formula"].append(path.name)
    site.record["formula_safety"]["check"] = result


def check_unsupported(site):
    record = site.record.setdefault("unsupported", {})
    fastq = site.save("unsupported", "reads.fastq", "@r1\nACGT\n+\nIIII\n")
    response, _ = site.upload([fastq])
    record["fastq"] = {"status": response.status_code, "message": re.findall(r"Raw sequencing reads[^.]*\.",
                                                                              words(response.text))[:1]}
    rows = 160_000  # 1.6 million cells (> MAX_WORKSHEET_CELLS 1.5 M on the server), about 3.5 MB
    big = "a,b,c,d,e,f,g,h,i,j\n" + "1,2,3,4,5,6,7,8,9,0\n" * rows
    path = site.save("unsupported", "too_many_cells.csv", big)
    response, seconds = site.upload([path])
    record["too_many_cells"] = {"status": response.status_code, "seconds": round(seconds, 1),
                                "message": re.findall(r"larger than this server can process safely[^.]*\.",
                                                      words(response.text))[:1]}
    path = site.save("unsupported", "too_large.csv", "x" * (11 * 1024 * 1024))
    response, _ = site.upload([path])
    record["over_10_mb"] = {"status": response.status_code, "message": words(response.text)[:200]}
    response = site.get("/health")
    record["health_after"] = response.status_code


def check_busy(site, mmc2):
    """Two runs at once: the second must be told the server is busy and keep its files."""

    record = site.record.setdefault("busy", {})
    tokens = []
    for _ in range(2):
        response, _ = site.upload([mmc2])
        page = response.text
        fields = form_fields(page)
        tokens.append(fields)
    results = {}

    def run(index):
        start = time.perf_counter()
        response = requests.post(site.base + "/run", data=tokens[index], timeout=TIMEOUT)
        results[index] = (response.status_code, time.perf_counter() - start, words(response.text))

    first = threading.Thread(target=run, args=(0,))
    first.start()
    time.sleep(3)
    run(1)
    first.join()
    for index, (status, seconds, text) in results.items():
        record[f"run_{index}"] = {"status": status, "seconds": round(seconds, 1),
                                  "busy_message": bool(re.search(r"Another analysis is running", text)),
                                  "results": "Major findings" in text}


def check_literature(site, ncbi):
    """Free literature discovery for the ten-record H5N1 example (no Claude request)."""

    record = site.record.setdefault("literature", {})
    response, _ = site.upload(ncbi, extra={"enrich": "on", "enrich_shown": "1"})
    page = response.text
    fields = form_fields(page)
    start = time.perf_counter()
    response = site.post("/enrich/start", data=fields)
    record["status"] = response.status_code
    record["seconds"] = round(time.perf_counter() - start, 1)
    site.save("literature", "findings.html", response.text)
    text = words(response.text)
    record["text"] = text[:6000]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--data", required=True)
    parser.add_argument("--only", nargs="*")
    parser.add_argument("--skip-network", action="store_true")
    args = parser.parse_args()
    site = Site(args.base, args.output)
    start = time.perf_counter()
    try:
        check_basics(site)
        run_scenarios(site, Path(args.data), args.only, args.skip_network)
    finally:
        site.record["total_seconds"] = round(time.perf_counter() - start, 1)
        (site.output / "record.json").write_text(json.dumps(site.record, indent=1, default=str), encoding="utf-8")
    print(json.dumps({key: {k: v for k, v in value.items() if k not in {"downloads", "text", "major_findings"}}
                      if isinstance(value, dict) else value for key, value in site.record.items()},
                     indent=1, default=str)[:6000])


if __name__ == "__main__":
    sys.exit(main())
