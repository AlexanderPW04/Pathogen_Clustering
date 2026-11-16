"""Hosted checks for the final-polish changes (genetic code, gene "NA", profiles, downloads, isolation).

Usage:
    python validation/hosted_polish_checks.py --base <url> --output <folder> --h5n1 <folder> [--profile <json>]

`--h5n1` must hold `variants/<RUN>_NA_variants.tsv`, `ref_NA.fasta`, and
`h5n1_all_segments_relabelled.gff3` (real andersen-lab H5N1 iVar tables with
segment-labelled NCBI annotation; see docs/RENDER_HOSTED_VALIDATION_REPORT.md).
`--profile` is a mapping profile saved by the earlier, validated version.
Everything else is SYNTHETIC. One request at a time; no credentials; no paid AI.
"""

import argparse
import io
import json
import re
import sys
import threading
import time
from pathlib import Path

import openpyxl
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hosted_smoke_test import Site, form_fields, words  # noqa: E402

# A 16-base SYNTHETIC genome with one CDS at 3..14: ATG TGG AAA TAA. Position 8 G>A turns TGG into TGA:
# a stop under the standard code (table 1), tryptophan under table 4 (W2* versus W2W; NCBI tables).
GENOME = "CC" + "ATGTGGAAATAA" + "GG"
CODE_FASTA = f">CODE.1 synthetic genetic-code test\n{GENOME}\n"
CODE_GFF3 = ("##gff-version 3\n##sequence-region CODE.1 1 16\n"
             "CODE.1\ttest\tCDS\t3\t14\t.\t+\t0\tID=cds-C1;gene=codeA;protein_id=C1\n")
CODE_SNVS = ("Sample_ID,Reference_Genome,Position,Ref,Alt,Alt_Frequency,Read_Depth,Alt_Depth,Collection_Date\n"
             "S1,CODE.1,8,G,A,0.2,100,20,2026-01-05\nS2,CODE.1,8,G,A,0.3,100,30,2026-01-12\n"
             "S3,CODE.1,8,G,A,0.4,100,40,2026-01-19\n")


def run_upload(site, name, paths, choices):
    text = site.scenario(name, paths, choices=choices)
    return text, site.record[name]


def table(site, name, filename):
    path = site.output / name / "downloads" / filename
    return pd.read_csv(path, keep_default_na=False, na_values=[""]) if path.exists() else None


def genetic_code(site):
    folder = site.output / "_inputs"
    folder.mkdir(parents=True, exist_ok=True)
    paths = [site.save("_inputs", "code_snvs.csv", CODE_SNVS), site.save("_inputs", "code.fasta", CODE_FASTA),
             site.save("_inputs", "code.gff3", CODE_GFF3)]
    roles = {"role_0": "snv", "role_1": "reference", "role_2": "annotation"}
    outcome = {}
    for label, code in [("not_confirmed", ""), ("confirmed_4", "4"), ("confirmed_1", "1")]:
        name = f"genetic_code_{label}"
        page, entry = run_upload(site, name, paths, {**roles, "genetic_code": code})
        mutations = table(site, name, "sprint5_annotated_mutations.csv")
        sources = table(site, name, "sprint5_annotation_sources.csv")
        row = mutations.iloc[0].to_dict() if mutations is not None and len(mutations) else {}
        outcome[label] = {
            "run_status": entry.get("run_status"),
            "review_warning": bool(re.search(r"state no genetic code", words(
                (site.output / name / "review.html").read_text(encoding="utf-8")))),
            "consequence": row.get("Consequence_Category"), "amino_acid_change": row.get("Amino_Acid_Changes"),
            "annotation_status": row.get("Annotation_Status"),
            "table_source": sources["Translation_Table_Source"].tolist() if sources is not None else None,
            "cds_status": sources["Status"].tolist() if sources is not None else None}
    site.record["genetic_code"] = outcome


def gene_na(site, h5n1):
    # A table with only insertions/deletions has no single-base variants and is refused at review
    # ("The processed SNV table has no records"); only tables with at least one SNV are uploaded.
    variants = [path for path in sorted((h5n1 / "variants").glob("*_NA_variants.tsv"))
                if (~pd.read_csv(path, sep="\t", dtype=str)["ALT"].str[0].isin(["+", "-"])).any()]
    site.record["influenza_na_skipped_indel_only"] = sorted(
        path.name for path in (h5n1 / "variants").glob("*_NA_variants.tsv") if path not in variants)
    paths = variants + [h5n1 / "ref_NA.fasta", h5n1 / "h5n1_all_segments_relabelled.gff3"]
    roles = {f"role_{index}": "snv" for index in range(len(variants))}
    roles.update({f"role_{len(variants)}": "reference", f"role_{len(variants) + 1}": "annotation",
                  "genetic_code": "1"})
    page, entry = run_upload(site, "influenza_na", paths, roles)
    mutations = table(site, "influenza_na", "sprint5_annotated_mutations.csv")
    record = {"run_status": entry.get("run_status"), "files": len(paths)}
    if mutations is not None:
        genes = mutations["Genes"]
        record.update(rows=len(mutations), csv_genes_na=int((genes == "NA").sum()),
                      csv_genes_missing=int(genes.isna().sum()),
                      examples=list(mutations.loc[genes == "NA", "Amino_Acid_Changes"].head(3)))
    workbook = site.output / "influenza_na" / "downloads" / "analysis_results.xlsx"
    if workbook.exists():
        book = openpyxl.load_workbook(workbook, read_only=True)
        rows = list(book["Annotated mutations"].iter_rows(values_only=True))
        column = rows[0].index("Genes")
        cells = [row[column] for row in rows[1:]]
        record.update(xlsx_genes_na=sum(cell == "NA" for cell in cells), xlsx_genes_blank=sum(cell is None for cell in cells))
    html = (site.output / "influenza_na" / "results.html").read_text(encoding="utf-8")
    record["website_mentions_NA_change"] = bool(re.search(r"\bNA:[A-Z*]\d+[A-Z*]", words(html)))
    site.record["influenza_na"] = record


def profiles(site, old_profile):
    record = {}
    csv = ("Sample_ID,Reference_Genome,Position,Ref,Alt,Alt_Frequency,Read_Depth\n"
           "S1,SYN.1,2,C,T,0.1,100\nS2,SYN.1,7,G,A,0.2,100\n")
    response, _ = site.upload([site.save("_inputs", "profile_test.csv", csv)])
    token = form_fields(response.text)["token"]
    bad = {"profile_format": "genomic-surveillance-mapping-profile", "profile_version": 1,
           "files": [{"role": "snv", "fields": {}, "options": {"header_row": -5}}]}
    page = site.post("/review", data={"token": token},
                     files={"mapping_profile": ("bad.json", json.dumps(bad).encode())}).text
    record["negative_header_row"] = re.findall(r"'header_row' is -5; it must be [^.]*\.", words(page))[:1]
    page = site.post("/review", data={"token": token},
                     files={"mapping_profile": ("bad.json", b"{not json")}).text
    record["invalid_json"] = re.findall(r"The mapping profile is not valid JSON\.", words(page))[:1]
    page = site.post("/review", data={"token": token},
                     files={"mapping_profile": ("evil.json", json.dumps({**bad, "files": [
                         {"role": "../../app", "fields": {}, "options": {}}]}).encode())}).text
    record["unknown_role"] = re.findall(r"unknown role '[^']*'", words(page))[:1]
    record["review_kept_after_bad_profile"] = token in page
    if old_profile and old_profile.exists():
        mmc2 = Path(old_profile).parent / "mmc2.xlsx"
        if mmc2.exists():
            response, _ = site.upload([mmc2])
            token = form_fields(response.text)["token"]
            page = site.post("/review", data={"token": token},
                             files={"mapping_profile": ("old.json", old_profile.read_bytes())}).text
            record["old_profile_report"] = re.findall(r"Profile entry 1 \([a-z_]+\) [^.]*\.", words(page))[:2]
    site.record["profiles"] = record


def downloads_and_isolation(site):
    record = {}
    page = (site.output / "genetic_code_confirmed_4" / "results.html").read_text(encoding="utf-8")
    link = re.search(r'href="(/results/[^"]+/analysis_results\.xlsx)"', page).group(1)
    results = {}

    def click(index):
        start = time.perf_counter()
        response = requests.get(site.base + link, timeout=300)
        results[index] = (response.status_code, len(response.content), round(time.perf_counter() - start, 1),
                          response.content[:2])

    threads = [threading.Thread(target=click, args=(index,)) for index in range(3)]
    for thread in threads:
        thread.start()
        time.sleep(0.2)
    for thread in threads:
        thread.join()
    record["three_concurrent_xlsx_clicks"] = results
    pdf = re.search(r'href="(/results/[^"]+/analysis_report\.pdf)"', page).group(1)
    first = site.get(pdf)
    again = site.get(pdf)
    record["pdf_twice"] = [first.status_code, again.status_code, first.content == again.content,
                           first.content[:5].decode("latin-1")]
    download_id = link.split("/")[2]
    record["other_token"] = site.get(f"/results/{'x' * len(download_id)}/analysis_results.xlsx").status_code
    record["traversal"] = [site.get(path).status_code for path in
                           [f"/results/{download_id}/..%2f..%2fapp.py", "/examples/..%2fapp.py",
                            "/examples/%2e%2e%2fDockerfile"]]
    record["unknown_file_in_valid_results"] = site.get(f"/results/{download_id}/not_a_result.csv").status_code
    garbage, _ = site.upload([site.save("_inputs", "garbage.csv", "\x00\x01\x02 not a table \xff")])
    garbage_text = words(garbage.text)
    record["garbage_csv"] = {"status": garbage.status_code,
                             "run_blocked": "disabled>Run analysis" in garbage.text,
                             "shown_as": re.findall(r"garbage\.csv (.{0,45})", garbage_text)[:1]}
    broken = site.save("_inputs", "broken.xlsx", b"PK\x03\x04 this is not a real workbook")
    response, _ = site.upload([broken])
    record["broken_xlsx"] = {"status": response.status_code,
                             "message": re.findall(r"Could not read this file[^.]*\.", words(response.text))[:1]}
    record["health_after"] = site.get("/health").status_code
    site.record["downloads_isolation"] = record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--h5n1", required=True)
    parser.add_argument("--profile")
    args = parser.parse_args()
    site = Site(args.base, args.output)
    start = time.perf_counter()
    try:
        genetic_code(site)
        gene_na(site, Path(args.h5n1))
        profiles(site, Path(args.profile) if args.profile else None)
        downloads_and_isolation(site)
    finally:
        site.record["total_seconds"] = round(time.perf_counter() - start, 1)
        (site.output / "record.json").write_text(json.dumps(site.record, indent=1, default=str), encoding="utf-8")
    keys = ["genetic_code", "influenza_na", "profiles", "downloads_isolation", "total_seconds"]
    print(json.dumps({key: site.record.get(key) for key in keys}, indent=1, default=str))


if __name__ == "__main__":
    sys.exit(main())
