"""Local UI preview of the homepage and every Sprint 4 state, with SYNTHETIC fixtures only.

Run it on your own computer (never deployed):

    python validation/local_ui_preview.py            # then open http://127.0.0.1:5058/preview

It serves the real application with one extra page, /preview. The results links run the real
analysis on the bundled synthetic example files (individual only, pooled only, combined, and
missing-data variants made from them in memory); the optional Fontenele links run the published
mmc2.xlsx workbook only when FONTENELE_MMC2 points to a local copy (it is never stored in Git).
Each Sprint 4 link stages the
synthetic NCBI-shaped influenza test records (tests/fixtures/ncbi_virus_structure_*.csv/.fasta,
ten invented accessions collected in September 2026) and opens the real Sprint 4 page. The only
thing replaced is the publication source: an in-memory set of invented articles from the test
suite answers instead of NCBI and PubTator. Claude is never called: the session is marked
synthetic, so the offline mock provider (fixed answers) is used, and the pages carry a
"LOCAL UI PREVIEW" banner and the application's own SYNTHETIC labels.

Nothing is sent to NCBI, PubTator, GeoNames, or Anthropic, and nothing is billed. "Continue to
Genomic Analysis" runs the real analysis on the synthetic records (MAFFT must be on PATH).
"""

import io
import json
import os
import re
import sys
import tempfile
import time
from contextlib import redirect_stdout
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.update(AI_ENRICHMENT_ENABLED="false", ANTHROPIC_API_KEY="", AI_ACCESS_CODE="", GEONAMES_USERNAME="")
sys.path[:0] = [str(ROOT), str(ROOT / "tests")]

from werkzeug.datastructures import FileStorage  # noqa: E402

import app as web  # noqa: E402
import literature_view  # noqa: E402
from literature.http_client import ServiceUnavailable  # noqa: E402
from test_literature_eligibility import FakeSource, old_study_article, passage  # noqa: E402

FIXTURES = ROOT / "tests" / "fixtures"
INTEGRATED = ROOT / "examples" / "integrated"
ANNOTATION = ROOT / "examples" / "annotation"
FILES = ["ncbi_virus_structure_sequences.fasta", "ncbi_virus_structure_metadata.csv"]
ACCESSIONS = [f"TST{number:06d}" for number in range(1, 11)]
DIRECT_SENTENCE = "Sample TST000001 was collected at a farm in Kingsbury County, South Dakota."


def invented_article(number, sentence, date="2026 Dec"):
    article_id = f"PMID:{number}"
    return article_id, {"Title": f"Invented article {number} (synthetic preview)", "Publication_Date": date,
                        "Passages": [passage(article_id, 1, "METHODS", sentence)]}


class UnreachableText(FakeSource):
    """The article is found, but its text service answers with an error."""

    def fetch(self, article_ids):
        raise ServiceUnavailable("PubTator3 returned HTTP 503 after 3 attempt(s) (synthetic preview).")


def study_links(article_id):
    return {accession: [(article_id, "study")] for accession in ACCESSIONS}


def sources():
    general = invented_article(9400002, "Samples for this surveillance program came from Maryland and Arizona.")
    direct = invented_article(9400008, DIRECT_SENTENCE)
    return {
        "one-publication": FakeSource({"PMID:9000001": old_study_article()}, study_links("PMID:9000001")),
        "general-mention": FakeSource(dict([general]), study_links(general[0])),
        "no-publication": FakeSource({}, {}),
        "eligible": FakeSource(dict([direct]), {"TST000001": [(direct[0], "direct")]}),
        "retrieval-error": UnreachableText(dict([direct]), {"TST000001": [(direct[0], "direct")]}),
    }


STATES = [
    ("before", "Before searching: the file review with the Sprint 4 step and the Find and Check Publications button"),
    ("one-publication", "One publication found (like the real H5N1 case): it predates the samples; nothing changes"),
    ("general-mention", "A publication names places, but not for these samples: no usable metadata improvement"),
    ("no-publication", "No publication found"),
    ("eligible", "Sample-specific evidence found: the Claude step is offered (offline mock; open it to see "
                 "the cost approval, then the proposals that need your approval)"),
    ("retrieval-error", "The article was found but its text could not be retrieved: an incomplete review"),
]
MOCK_FACTS = {"extraction": {"PMID:9400008": {"facts": [{
    "Sample_ID": "TST000001.1", "Field": "Location", "Value": "Kingsbury County, South Dakota, USA",
    "Location_Role": "Sample collection site", "Geographic_Level": "County", "Passage_ID": "PMID:9400008#p1",
    "Evidence_Quote": DIRECT_SENTENCE, "Evidence_Type": "Directly stated", "Ambiguity": "",
    "Reasoning_Summary": "The passage names the sample and its collection site."}],
    "insufficient_evidence_samples": [], "notes": ""}}}
MOCK_FILE = Path(tempfile.gettempdir()) / "genomic_surveillance_ui_preview_mock.json"
PREVIEW_TOKENS = {}  # staged upload token -> preview state
real_source = literature_view.literature_source


def preview_source(session):
    """The invented articles of this preview state; any other session keeps its usual source."""

    state = PREVIEW_TOKENS.get(session.get("staged_token"))
    if state is None:
        return real_source(session)
    session["mock"] = str(MOCK_FILE)
    return sources()[state]


literature_view.literature_source = preview_source
app = web.app
BANNER = ('<div style="background:#5b2a86;color:#fff;font-weight:bold;padding:0.4rem 1rem;text-align:center">'
          'LOCAL UI PREVIEW · synthetic fixtures · not real NCBI, PubTator, or Claude results · '
          '<a style="color:#fff" href="/preview">all preview states</a></div>')


@app.after_request
def add_banner(response):
    if response.mimetype == "text/html" and response.status_code == 200 and "LOCAL UI PREVIEW" not in             response.get_data(as_text=True):
        response.set_data(response.get_data(as_text=True).replace("<body>", "<body>" + BANNER, 1))
    return response


def example_bytes(name, drop=()):
    """A bundled synthetic example file, optionally without some metadata columns (missing-data states)."""

    path = INTEGRATED / name
    if not drop:
        return path.read_bytes()
    import pandas as pd
    return pd.read_csv(path).drop(columns=list(drop)).to_csv(index=False).encode()


def named_upload(fields):
    """Run the real analysis through the named-field upload (POST /) and return the results page."""

    data = {key: (io.BytesIO(content), name) for key, (name, content) in fields.items()}
    with app.test_client() as client, redirect_stdout(io.StringIO()):
        return client.post("/", data=data, content_type="multipart/form-data").get_data(as_text=True)


def guided_upload(files, mode):
    """Run the real guided flow (/review, then /run) and return the results page."""

    uploads = [(io.BytesIO(content), name) for name, content in files]
    with app.test_client() as client, redirect_stdout(io.StringIO()):
        review = client.post("/review", data={"files": uploads, "intended_mode": mode},
                             content_type="multipart/form-data").get_data(as_text=True)
        token = re.search(r'name="token" value="([^"]+)"', review)
        if token is None:
            return review
        return client.post("/run", data={"token": token.group(1), "intended_mode": mode}).get_data(as_text=True)


def example_page(button):
    with app.test_client() as client, redirect_stdout(io.StringIO()):
        return client.post("/", data={button: "1"}).get_data(as_text=True)


def fontenele_page(with_reference):
    workbook = Path(os.environ.get("FONTENELE_MMC2", ""))
    if not workbook.is_file():
        return None
    files = [(workbook.name, workbook.read_bytes())]
    if with_reference:
        files.append(("NC_045512.2.fasta", (ANNOTATION / "NC_045512.2.fasta").read_bytes()))
    files.append(("NC_045512.2.gb", (ANNOTATION / "NC_045512.2.gb").read_bytes()))
    return guided_upload(files, "pooled")


INDIVIDUAL = {"fasta_file": ("synthetic_individual_sequences.fasta", "synthetic_individual_sequences.fasta"),
              "metadata_file": ("synthetic_individual_metadata.csv", "synthetic_individual_metadata.csv")}


def individual_fields(reference=False, drop=()):
    fields = {key: (name, example_bytes(source, drop if key == "metadata_file" else ()))
              for key, (name, source) in INDIVIDUAL.items()}
    if reference:
        fields["reference_file"] = ("synthetic_reference.fasta", example_bytes("synthetic_reference.fasta"))
    return fields


RESULT_STATES = {
    "individual": ("Individual sequences only (Sprint 2): synthetic FASTA + metadata, no reference",
                   lambda: named_upload(individual_fields())),
    "individual-reference": ("Individual sequences with a reference genome (Sprint 2 + reference positions)",
                             lambda: named_upload(individual_fields(reference=True))),
    "pooled": ("Pooled sequencing only (Sprint 3 + 5): synthetic SNV table, no reference FASTA",
               lambda: named_upload({"snv_file": ("synthetic_pooled_snvs.csv",
                                                  example_bytes("synthetic_pooled_snvs.csv"))})),
    "combined": ("Combined (Sprint 2 + 3 + 5): the bundled synthetic integrated example",
                 lambda: example_page("run_integrated_example")),
    "pooled-no-dates": ("Missing data: pooled SNV table without collection dates or locations",
                        lambda: named_upload({"snv_file": ("synthetic_pooled_snvs_no_dates.csv", example_bytes(
                            "synthetic_pooled_snvs.csv", drop=("Collection_Date", "Location_Name")))})),
    "individual-sparse": ("Missing data: individual sequences whose metadata has no coordinates",
                          lambda: named_upload(individual_fields(drop=("Latitude", "Longitude")))),
    "sprint4": ("With Sprint 4: the synthetic literature-enrichment example (offline mock, no Claude call)",
                lambda: example_page("run_literature_example")),
    "fontenele": ("Real data (local only): Fontenele et al. 2021 mmc2.xlsx + NC_045512.2 FASTA + GenBank "
                  "(needs FONTENELE_MMC2)", lambda: fontenele_page(True)),
    "fontenele-no-fasta": ("Real data (local only): mmc2.xlsx + GenBank, no reference FASTA (needs FONTENELE_MMC2)",
                           lambda: fontenele_page(False)),
}
RESULT_CACHE = {}  # state -> (time, page); downloads live 10 minutes on the server, so pages are reused for 8


@app.route("/preview/results/<state>")
def preview_results(state):
    if state not in RESULT_STATES:
        return "Unknown preview state.", 404
    cached = RESULT_CACHE.get(state)
    if cached is None or time.time() - cached[0] > 480:
        page = RESULT_STATES[state][1]()
        if page is None:
            return ("Set FONTENELE_MMC2 to a local copy of the published mmc2.xlsx to preview the real "
                    "Fontenele et al. (2021) result. It is never stored in Git."), 404
        RESULT_CACHE[state] = (time.time(), page)
    return RESULT_CACHE[state][1]


@app.route("/preview")
def preview_index():
    items = "".join(f'<li><a href="/preview/sprint4/{key}">{escape(label)}</a></li>' for key, label in STATES)
    results = "".join(f'<li><a href="/preview/results/{key}">{escape(label)}</a></li>'
                      for key, (label, _) in RESULT_STATES.items())
    modes = "".join(f'<li><a href="/preview/home/{key}">{label}</a></li>' for key, label in
                    [("individual", "Individual Sequences"), ("pooled", "Pooled Sequencing"),
                     ("combined", "Combined Analysis"), ("auto", "Detect from my files (the default)")])
    return (f"<!doctype html><meta charset='utf-8'><title>UI preview</title><body style='font-family:sans-serif;"
            f"max-width:52rem;margin:1rem auto;padding:0 1rem'><h1>Local UI preview (synthetic)</h1>"
            f"<p>Every Sprint 4 state uses ten invented NCBI-shaped influenza records and invented articles. "
            f"Nothing is sent to NCBI, PubTator, or Claude.</p><h2>Homepage</h2><ul>{modes}</ul>"
            f"<h2>Results dashboard</h2><p>Each link runs the real analysis (the Fontenele links take about a "
            f"minute). Results use synthetic example files unless marked real data.</p><ul>{results}</ul>"
            f"<h2>Sprint 4 states</h2><ul>{items}</ul></body>")


@app.route("/preview/home/<mode>")
def preview_home(mode):
    """The homepage with one analysis card already chosen (through the real Back-to-upload route)."""

    if mode not in web.MODES:
        return "Unknown analysis mode.", 404
    return (f"<!doctype html><meta charset='utf-8'><body><form method='post' action='/back'>"
            f"<input type='hidden' name='target' value='upload'>"
            f"<input type='hidden' name='intended_mode' value='{mode}'>"
            f"<noscript><button>Open the preview</button></noscript></form>"
            f"<script>document.forms[0].submit()</script></body>")


@app.route("/preview/sprint4/<state>")
def preview_state(state):
    if state not in dict(STATES):
        return "Unknown preview state.", 404
    handles = [(FIXTURES / name).open("rb") for name in FILES]
    try:
        token = web.stage_uploads([(FileStorage(stream=handle, filename=name), None)
                                   for handle, name in zip(handles, FILES)], enrich=True, synthetic=True)
    finally:
        for handle in handles:
            handle.close()
    PREVIEW_TOKENS[token] = "one-publication" if state == "before" else state
    action = "/review" if state == "before" else "/enrich/start"
    # A self-submitting form: the real routes run exactly as after a visitor's own upload.
    return (f"<!doctype html><meta charset='utf-8'><body><form method='post' action='{action}'>"
            f"<input type='hidden' name='token' value='{escape(token)}'><input type='hidden' name='enrich' value='on'>"
            f"<input type='hidden' name='enrich_shown' value='1'>"
            f"<input type='hidden' name='intended_mode' value='individual'>"
            f"<noscript><button>Open the preview</button></noscript></form>"
            f"<script>document.forms[0].submit()</script></body>")


if __name__ == "__main__":
    MOCK_FILE.write_text(json.dumps(MOCK_FACTS), encoding="utf-8")
    app.run(host="127.0.0.1", port=int(os.getenv("PREVIEW_PORT", "5058")), threaded=True)
