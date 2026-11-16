"""Local UI preview of the homepage and every Sprint 4 state, with SYNTHETIC fixtures only.

Run it on your own computer (never deployed):

    python validation/local_ui_preview.py            # then open http://127.0.0.1:5058/preview

It serves the real application with one extra page, /preview. Each Sprint 4 link stages the
synthetic NCBI-shaped influenza test records (tests/fixtures/ncbi_virus_structure_*.csv/.fasta,
ten invented accessions collected in September 2026) and opens the real Sprint 4 page. The only
thing replaced is the publication source: an in-memory set of invented articles from the test
suite answers instead of NCBI and PubTator. Claude is never called: the session is marked
synthetic, so the offline mock provider (fixed answers) is used, and the pages carry a
"LOCAL UI PREVIEW" banner and the application's own SYNTHETIC labels.

Nothing is sent to NCBI, PubTator, GeoNames, or Anthropic, and nothing is billed. "Continue to
Genomic Analysis" runs the real analysis on the synthetic records (MAFFT must be on PATH).
"""

import json
import os
import sys
import tempfile
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
    if response.mimetype == "text/html" and response.status_code == 200:
        response.set_data(response.get_data(as_text=True).replace("<body>", "<body>" + BANNER, 1))
    return response


@app.route("/preview")
def preview_index():
    items = "".join(f'<li><a href="/preview/sprint4/{key}">{escape(label)}</a></li>' for key, label in STATES)
    modes = "".join(f'<li><a href="/preview/home/{key}">{label}</a></li>' for key, label in
                    [("individual", "Individual Sequences"), ("pooled", "Pooled Sequencing"),
                     ("combined", "Combined Analysis"), ("auto", "Detect from my files (the default)")])
    return (f"<!doctype html><meta charset='utf-8'><title>UI preview</title><body style='font-family:sans-serif;"
            f"max-width:52rem;margin:1rem auto;padding:0 1rem'><h1>Local UI preview (synthetic)</h1>"
            f"<p>Every Sprint 4 state uses ten invented NCBI-shaped influenza records and invented articles. "
            f"Nothing is sent to NCBI, PubTator, or Claude.</p><h2>Homepage</h2><ul>{modes}</ul>"
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
