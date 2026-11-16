"""Homepage + Sprint 4 UI pass (2026-11-16): what the new screens claim must match what was verified.

Mocked external services only (in-memory articles and the offline mock provider): no NCBI,
PubTator, GeoNames, or Claude request is sent, and nothing is billed. SYNTHETIC records only.
"""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_literature_eligibility import FakeSource, ten_records  # noqa: E402
from test_sprint4_interface import ACCESSION, article  # noqa: E402

import literature_view  # noqa: E402
from app import STAGED, app  # noqa: E402
from literature.providers import MockProvider  # noqa: E402
from literature.review import review_summary  # noqa: E402
from literature.workflow import find_literature, new_session, run_analysis  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
QUOTE = f"Sample {ACCESSION} was collected in Brookings, South Dakota."


def analyzed_session(synthetic=False):
    """One directly linked article, read by the offline mock provider, with one verified proposal."""

    identifier, item = article(9400011, QUOTE)
    session = new_session(individual=ten_records([ACCESSION]))
    session["synthetic"] = synthetic
    find_literature(session, FakeSource({identifier: item}, {ACCESSION: [(identifier, "direct")]}))
    with tempfile.TemporaryDirectory() as folder:
        mock = Path(folder) / "mock.json"
        mock.write_text(json.dumps({"extraction": {identifier: {"facts": [{
            "Sample_ID": f"{ACCESSION}.1", "Field": "Location", "Value": "Brookings, South Dakota, USA",
            "Location_Role": "Sample collection site", "Geographic_Level": "City or town",
            "Passage_ID": f"{identifier}#p1", "Evidence_Quote": QUOTE, "Evidence_Type": "Directly stated",
            "Ambiguity": "", "Reasoning_Summary": "Named."}], "insufficient_evidence_samples": [], "notes": ""}}}),
            encoding="utf-8")
        run_analysis(session, MockProvider(mock), {**literature_view.ai_settings(), "max_articles": 5,
                                                   "max_requests": 10, "max_input_tokens": 60000})
    return session


class PublicationAnswersTest(unittest.TestCase):
    def test_changes_applied_is_yes_only_after_a_value_was_applied(self):
        session = analyzed_session()
        card = review_summary(session)["cards"][0]
        self.assertEqual((card["claude"], card["improvement"], card["changes"]), ("Yes", "Yes", "No"))
        approvable = [fact["Fact_ID"] for fact in session["facts"] if fact["Approvable"]]
        self.assertEqual(len(approvable), 1)
        session["applied_facts"] = approvable  # what enrich_run records after applying it
        self.assertEqual(review_summary(session)["cards"][0]["changes"], "Yes")

    def test_offline_demonstration_is_never_reported_as_a_claude_analysis(self):
        card = review_summary(analyzed_session(synthetic=True))["cards"][0]
        self.assertEqual(card["claude"], "No (offline demonstration)")

    def test_reviewed_is_said_only_when_the_text_was_read(self):
        identifier = "PMID:9400012"
        empty = {"Title": "Invented article without text", "Publication_Date": "2026 Dec", "Passages": []}
        session = new_session(individual=ten_records([ACCESSION]))
        find_literature(session, FakeSource({identifier: empty}, {ACCESSION: [(identifier, "study")]}))
        summary = review_summary(session)
        self.assertFalse(summary["cards"][0]["text_checked"])
        self.assertEqual(summary["cards"][0]["location"], "Not checked")
        self.assertNotIn("reviewed", summary["headline"])
        self.assertIn("One publication was found.", summary["headline"])


class HomepageModesTest(unittest.TestCase):
    def setUp(self):
        os.environ["AI_ENRICHMENT_ENABLED"] = "false"
        self.client = app.test_client()
        self.addCleanup(STAGED.clear)

    def group(self, page, key):
        start = page.index(f'class="field-group g-{key}"')
        return page[start:page.index(">", start)]

    def test_a_chosen_mode_shows_only_its_field_groups(self):
        page = self.client.post("/back", data={"target": "upload", "intended_mode": "pooled"}).get_data(as_text=True)
        self.assertIn('value="pooled" checked', page)
        self.assertNotIn("hidden", self.group(page, "pooled"))
        self.assertNotIn("hidden", self.group(page, "shared"))
        self.assertIn("hidden", self.group(page, "individual"))
        page = self.client.post("/back", data={"target": "upload", "intended_mode": "combined"}).get_data(as_text=True)
        for key in ["individual", "pooled", "shared"]:
            self.assertNotIn("hidden", self.group(page, key))

    def test_every_upload_field_has_a_short_question_and_answer(self):
        page = self.client.get("/").get_data(as_text=True)
        self.assertEqual(page.count('<details class="field-help">'), 9)
        self.assertEqual(page.count("<summary>Technical details</summary>"), 9)
        self.assertIn("A known genome used to verify positions and compare mutations.", page)


class LocalPreviewTest(unittest.TestCase):
    """The preview script runs in its own process, so its extra routes never touch this app."""

    def test_every_preview_state_opens_and_is_labelled_synthetic(self):
        script = (
            "import re, sys\n"
            "sys.argv = ['preview']\n"
            "import runpy\n"
            "module = runpy.run_path(r'%s', run_name='preview')\n"
            "client = module['app'].test_client()\n"
            "for state, _ in module['STATES']:\n"
            "    form = client.get('/preview/sprint4/' + state).get_data(as_text=True)\n"
            "    action = re.search(r\"action='([^']+)'\", form).group(1)\n"
            "    data = dict(re.findall(r\"name='([^']+)' value='([^']*)'\", form))\n"
            "    page = client.post(action, data=data).get_data(as_text=True)\n"
            "    print(state, 'LOCAL UI PREVIEW' in page, 'SYNTHETIC' in page,\n"
            "          'Find and Check Publications' in page or 'review-outcome' in page)\n"
        ) % (ROOT / "validation" / "local_ui_preview.py")
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, timeout=300,
                                cwd=ROOT, env={**os.environ, "AI_ENRICHMENT_ENABLED": "false"})
        self.assertEqual(result.returncode, 0, result.stderr[-2000:])
        lines = [line.split() for line in result.stdout.splitlines()]
        self.assertEqual([line[0] for line in lines], ["before", "one-publication", "general-mention",
                                                       "no-publication", "eligible", "retrieval-error"])
        for line in lines:
            # Every page carries the preview banner; every Sprint 4 page also the app's own SYNTHETIC label.
            expected = ["True", "False" if line[0] == "before" else "True", "True"]
            self.assertEqual(line[1:], expected, line)


if __name__ == "__main__":
    unittest.main()
