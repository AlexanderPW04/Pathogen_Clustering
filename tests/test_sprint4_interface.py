"""The Sprint 4 page answers direct questions, in order, with few clicks.

Mocked external services only (the in-memory source from
test_literature_eligibility and the offline mock provider): no NCBI, PubTator,
or Claude request is sent, and nothing is billed.
"""

import io
import json
import os
import re
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402
from test_literature_eligibility import (  # noqa: E402
    OFFLINE_AI, FakeSource, old_study_article, passage, ten_record_source, ten_records)

import literature_view  # noqa: E402
from app import STAGED, app  # noqa: E402
from literature.http_client import ServiceUnavailable  # noqa: E402
from literature.providers import ClaudeProvider, MockProvider  # noqa: E402
from literature.review import plain_explanation, PREDATES, review_summary  # noqa: E402
from literature.workflow import INDIVIDUAL, find_literature, new_session, run_analysis  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
ACCESSION = "QB228092"
PREDATES_SENTENCE = ("This article was published before these samples were collected and does not identify "
                     "them. Its collection locations cannot be applied to your samples.")


def text(page):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", page)).replace("&#39;", "'").replace("&quot;", '"')


def primary_part(page):
    """Everything above the expandable technical details."""

    start = page.index('id="literature"')
    return page[start:page.index('id="technical-details"')]


def article(number, sentence, date="2026 Dec", section="METHODS"):
    article_id = f"PMID:{number}"
    return article_id, {"Title": f"Invented article {number}", "Publication_Date": date,
                        "Passages": [passage(article_id, 1, section, sentence)]}


def review_of(source, samples=None):
    session = new_session(individual=samples if samples is not None else ten_records([ACCESSION]))
    find_literature(session, source)
    return session, review_summary(session)


class AnswersTest(unittest.TestCase):
    def one(self, summary):
        self.assertEqual(len(summary["cards"]), 1)
        return summary["cards"][0]

    def test_no_publication_found(self):
        session, summary = review_of(FakeSource({}, {}))
        self.assertEqual(summary["searched"], "1 accession record searched.")
        self.assertEqual(summary["headline"], "No publication was found for the uploaded samples. No metadata was changed.")
        self.assertEqual(summary["claude"]["why"], "Claude was not called because no publication was found.")
        self.assertEqual(summary["samples"][0]["Status"], "No associated publication found")

    def test_article_with_no_location_information(self):
        identifier, item = article(9400001, f"Sample {ACCESSION} was sequenced twice.")
        _, summary = review_of(FakeSource({identifier: item}, {ACCESSION: [(identifier, "direct")]}))
        card = self.one(summary)
        self.assertEqual((card["location"], card["places"]), ("No", []))

    def test_place_named_in_an_article_is_not_a_location_for_the_sample(self):
        identifier, item = article(9400002, "Samples for this program came from Maryland and Arizona.")
        _, summary = review_of(FakeSource({identifier: item}, {ACCESSION: [(identifier, "study")]}))
        card = self.one(summary)
        self.assertEqual(card["location"], "Yes")
        self.assertEqual(card["places"], ["Maryland", "Arizona"])
        self.assertEqual((card["applies"], card["improvement"], card["claude"]), ("No", "No", "No"))
        self.assertIn("does not name your samples", card["plain"])

    def test_affiliations_do_not_count_as_location_information(self):
        identifier, item = article(9400003, "Department of Biology, University of Arizona, Tucson, Arizona, USA.",
                                   section="AFFILIATION")
        _, summary = review_of(FakeSource({identifier: item}, {ACCESSION: [(identifier, "direct")]}))
        self.assertEqual(self.one(summary)["location"], "No")

    def test_publication_predating_collection(self):
        _, summary = review_of(ten_record_source(), ten_records())
        card = self.one(summary)
        self.assertEqual(card["plain"], PREDATES_SENTENCE)
        self.assertEqual((card["year"], card["applies"], card["improvement"], card["claude"]), ("2022", "No", "No", "No"))
        # The exact wording asked for: reviewed (retrieved with text and assessed), predates, nothing changed.
        self.assertEqual(summary["searched"], "10 accession records searched.")
        self.assertEqual(summary["headline"], "One publication was found and reviewed. It predates the uploaded "
                                              "samples and does not provide verifiable collection locations for "
                                              "them. No metadata was changed.")
        self.assertEqual(summary["claude_line"], "Claude was not called because no eligible sample-specific "
                                                 "evidence was found.")
        self.assertTrue(card["text_checked"])
        self.assertEqual(card["changes"], "No")
        for line in ["Ten accessions represent two distinct isolates.", "One associated publication was retrieved.",
                     "The publication is from 2022; the samples were collected in September 2026.",
                     "No usable direct sample-specific evidence was established."]:
            self.assertIn(line, summary["why"])
        self.assertIn(("Claude not run", "skipped"), summary["stages"])
        self.assertIn(("No changes supported", "done"), summary["stages"])

    def test_retrieval_failure_says_could_not_check_never_no_location(self):
        class Failing(FakeSource):
            def fetch(self, article_ids):
                raise ServiceUnavailable("PubTator3 returned HTTP 503 after 3 attempt(s).")

        identifier, item = article(9400004, "Sample QB228092 was collected in Tempe, Arizona.")
        _, summary = review_of(Failing({identifier: item}, {ACCESSION: [(identifier, "direct")]}))
        card = self.one(summary)
        self.assertEqual((card["location"], card["applies"], card["improvement"]),
                         ("Not checked", "Uncertain", "Not checked"))
        self.assertFalse(card["text_checked"])
        self.assertNotIn("reviewed", summary["headline"])
        self.assertIn("could not be retrieved", summary["headline"])
        self.assertNotIn("No verified collection locations", summary["headline"])
        self.assertIn(("Article retrieval failed", "failed"), summary["stages"])
        self.assertEqual(summary["claude"]["why"], "Claude was not called because the evidence check could not be completed.")

    def test_several_publications_get_one_row_each(self):
        first = article(9400005, "Samples came from Ohio.")
        second = article(9400006, "Samples came from Texas.")
        _, summary = review_of(FakeSource(dict([first, second]), {ACCESSION: [(first[0], "study"), (second[0], "study")]}))
        self.assertEqual([card["article"] for card in summary["cards"]], [first[0], second[0]])
        self.assertEqual(summary["claude"]["why"], "Claude was not called because no eligible sample-specific "
                                                   "evidence was found.")
        self.assertEqual(summary["headline"], "Two publications were found and reviewed. They do not provide "
                                              "verifiable collection locations for the uploaded samples. No metadata "
                                              "was changed.")

    def test_plain_sentences_avoid_technical_terms(self):
        sentence = plain_explanation(PREDATES, [], 0)
        for word in ["bioproject", "eligib", "provenance", "identifier"]:
            self.assertNotIn(word, sentence.casefold())


class ClaudeStateTest(unittest.TestCase):
    def analyzed(self, responses):
        identifier, item = article(9400007, f"Sample {ACCESSION} was collected in Brookings, South Dakota.")
        session = new_session(individual=ten_records([ACCESSION]))
        find_literature(session, FakeSource({identifier: item}, {ACCESSION: [(identifier, "direct")]}))
        before = review_summary(session)
        self.assertEqual(before["cards"][0]["improvement"], "Needs Claude review")
        self.assertEqual(before["claude"]["label"], "Claude not run")
        with tempfile.TemporaryDirectory() as folder:
            mock = Path(folder) / "mock.json"
            mock.write_text(json.dumps({"extraction": {identifier: responses}}), encoding="utf-8")
            run_analysis(session, MockProvider(mock), {**literature_view.ai_settings(), "max_articles": 5,
                                                       "max_requests": 10, "max_input_tokens": 60000})
        return review_summary(session)

    def test_verified_collection_site(self):
        summary = self.analyzed({"facts": [{
            "Sample_ID": f"{ACCESSION}.1", "Field": "Location", "Value": "Brookings, South Dakota, USA",
            "Location_Role": "Sample collection site", "Geographic_Level": "City or town",
            "Passage_ID": "PMID:9400007#p1",
            "Evidence_Quote": f"Sample {ACCESSION} was collected in Brookings, South Dakota.",
            "Evidence_Type": "Directly stated", "Ambiguity": "", "Reasoning_Summary": "Named."}],
            "insufficient_evidence_samples": [], "notes": ""})
        card = summary["cards"][0]
        self.assertEqual((card["applies"], card["improvement"], card["claude"]), ("Yes", "Yes", "Yes"))
        self.assertIn(("Claude analysis completed", "done"), summary["stages"])
        self.assertIn(("Metadata changes proposed", "done"), summary["stages"])
        self.assertIn("Verified collection details were proposed for one sample", summary["headline"])

    def test_claude_found_nothing(self):
        summary = self.analyzed({"facts": [], "insufficient_evidence_samples": [f"{ACCESSION}.1"], "notes": ""})
        card = summary["cards"][0]
        self.assertEqual((card["applies"], card["improvement"], card["claude"]), ("No", "No", "Yes"))
        self.assertIn(("No changes supported", "done"), summary["stages"])


class WorkflowTest(unittest.TestCase):
    def setUp(self):
        patches = [patch.dict(os.environ, OFFLINE_AI), patch("literature.config.ENV_FILE", ROOT / "missing.env"),
                   patch("align_sequences.shutil.which", return_value="mafft"),
                   patch.dict(app.config, {"LITERATURE_INLINE_JOBS": True})]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        self.addCleanup(STAGED.clear)
        self.client = app.test_client()
        self.sequences = (FIXTURES / "ncbi_virus_structure_sequences.fasta").read_text(encoding="utf-8")
        self.metadata = (FIXTURES / "ncbi_virus_structure_metadata.csv").read_bytes()

    def post(self, url, data):
        with redirect_stdout(io.StringIO()):
            response = self.client.post(url, data=data)
        if response.status_code == 302:
            response = self.client.get(response.headers["Location"])
        return response.get_data(as_text=True)

    def staged(self):
        files = [(io.BytesIO(self.sequences.encode()), "s.fasta"), (io.BytesIO(self.metadata), "m.csv")]
        page = self.client.post("/review", data={"files": files, "enrich": "on", "intended_mode": "individual"},
                                content_type="multipart/form-data").get_data(as_text=True)
        return page, re.search(r'name="token" value="([^"]+)"', page).group(1)

    def source(self, sentence=None):
        accessions = [f"TST{number:06d}" for number in range(1, 11)]
        if sentence is None:
            return FakeSource({"PMID:9000001": old_study_article()},
                              {accession: [("PMID:9000001", "study")] for accession in accessions})
        identifier, item = article(9400008, sentence)
        return FakeSource({identifier: item}, {"TST000001": [(identifier, "direct")]})

    def test_one_click_reaches_the_conclusion_and_the_next_action(self):
        page, staged = self.staged()
        self.assertIn("Find and Check Publications", page)
        calls = []
        source = self.source()
        original = source.fetch
        source.fetch = lambda ids: calls.append(list(ids)) or original(ids)
        with patch("literature_view.literature_source", return_value=source), \
                patch.object(ClaudeProvider, "count_tokens") as count, patch.object(ClaudeProvider, "run") as run:
            page = self.post("/enrich/start", {"token": staged, "enrich": "on", "enrich_shown": "1"})
            again = self.post("/enrich/start", {"token": staged, "enrich": "on", "enrich_shown": "1"})
        primary = text(primary_part(page))
        self.assertIn("10 accession records searched.", primary)
        self.assertIn("One publication was found and reviewed. It predates the uploaded samples", primary)
        self.assertIn("Claude was not called because no eligible sample-specific evidence was found.", primary)
        self.assertIn(PREDATES_SENTENCE, primary)
        for question in ["Geographic information found?", "Linked to uploaded samples?",
                         "Metadata improvement supported?", "Claude analysis performed?", "Changes applied?"]:
            self.assertIn(question, primary)
        self.assertIn("Claude not run", primary)
        self.assertLess(primary.index("10 accession records searched"), primary.index("Geographic information found?"))
        self.assertNotIn("bioproject_pubmed", primary_part(page))
        self.assertNotIn("Link provenance", primary_part(page))
        self.assertIn("Link provenance", page)  # kept, in the details
        self.assertIn('class="action-bar"', page)
        self.assertLess(page.index('id="technical-details"'), page.index('class="action-bar"'))
        self.assertIn("Continue to Genomic Analysis", page)
        self.assertIn('form="back-to-review"', page)
        self.assertNotIn("Analyze Literature with Claude AI", page)
        count.assert_not_called()
        run.assert_not_called()
        # Returning from the file review reopens the findings instead of searching again.
        self.assertEqual(len(calls), 1)
        self.assertIn("10 accession records searched", again)

    def test_continue_with_original_metadata_reports_zero_changes(self):
        _, staged = self.staged()
        first = self.sequences.split("\n")[1]
        with patch("literature_view.literature_source", return_value=self.source()), \
                patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(first)):
            page = self.post("/enrich/start", {"token": staged, "enrich": "on", "enrich_shown": "1"})
            token = re.search(r'name="token" value="([^"]+)"', page[page.index('id="literature"'):]).group(1)
            result = self.post("/enrich/run", {"token": token, "skip": "1"})
        self.assertIn("Pairwise comparisons:</strong> 2", result)
        self.assertIn("Sprint 4: 10 accession records searched. One publication was found and reviewed.",
                      text(result).replace("  ", " "))
        self.assertIn("No metadata was changed.", text(result))

    def proposals(self):
        _, staged = self.staged()
        sentence = "Sample TST000001 was collected at a farm in Kingsbury County, South Dakota."
        mock_facts = {"extraction": {"PMID:9400008": {"facts": [{
            "Sample_ID": "TST000001.1", "Field": "Location", "Value": "Kingsbury County, South Dakota, USA",
            "Location_Role": "Sample collection site", "Geographic_Level": "County",
            "Passage_ID": "PMID:9400008#p1", "Evidence_Quote": sentence, "Evidence_Type": "Directly stated",
            "Ambiguity": "", "Reasoning_Summary": "Named."}], "insufficient_evidence_samples": [], "notes": ""}}}
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        mock = Path(folder.name) / "mock.json"
        mock.write_text(json.dumps(mock_facts), encoding="utf-8")
        with patch("literature_view.literature_source", return_value=self.source(sentence)), \
                patch("literature_view.choose_provider", return_value=(MockProvider(mock), "Mock (test)")):
            page = self.post("/enrich/start", {"token": staged, "enrich": "on", "enrich_shown": "1"})
            token = re.search(r'name="token" value="([^"]+)"', page[page.index('id="literature"'):]).group(1)
            self.assertIn("Analyze Literature with Claude AI", page)
            self.assertIn("Needs Claude review", page)
            estimate = self.post("/literature/analyze", {"token": token})
            self.assertIn("Before you start", estimate)
            reviewed = self.post("/literature/analyze", {"token": token, "confirm": "1"})
        return token, reviewed

    def finish(self, token, choice):
        session = literature_view.get_session(token)
        fact = session["facts"][0]["Fact_ID"]
        first = self.sequences.split("\n")[1]
        with patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(first)):
            return self.post("/enrich/run", {"token": token, f"review_{fact}": choice})

    def test_authorized_run_and_approval(self):
        token, reviewed = self.proposals()
        self.assertIn("Claude analysis completed", reviewed)
        self.assertIn("Approve all verified proposals", reviewed)
        self.assertIn('form="review-form"', reviewed)
        result = self.finish(token, "approve")
        self.assertIn("Kingsbury County, South Dakota, USA", result)
        self.assertIn("Values changed:</strong> 1", result)

    def test_rejected_proposal_changes_nothing(self):
        token, _ = self.proposals()
        result = self.finish(token, "reject")
        self.assertIn("Values changed:</strong> 0", result)
        self.assertIn("Pairwise comparisons:</strong> 2", result)

    def test_samples_with_different_results_get_a_table(self):
        session = new_session(individual=ten_records(["QB228092", "QB228093"]))
        session["samples"]["Sample_Source"] = INDIVIDUAL
        identifier, item = article(9400009, "Sample QB228092 came from Tempe, Arizona.")
        token = literature_view.store_session(session)
        source = FakeSource({identifier: item}, {"QB228092": [(identifier, "direct")]})
        with patch("literature_view.literature_source", return_value=source):
            page = self.client.post("/literature/find", data={"token": token}).get_data(as_text=True)
        self.assertIn("Result for each sample (2)", page)
        self.assertIn("No associated publication found", page)


if __name__ == "__main__":
    unittest.main()
