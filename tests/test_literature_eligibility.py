"""Which publications may be analyzed as evidence about a sample, and safe preparation errors.

Reproduces the real ten-record NCBI case (accessions QB228092.1-QB228101.1,
collected 2026-09-07, linked through a BioProject to a 2022 paper that names
none of them) with an in-memory literature source: invented article text,
no NCBI, PubTator, or Claude request. Claude calls are replaced by mocks that
must never be reached when no article is usable.
"""

import io
import os
import re
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402

import literature_view  # noqa: E402
from app import STAGED, app  # noqa: E402
from literature.chronology import chronology, earliest_day, latest_day  # noqa: E402
from literature.config import SDK_MISSING, ai_settings, ai_status  # noqa: E402
from literature.discovery import DIRECT, STUDY, article_overview  # noqa: E402
from literature.passages import select_passages  # noqa: E402
from literature.providers import ClaudeProvider, ProviderError  # noqa: E402
from literature.schema import EXTRACTION_SCHEMA  # noqa: E402
from literature.workflow import INDIVIDUAL, find_literature, new_session  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
ACCESSIONS = [f"QB228{number:03d}" for number in range(92, 102)]
OLD_STUDY = "PMID:9000001"
NO_EVIDENCE = ("One study-associated publication was found, but it does not provide verifiable collection "
               "information for these samples. No metadata changes are proposed")
REVIEWED_NO_CHANGE = ("10 accession records searched. One publication found. No verified collection locations could be "
                      "assigned to these samples. No metadata was changed.")
OFFLINE_AI = {"AI_ENRICHMENT_ENABLED": "true", "ANTHROPIC_API_KEY": "test-key-not-real",
              "ANTHROPIC_MODEL": "claude-opus-5-5", "AI_ACCESS_CODE": "", "GEONAMES_USERNAME": ""}


def sdk(installed):
    """Settings as if the 'anthropic' package were (or were not) installed; the Claude calls are mocked."""

    return patch("literature_view.ai_settings", side_effect=lambda: {**ai_settings(), "sdk_installed": installed})


def passage(article_id, number, section, text, kind="paragraph"):
    return {"Article_ID": article_id, "Passage_ID": f"{article_id}#p{number}", "Section": section, "Type": kind,
            "Offset": 0, "Text": text, "Annotations": []}


class FakeSource:
    """Answers like LiveSource from in-memory records (invented text, no network)."""

    label = "In-memory test articles"
    synthetic = False

    def __init__(self, articles, links, aliases=None):
        self.articles, self.links, self.aliases = articles, links, aliases or {}

    def record_links(self, kind, value):
        return [(article_id, relationship, f"Test record {value} links {article_id}.", "Test record link")
                for article_id, relationship in self.links.get(str(value).split(".")[0], [])]

    def record_aliases(self, kind, value):
        return self.aliases.get(str(value).split(".")[0], [])

    def searchable(self, alias):
        return False

    def text_search(self, alias):
        return []

    def describe(self, article_ids):
        return {article_id: {"Title": self.articles[article_id]["Title"],
                             "Publication_Date": self.articles[article_id]["Publication_Date"],
                             "Epub_Date": self.articles[article_id].get("Epub_Date", ""),
                             "PMID": article_id.split(":")[1], "PMCID": "", "DOI": "",
                             "Source_URL": f"https://pubmed.ncbi.nlm.nih.gov/{article_id.split(':')[1]}/"}
                for article_id in article_ids if article_id in self.articles}

    def fetch(self, article_ids):
        return {article_id: {"Passages": self.articles[article_id]["Passages"], "Text_Status": "Full text (test)"}
                for article_id in article_ids if article_id in self.articles}


def old_study_article():
    """A 2022 surveillance paper that names no 2026 sample (like PMID 35302933)."""

    return {"Title": "Invented wild bird surveillance study, 2021", "Publication_Date": "2022 May",
            "Epub_Date": "2022 Mar 18", "Passages": [
                passage(OLD_STUDY, 1, "METHODS", "Wild bird samples were collected in December 2021 in two "
                                                 "coastal states and screened for influenza A virus."),
                passage(OLD_STUDY, 2, "RESULTS", "A sample collected on 30 December 2021 from a duck in a "
                                                 "coastal county was positive for H5N1.")]}


def ten_records(accessions=ACCESSIONS, date="2026-09-07"):
    """Like the real file: QB228092-099 are isolate 1 (8 segments), QB228100-101 isolate 2."""

    return pd.DataFrame({"Sample_ID": [f"{accession}.1" for accession in accessions],
                         "Accession": [f"{accession}.1" for accession in accessions],
                         "Isolate": [f"A/Turkey/SD/TEST-{1 if accession < 'QB228100' else 2}/2026"
                                     for accession in accessions],
                         "Location_Name": ["USA: SD"] * len(accessions),
                         "Collection_Date": [date] * len(accessions)})


def discover(source, samples=None):
    session = new_session(individual=samples if samples is not None else ten_records())
    find_literature(session, source)
    return session


def ten_record_source():
    return FakeSource({OLD_STUDY: old_study_article()},
                      {accession: [(OLD_STUDY, "study")] for accession in ACCESSIONS})


class ChronologyTest(unittest.TestCase):
    def test_partial_pubmed_dates_mean_their_latest_day(self):
        self.assertEqual(latest_day("2022 Mar 18"), pd.Timestamp(2022, 3, 18))
        self.assertEqual(latest_day("2022 Mar-Apr"), pd.Timestamp(2022, 4, 30))
        self.assertEqual(latest_day("2022 May"), pd.Timestamp(2022, 5, 31))
        self.assertEqual(latest_day("2022 Spring"), pd.Timestamp(2022, 12, 31))
        self.assertIsNone(latest_day(""))
        self.assertEqual(earliest_day("2026"), pd.Timestamp(2026, 1, 1))
        self.assertEqual(earliest_day("2026-09"), pd.Timestamp(2026, 9, 1))

    def test_only_a_certainly_earlier_publication_is_ruled_out(self):
        self.assertTrue(chronology(["2022 May", "2022 Mar 18"], "2026-09-07")[0])
        self.assertIn("Published 2022 Mar 18", chronology(["2022 May", "2022 Mar 18"], "2026-09-07")[1])
        self.assertFalse(chronology(["2026 Nov"], "2026-09-07")[0])
        self.assertFalse(chronology(["2026"], "2026-09-07")[0])  # same year: cannot be ruled out
        self.assertIsNone(chronology(["2022 May"], None)[0])
        self.assertIsNone(chronology([""], "2026-09-07")[0])


class EligibilityTest(unittest.TestCase):
    def test_ten_accessions_linked_only_to_an_older_study_paper(self):
        session = discover(ten_record_source())
        table = session["publications"]
        self.assertEqual(len(table), 10)
        self.assertTrue((table["Relationship_Type"] == STUDY).all())
        self.assertFalse(table["Eligible_For_Enrichment"].any())
        self.assertIn("before this sample was collected (2026-09-07)", table["Eligibility_Reason"].iloc[0])
        self.assertIn("Linked only through the broader study", table["Eligibility_Reason"].iloc[0])
        self.assertEqual(table["Epub_Date"].iloc[0], "2022 Mar 18")  # publication metadata kept
        self.assertEqual(table["Association_Evidence"].iloc[0], "Test record QB228092.1 links PMID:9000001.")
        self.assertEqual(session["plan"], [])
        self.assertTrue(session["no_evidence"]["message"].startswith(NO_EVIDENCE))
        self.assertEqual(session["no_evidence"]["reasons"][0]["links"], "10 sample links")
        overview = article_overview(table, 10)
        self.assertEqual((overview["articles"], overview["direct"], overview["study_only"], overview["named"],
                          overview["usable"], overview["unsuitable"], overview["samples_without"]),
                         (1, 0, 1, 0, 0, 1, 10))

    def test_paper_directly_linked_to_a_sample(self):
        article = {"Title": "Invented outbreak report", "Publication_Date": "2026 Nov", "Passages": [
            passage("PMID:9000002", 1, "METHODS", "Turkeys were sampled on a farm in Kingsbury County, "
                                                  "South Dakota, on 7 September 2026.")]}
        source = FakeSource({"PMID:9000002": article}, {ACCESSIONS[0]: [("PMID:9000002", "direct")]})
        session = discover(source)
        row = session["publications"].iloc[0]
        self.assertEqual(row["Relationship_Type"], DIRECT)
        self.assertTrue(row["Eligible_For_Enrichment"])
        self.assertEqual(row["Eligibility_Reason"], "Usable: the sample's own record links this article.")
        self.assertEqual([item["Article_ID"] for item in session["plan"]], ["PMID:9000002"])
        self.assertIsNone(session["no_evidence"])

    def test_study_paper_with_a_sample_to_site_table_is_used_only_for_the_named_sample(self):
        article = {"Title": "Invented study with a sample table", "Publication_Date": "2026 Dec", "Passages": [
            passage("PMID:9000003", 1, "METHODS", "Samples for this study came from several farms in South Dakota."),
            passage("PMID:9000003", 2, "TABLE", "Isolate | Farm site | Collected\nA/Turkey/SD/TEST-1/2026 | "
                                                "Farm 4, Kingsbury County | 2026-09-07", kind="table")]}
        samples = ten_records(ACCESSIONS[:1] + ACCESSIONS[8:9])  # isolates TEST-1 and TEST-2
        source = FakeSource({"PMID:9000003": article},
                            {accession: [("PMID:9000003", "study")] for accession in ACCESSIONS})
        session = discover(source, samples)
        table = session["publications"].set_index("Sample_ID")
        self.assertTrue(table.loc["QB228092.1", "Eligible_For_Enrichment"])  # named in the table
        self.assertEqual(table.loc["QB228092.1", "Eligibility_Reason"], "Usable: the article text names this sample.")
        self.assertFalse(table.loc["QB228100.1", "Eligible_For_Enrichment"])  # same study, not named
        plan = session["plan"][0]
        self.assertEqual([sample["Sample_ID"] for sample in plan["Samples"]], ["QB228092.1"])
        self.assertIn("PMID:9000003#p2", [item["Passage_ID"] for item in plan["Passages"]])
        self.assertTrue(any("1 other linked sample(s) are left out (QB228100.1)" in note
                            for note in session["plan_notes"]))

    def test_study_paper_without_sample_identifiers_is_not_used(self):
        article = {"Title": "Invented study overview", "Publication_Date": "2026 Dec", "Passages": [
            passage("PMID:9000004", 1, "METHODS", "Samples were collected from turkey farms in South Dakota "
                                                  "during September 2026.")]}
        source = FakeSource({"PMID:9000004": article},
                            {accession: [("PMID:9000004", "study")] for accession in ACCESSIONS})
        session = discover(source)
        reason = session["publications"]["Eligibility_Reason"].iloc[0]
        self.assertNotIn("before this sample was collected", reason)  # published after collection
        self.assertIn("no identifier of this sample", reason)
        self.assertEqual(session["plan"], [])  # the study location is not given to every sample
        self.assertTrue(session["no_evidence"]["message"].startswith(NO_EVIDENCE))

    def test_biosample_accession_from_ncbi_records_matches_a_table(self):
        article = {"Title": "Invented BioSample table", "Publication_Date": "2026 Dec", "Passages": [
            passage("PMID:9000005", 1, "SUPPL", "BioSample | Site\nSAMN99000001 | Farm 7, Brookings County, "
                                                "South Dakota", kind="table")]}
        source = FakeSource({"PMID:9000005": article}, {ACCESSIONS[0]: [("PMID:9000005", "study")]},
                            aliases={ACCESSIONS[0]: ["SAMN99000001"]})
        session = discover(source, ten_records(ACCESSIONS[:1]))
        self.assertTrue(session["publications"]["Eligible_For_Enrichment"].iloc[0])
        self.assertIn("SAMN99000001", session["aliases"]["QB228092.1"])
        self.assertEqual(len(session["plan"]), 1)

    def test_direct_link_published_before_collection_is_not_used(self):
        article = {"Title": "Invented earlier report", "Publication_Date": "2023 Jan", "Passages": [
            passage("PMID:9000006", 1, "METHODS", "Samples were collected in Brookings County in 2022.")]}
        source = FakeSource({"PMID:9000006": article}, {ACCESSIONS[0]: [("PMID:9000006", "direct")]})
        session = discover(source, ten_records(ACCESSIONS[:1]))
        row = session["publications"].iloc[0]
        self.assertFalse(row["Eligible_For_Enrichment"])
        self.assertIn("Published 2023 Jan, before this sample was collected", row["Eligibility_Reason"])
        self.assertTrue(session["no_evidence"]["message"].startswith("One linked publication was found"))

    def test_passages_that_name_the_sample_are_preferred(self):
        article = {"Passages": [passage("A", 1, "METHODS", "Samples were collected at three sites."),
                                passage("A", 2, "TABLE", "QB228092 | Farm 4", kind="table")]}
        chosen = select_passages(article, ["QB228092"])
        self.assertIn("A#p2", [item["Passage_ID"] for item in chosen])


class RouteTest(unittest.TestCase):
    """The Sprint 4 pages with mocked Claude calls (never reached for unusable articles)."""

    def setUp(self):
        patches = [patch.dict(os.environ, OFFLINE_AI),
                   patch("literature.config.ENV_FILE", ROOT / "missing.env"),
                   patch("align_sequences.shutil.which", return_value="mafft"),
                   patch.dict(app.config, {"LITERATURE_INLINE_JOBS": True})]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        self.addCleanup(STAGED.clear)
        self.client = app.test_client()

    def session_token(self, source, samples=None):
        session = new_session(individual=samples if samples is not None else ten_records())
        session["samples"]["Sample_Source"] = INDIVIDUAL
        token = literature_view.store_session(session)
        with patch("literature_view.literature_source", return_value=source):
            page = self.client.post("/literature/find", data={"token": token}).get_data(as_text=True)
        return token, page

    def eligible_source(self):
        article = {"Title": "Invented outbreak report", "Publication_Date": "2026 Nov", "Passages": [
            passage("PMID:9000002", 1, "METHODS", "Sample QB228092 was collected on a turkey farm in "
                                                  "Kingsbury County, South Dakota.")]}
        return FakeSource({"PMID:9000002": article}, {ACCESSIONS[0]: [("PMID:9000002", "direct")]})

    def test_no_eligible_publication_shows_a_clear_result_and_charges_nothing(self):
        with patch.object(ClaudeProvider, "count_tokens") as count, patch.object(ClaudeProvider, "run") as run:
            token, page = self.session_token(ten_record_source())
            analyzed = self.client.post("/literature/analyze", data={"token": token, "confirm": "1",
                                                                       "acknowledge": "on"}).get_data(as_text=True)
        self.assertIn(REVIEWED_NO_CHANGE, page)
        self.assertIn("Publications with study-level links only:</strong> 1", page)
        self.assertIn("Publications with usable sample-specific evidence:</strong> 0", page)
        self.assertIn("Samples without a supported improvement:</strong> 10", page)
        self.assertNotIn("Analyze Literature with Claude AI", page)
        self.assertNotIn("Analysis could not be completed", analyzed)
        self.assertIn(REVIEWED_NO_CHANGE, analyzed)
        count.assert_not_called()
        run.assert_not_called()

    def test_missing_anthropic_sdk_is_explained_before_any_request(self):
        """The root cause of the real failure: the server's Python had no 'anthropic' package."""

        token, _ = self.session_token(self.eligible_source())
        with patch("literature.config.importlib.util.find_spec", return_value=None):
            self.assertEqual(ai_status()["message"], SDK_MISSING)
            self.assertFalse(ai_settings()["sdk_installed"])
        with sdk(installed=False):
            page = self.client.post("/literature/analyze", data={"token": token}).get_data(as_text=True)
        self.assertIn("The Anthropic Python SDK (the &#39;anthropic&#39; package listed in requirements.txt) is "
                      "not installed", page)
        self.assertNotIn("The analysis could not be prepared.", page)

    def test_preparation_errors_are_specific_and_never_leak_details(self):
        token, _ = self.session_token(self.eligible_source())
        problem = ProviderError("The Claude API rejected the configured key while counting tokens (HTTP 401).")
        with patch.object(ClaudeProvider, "count_tokens", side_effect=problem), sdk(installed=True):
            page = self.client.post("/literature/analyze", data={"token": token}).get_data(as_text=True)
        self.assertIn("The analysis could not be prepared: The Claude API rejected the configured key while "
                      "counting tokens (HTTP 401). Nothing was sent for extraction and nothing was billed.", page)
        with patch("app.prepare", side_effect=RuntimeError("secret internal detail")), \
                self.assertLogs(app.logger, "ERROR"):
            page = self.client.post("/literature/analyze", data={"token": token}).get_data(as_text=True)
        self.assertIn("unexpected server error (RuntimeError). Nothing was sent to Claude and nothing was billed",
                      page)
        self.assertNotIn("secret internal detail", page)

    def test_eligible_evidence_still_reaches_the_cost_confirmation(self):
        with patch.object(ClaudeProvider, "count_tokens", return_value=(1500, "Counted (test)")) as count, \
                patch.object(ClaudeProvider, "run") as run, sdk(installed=True):
            token, page = self.session_token(self.eligible_source())
            self.assertIn("Analyze Literature with Claude AI", page)
            confirm = self.client.post("/literature/analyze", data={"token": token}).get_data(as_text=True)
        self.assertNotIn("Analysis could not be completed", confirm)
        self.assertIn("Input tokens:</strong> about 1,500", confirm)  # shown before anything is billed
        count.assert_called()
        run.assert_not_called()  # nothing runs without explicit confirmation

    def test_continuing_with_the_original_metadata_runs_sprint_2(self):
        metadata = (FIXTURES / "ncbi_virus_structure_metadata.csv").read_bytes()
        sequences = (FIXTURES / "ncbi_virus_structure_sequences.fasta").read_text(encoding="utf-8")
        first = sequences.split("\n")[1]
        source = FakeSource({OLD_STUDY: old_study_article()},
                            {f"TST{number:06d}": [(OLD_STUDY, "study")] for number in range(1, 11)})
        files = [(io.BytesIO(sequences.encode()), "sequences.fasta"), (io.BytesIO(metadata), "metadata.csv")]
        page = self.client.post("/review", data={"files": files, "enrich": "on", "intended_mode": "individual"},
                                content_type="multipart/form-data").get_data(as_text=True)
        staged = re.search(r'name="token" value="([^"]+)"', page).group(1)
        with redirect_stdout(io.StringIO()), patch("literature_view.literature_source", return_value=source), \
                patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(first)):
            page = self.client.post("/enrich/start", data={"token": staged, "enrich": "on",
                                                           "enrich_shown": "1"}).get_data(as_text=True)
            token = re.search(r'name="token" value="([^"]+)"', page[page.index('id="literature"'):]).group(1)
            found = self.client.post("/literature/find", data={"token": token}).get_data(as_text=True)
            result = self.client.post("/enrich/run", data={"token": token, "skip": "1"}).get_data(as_text=True)
        self.assertIn(REVIEWED_NO_CHANGE, found)
        self.assertIn("Continue with the original metadata", found)
        self.assertIn("Pairwise comparisons:</strong> 2", result)
        self.assertNotIn('id="enrichment-report"', result)  # nothing was applied


@unittest.skipUnless(__import__("importlib").util.find_spec("anthropic"), "anthropic SDK not installed")
class TokenCountingErrorsTest(unittest.TestCase):
    def provider(self, error):
        client = MagicMock()
        client.messages.count_tokens.side_effect = error
        settings = {"model": "claude-opus-5-5", "effort": "low", "refusal_fallback": True}
        return ClaudeProvider(settings, client)

    def test_problems_the_paid_request_would_hit_stop_before_billing(self):
        import anthropic
        import httpx2
        request = httpx2.Request("POST", "https://api.anthropic.com/v1/messages/count_tokens")
        for error, words in [
                (anthropic.AuthenticationError("secret detail", response=httpx2.Response(401, request=request),
                                               body=None), "rejected the configured key"),
                (anthropic.NotFoundError("secret detail", response=httpx2.Response(404, request=request),
                                         body=None), "not available to this API key"),
                (anthropic.BadRequestError("secret detail", response=httpx2.Response(400, request=request),
                                           body=None), "token-count request as invalid"),
                (anthropic.AnthropicError("secret detail"), "client could not be created")]:
            with self.assertRaises(ProviderError) as raised:
                self.provider(error).count_tokens("s", "u", EXTRACTION_SCHEMA)
            self.assertIn(words, str(raised.exception))
            self.assertNotIn("secret detail", str(raised.exception))

    def test_temporary_outages_fall_back_to_an_approximation(self):
        import anthropic
        import httpx2
        request = httpx2.Request("POST", "https://api.anthropic.com/v1/messages/count_tokens")
        tokens, method = self.provider(anthropic.APIConnectionError(request=request)).count_tokens(
            "system", "user text", EXTRACTION_SCHEMA)
        self.assertGreater(tokens, 0)
        self.assertIn("Approximate (token counting unavailable: APIConnectionError)", method)


if __name__ == "__main__":
    unittest.main()
