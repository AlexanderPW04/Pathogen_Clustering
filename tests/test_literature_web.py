"""Flask tests for the optional literature-enrichment section (Step 4).

MAFFT is mocked with tests/fake_mafft.py, GeoNames is disabled, literature
comes from the synthetic corpus (or a fake live source), and the Claude
provider is replaced by a fake client. No network request and no paid API
request is made, so these tests do not show that live services work.
"""

import io
import os
import re
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import MagicMock, patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402

import literature_view  # noqa: E402
from app import STATUS_CACHE, app  # noqa: E402
from literature import budget as budget_module  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "literature"
FIRST_SEQUENCE = (EXAMPLE / "synthetic_literature_sequences.fasta").read_text(encoding="utf-8").split("\n")[1]
SECRET = "fake-test-key-NEVER-SHOW-THIS"
EMPTY = '{"facts": [], "insufficient_evidence_samples": [], "notes": ""}'


def token_of(page):
    match = re.search(r'name="token" value="([^"]+)"', page)
    return match.group(1) if match else None


class LiteratureWebTest(unittest.TestCase):
    def setUp(self):
        patches = [patch.dict(os.environ, {"GEONAMES_USERNAME": "", "AI_ENRICHMENT_ENABLED": "",
                                           "ANTHROPIC_API_KEY": "", "AI_ACCESS_CODE": "", "LOCAL_MODEL_URL": ""}),
                   patch("literature.config.ENV_FILE", ROOT / "missing.env"),
                   patch("align_sequences.shutil.which", return_value="mafft"),
                   patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(FIRST_SEQUENCE)),
                   patch.dict(app.config, {"LITERATURE_INLINE_JOBS": True})]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        budget_module.DAILY_LEDGER.update(date=None, spent=0.0)
        self.client = app.test_client()

    def run_example(self):
        with redirect_stdout(io.StringIO()):
            page = self.client.post("/", data={"run_literature_example": "1"}).get_data(as_text=True)
        return page, token_of(page)

    def full_demo(self):
        page, token = self.run_example()
        self.client.post("/literature/find", data={"token": token})
        response = self.client.post("/literature/analyze", data={"token": token, "confirm": "1"})
        page = self.client.get(response.headers["Location"]).get_data(as_text=True)
        return page, token

    def test_example_button_and_section_appear_without_any_model_call(self):
        home = self.client.get("/").get_data(as_text=True)
        self.assertIn("Run Literature Enrichment Example", home)
        self.assertIn("4 · Metadata enrichment (Sprint 4, optional)", home)
        with patch("literature.providers.MockProvider.run") as mock_run, \
                patch("literature.providers.ClaudeProvider.run") as claude_run:
            page, token = self.run_example()
        mock_run.assert_not_called()
        claude_run.assert_not_called()
        self.assertIsNotNone(token)
        self.assertIn("Sprint 4 — AI Metadata Enrichment", page)
        self.assertIn("Pairwise comparisons:", page)  # the Sprint 2 analysis still ran
        self.assertIn("SYNTHETIC demonstration", page)

    def test_find_analyze_review_apply(self):
        page, token = self.run_example()
        page = self.client.post("/literature/find", data={"token": token}).get_data(as_text=True)
        self.assertIn("Publications found:</strong> 11", page)
        self.assertIn("Exact accession mentioned in article text", page)
        page = self.client.post("/literature/analyze", data={"token": token}).get_data(as_text=True)
        self.assertIn("Before you start", page)
        self.assertIn("Run offline demonstration", page)
        self.assertNotIn("C. Review proposed metadata", page)  # nothing ran without confirmation
        response = self.client.post("/literature/analyze", data={"token": token, "confirm": "1"})
        self.assertEqual(response.status_code, 302)
        page = self.client.get(response.headers["Location"]).get_data(as_text=True)
        self.assertIn("C. Review proposed metadata", page)
        approvable = re.findall(r'name="review_(F\d+)" value="approve"', page)
        self.assertEqual(len(approvable), 7)
        self.assertNotIn('name="review_F15"', page)  # Tucson was lowered by the second check
        self.assertIn("literature_ai_proposals.csv", page)
        form = {"token": token, "review_F15": "approve", "review_F4": "approve",
                **{f"review_{fact}": "approve" for fact in approvable}}
        page = self.client.post("/literature/apply", data=form).get_data(as_text=True)
        self.assertIn("D. Original versus enriched metadata", page)
        self.assertIn("Newly meets all screening criteria", page)
        self.assertIn("not confirmed transmission", page)
        enriched = re.search(r'href="(/results/[^"]+/literature_enriched_metadata\.csv)"', page).group(1)
        table = self.client.get(enriched).get_data(as_text=True)
        self.assertIn("Original_Location_Name", table)
        self.assertNotIn("Tucson", table)
        self.assertNotIn("Albuquerque", table)
        report = self.client.get(re.search(r'href="(/results/[^"]+/literature_enrichment_report\.json)"',
                                           page).group(1)).get_data(as_text=True)
        self.assertIn("not calibrated probabilities", report)

    def test_article_text_is_escaped(self):
        page, _ = self.full_demo()
        self.assertNotIn("<title>", page.split("</head>", 1)[1])
        self.assertIn("Ignore all previous instructions", self.client.get(
            re.search(r'href="(/results/[^"]+/literature_evidence_passages\.csv)"', page).group(1)
        ).get_data(as_text=True))

    def test_expired_or_unknown_session(self):
        page = self.client.post("/literature/find", data={"token": "nope"}).get_data(as_text=True)
        self.assertIn("This enrichment session has expired", page)
        self.assertIn("This enrichment session has expired",
                      self.client.get("/literature/session/nope").get_data(as_text=True))

    def test_real_data_without_a_key_keeps_working_and_offers_no_paid_call(self):
        fake_source = MagicMock()
        fake_source.label = "fake live source"
        fake_source.synthetic = False
        fake_source.record_links.side_effect = lambda kind, value: (
            [("PMID:1", "direct", "record link", "fake")] if value == "SYNLIT_01" else [])
        fake_source.searchable.return_value = False
        fake_source.describe.return_value = {"PMID:1": {"Title": "T", "Source_URL": "javascript:alert(1)"}}
        fake_source.fetch.return_value = {"PMID:1": {"Passages": [
            {"Article_ID": "PMID:1", "Passage_ID": "PMID:1#p1", "Section": "METHODS", "Type": "paragraph",
             "Offset": 0, "Text": "Sample SYNLIT_01 was collected in Tempe, Arizona.", "Annotations": []}],
            "Text_Status": "Title and abstract only (via PubTator3)"}}
        files = {"files": [(io.BytesIO((EXAMPLE / "synthetic_literature_sequences.fasta").read_bytes()), "s.fasta"),
                           (io.BytesIO((EXAMPLE / "synthetic_literature_metadata.csv").read_bytes()), "m.csv")]}
        with redirect_stdout(io.StringIO()):
            review = self.client.post("/review", data=files, content_type="multipart/form-data").get_data(as_text=True)
            page = self.client.post("/run", data={"token": token_of(review)}).get_data(as_text=True)
        self.assertIn("Pairwise comparisons:", page)
        token = re.findall(r'name="token" value="([^"]+)"', page)[-1]
        with patch("literature_view.LiveSource", return_value=fake_source), \
                patch("literature.providers.ClaudeProvider.run") as claude_run:
            page = self.client.post("/literature/find", data={"token": token}).get_data(as_text=True)
            self.assertNotIn("javascript:", page)
            self.assertIn("AI analysis is disabled on this server", page)
            self.assertNotIn("Analyze Literature with Claude AI", page)
            page = self.client.post("/literature/analyze", data={"token": token, "confirm": "1"}).get_data(
                as_text=True)
        claude_run.assert_not_called()
        self.assertIn("AI analysis is disabled on this server", page)

    def test_live_run_needs_acknowledgement_access_and_stays_within_limits(self):
        fake_client = MagicMock()
        fake_client.messages.count_tokens.return_value = MagicMock(input_tokens=800)
        message = MagicMock(stop_reason="end_turn", model="claude-opus-5-5")
        message.content = [MagicMock(type="text", text=EMPTY)]
        message.usage = MagicMock(input_tokens=800, output_tokens=200, cache_creation_input_tokens=0,
                                  cache_read_input_tokens=0, iterations=None)
        fake_client.beta.messages.create.return_value = message
        page, token = self.run_example()
        session = literature_view.get_session(token)
        session["synthetic"] = False  # treat the demo data as real uploads for this test
        literature_view.find_literature(session, literature_view.LocalCorpusSource(literature_view.SYNTHETIC_CORPUS),
                                        max_articles=3)
        cache = tempfile.TemporaryDirectory()  # never the real result cache
        self.addCleanup(cache.cleanup)
        settings = {"AI_ENRICHMENT_ENABLED": "true", "ANTHROPIC_API_KEY": SECRET, "AI_ACCESS_CODE": "letmein",
                    "AI_MAX_REQUESTS_PER_RUN": "6", "AI_CACHE_DIR": cache.name}
        with patch.dict(os.environ, settings), patch.object(literature_view.ClaudeProvider, "client",
                                                            return_value=fake_client):
            page = self.client.post("/literature/analyze", data={"token": token}).get_data(as_text=True)
            self.assertIn("Worst-case cost", page)
            self.assertIn("Data sharing", page)
            self.assertIn("Access code", page)
            fake_client.beta.messages.create.assert_not_called()
            page = self.client.post("/literature/analyze", data={"token": token, "confirm": "1",
                                                                 "access_code": "letmein"}).get_data(as_text=True)
            self.assertIn("Tick the box", page)
            page = self.client.post("/literature/analyze", data={"token": token, "confirm": "1",
                                                                 "acknowledge": "on", "access_code": "wrong"}
                                    ).get_data(as_text=True)
            self.assertIn("access code", page)
            fake_client.beta.messages.create.assert_not_called()
            response = self.client.post("/literature/analyze", data={"token": token, "confirm": "1",
                                                                     "acknowledge": "on", "access_code": "letmein"})
            page = self.client.get(response.headers["Location"]).get_data(as_text=True)
        self.assertEqual(fake_client.beta.messages.create.call_count, 3)
        self.assertEqual(len(list(Path(cache.name).glob("*.json"))), 3)  # results cached for reuse
        self.assertIn("Requests sent:</strong> 3", page)
        self.assertNotIn(SECRET, page)
        for call in fake_client.beta.messages.create.call_args_list:
            self.assertNotIn("Sequence", call.kwargs["messages"][0]["content"])
            self.assertNotIn(FIRST_SEQUENCE, call.kwargs["messages"][0]["content"])

    def test_paid_runs_from_other_computers_need_an_access_code(self):
        with patch.dict(os.environ, {"AI_ACCESS_CODE": ""}):
            self.assertTrue(literature_view.live_request_allowed("127.0.0.1", None))
            self.assertFalse(literature_view.live_request_allowed("203.0.113.9", None))
        with patch.dict(os.environ, {"AI_ACCESS_CODE": "abc"}):
            self.assertFalse(literature_view.live_request_allowed("127.0.0.1", "wrong"))
            self.assertTrue(literature_view.live_request_allowed("203.0.113.9", "abc"))

    def test_diagnostics_show_ai_settings_without_calling_claude(self):
        mocked = {"status": "Operational", "message": "ok"}
        with patch.dict(STATUS_CACHE, {"checked_at": 0.0, "results": None}), \
                patch("app.check_mafft", return_value=mocked), patch("app.check_geonames", return_value=mocked), \
                patch("app.check_pubtator", return_value=mocked) as pubtator, \
                patch("literature_view.check_claude") as claude, \
                patch.dict(os.environ, {"AI_ENRICHMENT_ENABLED": "true", "ANTHROPIC_API_KEY": SECRET}):
            home = self.client.get("/").get_data(as_text=True)
            pubtator.assert_not_called()
            page = self.client.post("/system-status").get_data(as_text=True)
            claude.assert_not_called()
        self.assertIn("PubTator3: Operational", page)
        self.assertIn("Claude API configuration: Configured", page)
        self.assertIn("Check Claude API access", page)
        self.assertNotIn(SECRET, page + home)

    def test_explicit_claude_check_uses_token_counting_only(self):
        fake = MagicMock()
        fake.messages.count_tokens.return_value = MagicMock(input_tokens=5)
        with patch.dict(os.environ, {"AI_ENRICHMENT_ENABLED": "true", "ANTHROPIC_API_KEY": SECRET}), \
                patch("anthropic.Anthropic", return_value=fake):
            result = literature_view.check_claude("127.0.0.1", None)
        self.assertEqual(result["status"], "Operational")
        fake.messages.create.assert_not_called()
        with patch.dict(os.environ, {"AI_ENRICHMENT_ENABLED": ""}):
            self.assertEqual(literature_view.check_claude("127.0.0.1", None)["status"], "Disabled")


if __name__ == "__main__":
    unittest.main()
