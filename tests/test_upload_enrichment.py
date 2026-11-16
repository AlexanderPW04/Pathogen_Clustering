"""Sprint 4 started from the upload page and applied before the genomic analysis.

MAFFT is mocked (tests/fake_mafft.py), GeoNames is disabled, literature comes
from the synthetic article collections, and the extraction "model" is the
offline mock provider. No NCBI, PubTator, or Claude request is sent, so these
tests show the workflow, not that live services work.
"""

import io
import os
import re
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402

import literature_view  # noqa: E402
from app import STAGED, app, staged_item  # noqa: E402
from literature.discovery import STUDY, discover_publications  # noqa: E402
from literature.http_client import ServiceUnavailable  # noqa: E402
from literature.providers import MockProvider  # noqa: E402
from literature.sources import LocalCorpusSource  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
INTEGRATED = ROOT / "examples" / "integrated"
LITERATURE = ROOT / "examples" / "literature"
INT_CORPUS = INTEGRATED / "synthetic_integrated_articles.json"
INT_MOCK = INTEGRATED / "synthetic_integrated_mock_responses.json"
REFERENCE = (INTEGRATED / "synthetic_reference.fasta").read_text(encoding="utf-8").split("\n", 1)[1].replace("\n", "")
LIT_FIRST = (LITERATURE / "synthetic_literature_sequences.fasta").read_text(encoding="utf-8").split("\n")[1]
ALL_INTEGRATED = ["synthetic_individual_sequences.fasta", "synthetic_individual_metadata.csv",
                  "synthetic_reference.fasta", "synthetic_pooled_snvs.csv", "synthetic_pooled_base_counts.csv",
                  "synthetic_variant_signatures.csv"]


def upload(folder, name):
    return (io.BytesIO((folder / name).read_bytes()), name)


def token_of(page, after=None):
    if after:
        page = page[page.index(after):]
    match = re.search(r'name="token" value="([^"]+)"', page)
    return match.group(1) if match else None


def text(page):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", page))


class UploadEnrichmentTest(unittest.TestCase):
    reference = REFERENCE

    def setUp(self):
        patches = [patch.dict(os.environ, {"GEONAMES_USERNAME": "", "AI_ENRICHMENT_ENABLED": "",
                                           "ANTHROPIC_API_KEY": ""}),
                   patch("literature.config.ENV_FILE", ROOT / "missing.env"),
                   patch("align_sequences.shutil.which", return_value="mafft"),
                   patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(self.reference)),
                   patch.dict(app.config, {"LITERATURE_INLINE_JOBS": True})]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        self.addCleanup(STAGED.clear)
        self.client = app.test_client()

    # helpers ---------------------------------------------------------------------------
    def review(self, files, **form):
        response = self.client.post("/review", data={"files": files, **form}, content_type="multipart/form-data")
        page = response.get_data(as_text=True)
        return page, token_of(page)

    def post(self, url, data):
        with redirect_stdout(io.StringIO()):
            return self.client.post(url, data=data).get_data(as_text=True)

    def enrichment_example(self):
        page = self.post("/", {"run_integrated_example": "1", "enrich_example": "on"})
        return page, token_of(page, 'id="literature"')

    def analyze_all(self, token):
        self.client.post("/literature/find", data={"token": token})
        response = self.client.post("/literature/analyze", data={"token": token, "confirm": "1"})
        page = self.client.get(response.headers["Location"]).get_data(as_text=True)
        return page, re.findall(r'name="review_(F\d+)" value="approve"', page)

    def download(self, page, name):
        link = re.search(rf'href="(/results/[^"]+/{re.escape(name)})"', page).group(1)
        return pd.read_csv(io.StringIO(self.client.get(link).get_data(as_text=True)))

    # Toggle OFF ---------------------------------------------------------------------------
    def test_toggle_off_keeps_the_normal_workflow(self):
        home = self.client.get("/").get_data(as_text=True)
        self.assertIn("Use PubTator + Claude AI to find missing or more precise sample information", home)
        self.assertRegex(home, r'name="enrich" id="enrich" >')  # unchecked by default
        page, token = self.review([upload(INTEGRATED, name) for name in ALL_INTEGRATED])
        self.assertIn(">Run analysis</button>", page)
        self.assertNotIn("Find and Check Publications", page)
        self.assertFalse(staged_item(token)["enrich"])
        result = self.post("/run", {"token": token})
        self.assertNotIn('id="enrichment-report"', result)
        self.assertNotIn('id="literature"', result)

    # Toggle ON: review preview ----------------------------------------------------------
    def test_toggle_on_shows_gaps_on_review_without_any_external_request(self):
        with patch("literature_view.LiveSource") as live:
            page, token = self.review([upload(INTEGRATED, name) for name in ALL_INTEGRATED], enrich="on",
                                      publication_ids="PMID 123, SYNPUB-INT-2")
        live.assert_not_called()
        self.assertTrue(staged_item(token)["enrich"])
        self.assertIn("10 of 14 samples</strong> have collection details", page)
        self.assertIn("State-only location: 1", page)
        self.assertIn("Find and Check Publications", page)
        self.assertIn("Skip enrichment and run analysis", page)
        self.assertIn("Nothing has been sent anywhere yet", page)

    # Combined (the synthetic integrated example) ---------------------------------------
    def test_integrated_example_with_enrichment_end_to_end(self):
        page, token = self.enrichment_example()
        self.assertIn("SYNTHETIC offline demonstration", page)
        self.assertIn("Metadata that could be improved", page)
        found = self.client.post("/literature/find", data={"token": token}).get_data(as_text=True)
        self.assertIn("Direct sample links:</strong> 1", found)
        self.assertIn("Study-level links:</strong> 8", found)
        page, approvable = self.analyze_all(token)
        self.assertEqual(len(approvable), 9)
        self.assertIn('action="/enrich/run"', page)
        self.assertIn("Unsupported", page)  # the planted fabricated quotation for WW_S1
        result = self.post("/enrich/run", {"token": token, **{f"review_{fact}": "approve" for fact in approvable}})
        self.assertIn('id="enrichment-report"', result)
        self.assertIn("Approved:</strong> 9", result)
        self.assertIn("Sprint 4 metadata enrichment was enabled: 9 of 9 proposed improvement(s) approved", result)
        self.assertIn("Newly meets all screening criteria", result)
        comparison = self.download(result, "cross_source_comparison.csv").set_index("Variant_Key")
        self.assertAlmostEqual(comparison.loc["SYNREF.1:15:T>C", "Nearest_Km_Between_Carriers"], 2.495, places=2)
        self.assertEqual(int((comparison["Category"] == "Shared: observed in both sources").sum()), 2)
        enriched = self.download(result, "literature_enriched_metadata.csv").set_index("Sample_ID")
        self.assertEqual(enriched.loc["IND_D", "Original_Location_Name"], "USA: AZ")
        self.assertEqual(enriched.loc["IND_D", "Location_Name"], "Tempe, Arizona")
        individual = self.download(result, "individual_samples.csv").set_index("Sample_ID")
        self.assertEqual(individual.loc["IND_D", "Location_Name"], "Tempe, Arizona")
        pooled = self.download(result, "pooled_samples_with_qc.csv").set_index("Sample_ID")
        self.assertAlmostEqual(pooled.loc["WW_N1", "Latitude"], 33.415)
        self.assertTrue(pd.isna(pooled.loc["WW_X1", "Latitude"]))  # no evidence, unchanged
        self.assertIn("Revisit metadata enrichment", result)
        # The gaps listed at review are reported again: resolved, improved, or unresolved.
        gaps = text(result[result.index("Metadata gaps after enrichment"):result.index("Values used in this analysis")])
        self.assertIn("Collection location 10 8 1 1 Improved: IND_D (still: named place without coordinates) "
                      "Unresolved: WW_X1", gaps)

    def test_genetic_outputs_are_unchanged_by_enrichment(self):
        plain = self.post("/", {"run_integrated_example": "1"})
        page, token = self.enrichment_example()
        _, approvable = self.analyze_all(token)
        enriched = self.post("/enrich/run", {"token": token, **{f"review_{fact}": "approve" for fact in approvable}})
        before = self.download(plain, "pairwise_distances.csv")
        after = self.download(enriched, "pairwise_distances.csv")
        pd.testing.assert_frame_equal(before[["Sample_ID_1", "Sample_ID_2", "P_Distance", "Comparable_Positions"]],
                                      after[["Sample_ID_1", "Sample_ID_2", "P_Distance", "Comparable_Positions"]])
        self.assertEqual(set(self.download(plain, "variant_positions.csv")["Variant_Key"]),
                         set(self.download(enriched, "variant_positions.csv")["Variant_Key"]))

    def test_rejected_values_are_not_applied(self):
        page, token = self.enrichment_example()
        _, approvable = self.analyze_all(token)
        choices = {f"review_{fact}": "approve" for fact in approvable}
        choices[f"review_{approvable[0]}"] = "reject"  # IND_D's location
        result = self.post("/enrich/run", {"token": token, **choices})
        self.assertIn("Rejected or left unresolved:</strong> 1", result)
        individual = self.download(result, "individual_samples.csv").set_index("Sample_ID")
        self.assertEqual(individual.loc["IND_D", "Location_Name"], "USA: AZ")

    def test_skipping_runs_with_the_original_metadata_and_deletes_uploads(self):
        page, token = self.enrichment_example()
        staged = literature_view.get_session(token)["staged_token"]
        result = self.post("/enrich/run", {"token": token, "skip": "1"})
        self.assertIn("Variants observed in both sources:</strong> 2", result)
        self.assertNotIn('id="enrichment-report"', result)
        self.assertIsNone(staged_item(staged))
        individual = self.download(result, "individual_samples.csv").set_index("Sample_ID")
        self.assertEqual(individual.loc["IND_D", "Location_Name"], "USA: AZ")

    # Uploaded individual or pooled data (mock provider, synthetic corpora) -------------------
    def run_uploaded(self, files, corpus, mock, **form):
        page, staged = self.review(files, enrich="on", **form)
        with patch("literature_view.literature_source", return_value=LocalCorpusSource(corpus)), \
                patch("literature_view.choose_provider", return_value=(MockProvider(mock), "Mock (test)")):
            page = self.post("/enrich/start", {"token": staged, "enrich": "on", "enrich_shown": "1"})
            token = token_of(page, 'id="literature"')
            found = self.client.post("/literature/find", data={"token": token}).get_data(as_text=True)
            page, approvable = self.analyze_all(token)
        result = self.post("/enrich/run", {"token": token, **{f"review_{fact}": "approve" for fact in approvable}})
        return found, page, result

    def test_toggle_on_with_pooled_data(self):
        found, _, result = self.run_uploaded([upload(INTEGRATED, "synthetic_pooled_snvs.csv")], INT_CORPUS, INT_MOCK)
        self.assertIn("Study-level links:</strong> 8", found)
        self.assertIn("Values changed:</strong> 8", result)
        pooled = self.download(result, "pooled_samples_with_qc.csv").set_index("Sample_ID")
        self.assertAlmostEqual(pooled.loc["WW_S2", "Longitude"], -111.94)
        self.assertEqual(pooled.loc["WW_S2", "Coordinate_Source"], "Publication (approved literature evidence)")

    def test_toggle_on_with_individual_data(self):
        files = [upload(LITERATURE, "synthetic_literature_sequences.fasta"),
                 upload(LITERATURE, "synthetic_literature_metadata.csv")]
        with patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(LIT_FIRST)):
            found, page, result = self.run_uploaded(files, LITERATURE / "synthetic_articles.json",
                                                    LITERATURE / "mock_responses.json", intended_mode="individual")
        # Links come from identifiers already in the files (accessions, Publication_ID, BioProject).
        self.assertIn("Direct sample links:</strong> 10", found)
        self.assertIn("Study-level links:</strong> 1", found)
        # Every linked article is retrieved (free); the paid limit of 5 articles per run splits
        # the Claude analysis into batches, each approved separately.
        self.assertNotIn("not retrieved", found)
        self.assertIn("batches of at most 5 (the per-run limit)", found)
        self.assertIn('id="enrichment-report"', result)
        self.assertIn("Newly meets all screening criteria", result)  # SYNLIT_01 – SYNLIT_09
        enriched = self.download(result, "literature_enriched_metadata.csv").set_index("Sample_ID")
        self.assertEqual(enriched.loc["SYNLIT_09", "Location_Name"], "Tempe, Arizona, USA")
        self.assertEqual(enriched.loc["SYNLIT_11", "Location_Name"], "USA: Arizona")  # contradiction kept out
        # The Sprint 2 run itself used the approved date (2026-03-05): the pair is 2 days apart.
        name = re.search(r'href="/results/[^"]+/([^"]+_pairwise_comparisons_temporal_geographic\.csv)"', result).group(1)
        pairs = self.download(result, name).set_index(["Sample_ID_1", "Sample_ID_2"])
        self.assertEqual(pairs.loc[("SYNLIT_01", "SYNLIT_09"), "Days_Apart"], 2)
        self.assertIn("Unsupported", page)  # fabricated quotation rejected

    def test_manual_dataset_identifier_links_the_study(self):
        samples = pd.DataFrame({"Sample_ID": ["WW_X1", "WW_N1"], "Location_Name": ["Site", "Site"]})
        source = LocalCorpusSource(INT_CORPUS)
        table, _, notes, _ = discover_publications(samples, source, dataset_identifiers=[("Publication", "SYNPUB-INT-1")])
        rows = table[table["Article_ID"] == "SYNPUB-INT-1"]
        self.assertEqual(set(rows["Sample_ID"]), {"WW_X1", "WW_N1"})
        self.assertTrue((rows["Relationship_Type"] == STUDY).all())
        self.assertTrue(rows["Retrieval_Method"].str.contains("entered by the user").all())
        _, _, notes, _ = discover_publications(samples, source, dataset_identifiers=[("Publication", "NOPE")])
        self.assertIn("The dataset identifier(s) entered did not resolve to a publication.", notes)

    # No links, unavailable services, Claude disabled ------------------------------------------
    def test_no_publication_found_still_continues(self):
        page, staged = self.review([upload(INTEGRATED, "synthetic_individual_sequences.fasta"),
                                    upload(INTEGRATED, "synthetic_individual_metadata.csv")], enrich="on")
        with patch("literature_view.literature_source", return_value=LocalCorpusSource(LITERATURE / "synthetic_articles.json")):
            page = self.post("/enrich/start", {"token": staged})
            token = token_of(page, 'id="literature"')
            found = self.client.post("/literature/find", data={"token": token}).get_data(as_text=True)
        self.assertIn("No publication was found", found)
        self.assertIn("Continue to Genomic Analysis", found)
        result = self.post("/enrich/run", {"token": token, "skip": "1"})
        self.assertIn("Pairwise comparisons:", result)

    def test_pubtator_unavailable_and_claude_disabled(self):
        page, staged = self.review([upload(INTEGRATED, "synthetic_individual_sequences.fasta"),
                                    upload(INTEGRATED, "synthetic_individual_metadata.csv")], enrich="on")
        failing = LocalCorpusSource(INT_CORPUS)
        with patch.object(failing, "record_links", side_effect=ServiceUnavailable("PubTator3 could not be reached.")), \
                patch("literature_view.literature_source", return_value=failing):
            page = self.post("/enrich/start", {"token": staged})
            token = token_of(page, 'id="literature"')
            found = self.client.post("/literature/find", data={"token": token}).get_data(as_text=True)
        self.assertIn("PubTator3 could not be reached.", found)
        self.assertIn("Continue to Genomic Analysis", found)
        # With a link found but Claude disabled, analysis is refused and the skip remains.
        page, staged = self.review([upload(INTEGRATED, "synthetic_individual_sequences.fasta"),
                                    upload(INTEGRATED, "synthetic_individual_metadata.csv")], enrich="on")
        with patch("literature_view.literature_source", return_value=LocalCorpusSource(INT_CORPUS)), \
                patch("literature.providers.ClaudeProvider.run") as claude:
            page = self.post("/enrich/start", {"token": staged})
            token = token_of(page, 'id="literature"')
            found = self.client.post("/literature/find", data={"token": token}).get_data(as_text=True)
            refused = self.client.post("/literature/analyze", data={"token": token, "confirm": "1"}).get_data(as_text=True)
        claude.assert_not_called()
        self.assertIn("AI analysis is disabled on this server", found)
        self.assertIn("AI analysis is disabled on this server", refused)
        self.assertIn("The genomic analysis can run without it.", found)

    # State, navigation ---------------------------------------------------------------------------
    def test_sessions_expire_are_isolated_and_need_the_staged_files(self):
        page, token = self.enrichment_example()
        self.assertIn("This enrichment session has expired", self.post("/enrich/run", {"token": "guess"}))
        session = literature_view.get_session(token)
        with patch("literature_view.monotonic", return_value=10 ** 9):
            self.assertIsNone(literature_view.get_session(token))
        STAGED.pop(session["staged_token"])  # the uploads expired
        page = self.post("/enrich/run", {"token": token, "skip": "1"})
        self.assertIn("The uploaded files are no longer on the server", page)

    def test_enrichment_page_navigation_returns_to_the_review(self):
        page, token = self.enrichment_example()
        # Phase 4A renumbered the steps (3 is now "Review available analyses").
        self.assertIn("4 · Run · Metadata enrichment (Sprint 4, optional)", page)
        staged = literature_view.get_session(token)["staged_token"]
        names = re.findall(r'name="previous_file" value="([^"]+)"', page)
        back = self.client.post("/back", data={"target": "review", "token": staged, "previous_file": names})
        self.assertIn("Step 2 · Review your dataset", back.get_data(as_text=True))


if __name__ == "__main__":
    unittest.main()
