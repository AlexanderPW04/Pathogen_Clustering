"""Usability tests: navigation, analysis modes, progressive results, and status text.

MAFFT is mocked with tests/fake_mafft.py and GeoNames is disabled. External
service checks are replaced with fixed answers, so these tests check the
page behavior, not whether live services work.
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

from app import STAGED, STATUS_CACHE, app, staged_item  # noqa: E402
from integrated_view import table_view  # noqa: E402

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "integrated"
REFERENCE = (EXAMPLE / "synthetic_reference.fasta").read_text(encoding="utf-8").split("\n", 1)[1].replace("\n", "")
INDIVIDUAL = ["synthetic_individual_sequences.fasta", "synthetic_individual_metadata.csv"]
ALL_FILES = INDIVIDUAL + ["synthetic_reference.fasta", "synthetic_pooled_snvs.csv",
                          "synthetic_pooled_base_counts.csv", "synthetic_variant_signatures.csv"]


def files(names):
    return [(io.BytesIO((EXAMPLE / name).read_bytes()), name) for name in names]


def token_of(page):
    match = re.search(r'name="token" value="([^"]+)"', page)
    return match.group(1) if match else None


class UiNavigationTest(unittest.TestCase):
    def setUp(self):
        patches = [patch.dict(os.environ, {"GEONAMES_USERNAME": ""}),
                   patch("align_sequences.shutil.which", return_value="mafft"),
                   patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(REFERENCE))]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        self.client = app.test_client()

    def review(self, names, mode="auto"):
        response = self.client.post("/review", data={"files": files(names), "intended_mode": mode},
                                    content_type="multipart/form-data")
        page = response.get_data(as_text=True)
        return page, token_of(page)

    def run_reviewed(self, token, **form):
        with redirect_stdout(io.StringIO()):
            return self.client.post("/run", data={"token": token, **form}).get_data(as_text=True)

    def integrated_example(self):
        with redirect_stdout(io.StringIO()):
            return self.client.post("/", data={"run_integrated_example": "1"}).get_data(as_text=True)

    # Home page -----------------------------------------------------------------
    def test_home_explains_the_three_sprints_and_offers_modes(self):
        page = self.client.get("/").get_data(as_text=True)
        for text in ["How this works", "Sprint 2 · Individual sequences", "Sprint 3 · Pooled sequencing",
                     "Shared reference genome = the connection", "Sprint 4 · Optional AI metadata enrichment",
                     "done before upload, offline", "Processed SNV table <em>(you upload this)</em>",
                     "Detect from my files", "Individual sequences (Sprint 2)",
                     "Pooled sequencing results (Sprint 3)", "Combined (Sprint 2 + Sprint 3)",
                     "What optional information unlocks", "Per-position A/C/G/T counts",
                     "MAFFT is not used for pooled reads", "Genetic p-distance", "Unresolved evidence"]:
            self.assertIn(text, page)
        self.assertRegex(page, r'value="auto" checked')
        # Example facts come from the bundled files, not from fixed text.
        self.assertIn("5 individual sequence records, 9 pooled samples from 2 sites, and a 120-base", page)
        # Legacy pages stay reachable but are not prominent.
        self.assertIn("<summary>Legacy single-source pages</summary>", page)

    def test_example_files_download_and_unknown_names_404(self):
        response = self.client.get("/examples/synthetic_pooled_snvs.csv")
        self.assertEqual(response.status_code, 200)
        self.assertIn("attachment", response.headers["Content-Disposition"])
        self.assertIn(b"Sample_ID", response.data)
        self.assertEqual(self.client.get("/examples/app.py").status_code, 404)
        self.assertEqual(self.client.get("/examples/..%2Fapp.py").status_code, 404)

    # Mode detection ------------------------------------------------------------
    def test_detected_mode_and_missing_files_for_the_chosen_mode(self):
        page, _ = self.review(INDIVIDUAL, "combined")
        self.assertIn("Detected analysis:</strong> Individual sequences (Sprint 2)", page)
        self.assertIn("You chose:</strong> Combined (Sprint 2 + Sprint 3)", page)
        self.assertIn("reference genome FASTA, processed pooled SNV table", page)
        self.assertIn('name="intended_mode" value="combined"', page)  # kept for Update review / Run

    def test_combined_and_pooled_files_are_recognized_without_choosing_a_mode(self):
        page, _ = self.review(ALL_FILES)
        self.assertIn("Detected analysis:</strong> Combined (Sprint 2 + Sprint 3)", page)
        self.assertNotIn("You chose:", page)
        page, _ = self.review(["synthetic_pooled_snvs.csv"], "pooled")
        self.assertIn("Detected analysis:</strong> Pooled sequencing results (Sprint 3)", page)
        self.assertIn("the files support it", page)
        page, _ = self.review(INDIVIDUAL + ["synthetic_pooled_snvs.csv"])
        self.assertIn("add a shared reference FASTA to connect them", page)

    # Navigation ----------------------------------------------------------------
    def test_review_page_has_navigation_and_back_to_upload_discards_staged_files(self):
        page, token = self.review(INDIVIDUAL)
        self.assertIn('aria-label="Page navigation"', page)
        self.assertIn("← Back to upload", page)
        self.assertIn("Start new analysis", page)
        self.assertIsNotNone(staged_item(token))
        back = self.client.post("/back", data={"target": "upload", "token": token,
                                               "previous_file": INDIVIDUAL, "max_p_distance": "0.02"})
        page = back.get_data(as_text=True)
        self.assertEqual(back.status_code, 200)
        self.assertIsNone(staged_item(token))
        self.assertIn("were removed from the server", page)
        self.assertIn("synthetic_individual_sequences.fasta, synthetic_individual_metadata.csv", page)
        self.assertIn('name="max_p_distance" min="0" step="any" value="0.02"', page)  # settings kept

    def test_back_to_review_works_only_while_files_still_exist(self):
        page, token = self.review(INDIVIDUAL)
        back = self.client.post("/back", data={"target": "review", "token": token}).get_data(as_text=True)
        self.assertIn("Step 2 · Review your dataset", back)
        result = self.run_reviewed(token)
        self.assertIn('href="#overview">Results overview</a>', result)
        self.assertIn("← Back to file review", result)
        self.assertIsNone(staged_item(token))  # deleted after the run, as before
        names = re.findall(r'name="previous_file" value="([^"]+)"', result)
        self.assertEqual(sorted(set(names)), sorted(INDIVIDUAL))
        back = self.client.post("/back", data={"target": "review", "previous_file": names}).get_data(as_text=True)
        self.assertIn("Your uploaded files are no longer on the server", back)
        self.assertIn("Files you used before (select them again)", back)
        self.assertNotIn("Step 2 · Review your dataset", back)
        self.assertNotIn('role="alert"', back)

    def test_back_after_an_example_explains_there_is_nothing_to_review(self):
        page = self.client.post("/back", data={"target": "review", "example": "1",
                                               "previous_file": ["synthetic_reference.fasta"]}).get_data(as_text=True)
        self.assertIn("The example used bundled files", page)
        self.assertNotIn("Files you used before", page)

    # Results -------------------------------------------------------------------
    def test_integrated_example_overview_story_and_progressive_sections(self):
        page = self.integrated_example()
        self.assertIn("Combined analysis (Sprint 2 + Sprint 3) on one verified reference", page)
        for finding in ["5 individual sequence records were analyzed; 4 could be placed on the reference genome.",
                        "9 pooled sample records were analyzed; 8 were kept after quality control.",
                        "2 variant(s) were observed in both sources",
                        "Population structure: Paper method: Yue-Clayton PCoA"]:
            self.assertIn(finding, page)
        self.assertIn("The invented inputs:</strong> 5 individual sequence records (Sprint 2), 9 pooled sample records", page)
        story = page[page.index('id="shared-variant-example"'):]
        story = story[:story.index("</div>")]
        for text in ["SYNREF.1:15:T&gt;C", "Position 15 of reference SYNREF.1 is <code>T</code>",
                     "IND_A (2026-03-02, Tempe, Arizona)", "WW_N1 (2026-03-01, North Synthetic Sewershed)",
                     "0.05 to 0.7 (median 0.175)", "Shared: observed in both sources",
                     "pooled 1 day(s) earlier", "no kilometer distance can be computed"]:
            self.assertIn(text, story)
        sections = re.findall(r'<section id="([a-z-]+)" class="result-section">', page)
        self.assertEqual(sections, ["overview", "cross-source", "variation", "temporal", "geography",
                                    "quality", "population", "relationships"])
        self.assertEqual(page.count('<p class="sprint-tag">'), 8)
        self.assertEqual(page.count("↑ Back to results overview"), 8)  # 7 sections + provenance
        figures = page.count('<figure>')
        self.assertEqual(page.count("How to read this chart"), figures)
        self.assertIn("In these results:</strong> 2 of 8 variants were observed in both sources", page)
        self.assertIn("More figures in this section", page)
        self.assertIn("Detailed tables", page)
        self.assertIn('<span class="badge origin">Faithful to Fontenele et al. (2021) Fig 2A</span>', page)

    def test_pooled_only_results_label_the_mode_and_exploratory_pcoa(self):
        page, token = self.review(["synthetic_pooled_snvs.csv"])
        result = self.run_reviewed(token)
        self.assertIn("Pooled sequencing analysis (Sprint 3)", result)
        self.assertNotIn('id="shared-variant-example"', result)
        self.assertIn("Cross-source comparison needs individual sequences", result)
        self.assertIn("Exploratory alternative (not the paper", result)

    def test_individual_only_results_label_the_mode(self):
        page, token = self.review(INDIVIDUAL)
        result = self.run_reviewed(token)
        self.assertIn("Individual sequence analysis (Sprint 2)", result)
        self.assertIn("Sprint 2 screening:", result)

    def test_missing_optional_metadata_is_not_shown_as_zero(self):
        metadata = b"Accession,Pathogen\nIND_A,V\nIND_B,V\nIND_C,V\nIND_D,V\nIND_E,V\n"
        response = self.client.post("/review", data={"files": files(INDIVIDUAL[:1]) + [(io.BytesIO(metadata), "meta.csv")]},
                                    content_type="multipart/form-data")
        result = self.run_reviewed(token_of(response.get_data(as_text=True)))
        self.assertIn("Collection dates:</strong> Not available", result)
        self.assertIn("Pairwise comparisons:", result)

    def test_long_tables_show_a_few_rows_and_point_to_full_downloads(self):
        table = pd.DataFrame({"Value": [index / 3 for index in range(300)]})
        view = table_view(table, limit=25, expand=200, title="Long", filename="long.csv")
        self.assertEqual((len(view["rows"]), len(view["more_rows"]), view["total"], view["shown"]), (25, 200, 300, 225))
        page = self.integrated_example()
        link = re.search(r'href="(/results/[^"]+/pooled_variants_with_qc\.csv)"', page).group(1)
        downloaded = pd.read_csv(io.StringIO(self.client.get(link).get_data(as_text=True)))
        self.assertEqual(len(downloaded), 22)  # every pooled record, not the page sample

    def test_sprint4_is_optional_and_has_no_article_summaries(self):
        page = self.integrated_example()
        # Enrichment was not chosen: no empty Sprint 4 section, only a link to revisit it.
        self.assertNotIn('id="literature"', page)
        self.assertIn("Open metadata enrichment for these results", page)
        self.assertNotIn("Surveillance context", page)
        self.assertNotIn("literature_surveillance_context.csv", page)

    # External dependencies -----------------------------------------------------
    def test_status_explains_each_service_and_failures(self):
        failed = {"status": "Check failed", "message": "Not found."}
        working = {"status": "Operational", "message": "ok"}
        with patch.dict(STATUS_CACHE, {"checked_at": 0.0, "results": None}), \
                patch("app.check_mafft", return_value=failed), patch("app.check_geonames", return_value=failed), \
                patch("app.check_pubtator", return_value=working), \
                patch.dict(os.environ, {"AI_ENRICHMENT_ENABLED": ""}), \
                patch("literature.config.ENV_FILE", Path("missing.env")):
            home = self.client.get("/").get_data(as_text=True)
            page = self.client.post("/system-status").get_data(as_text=True)
        self.assertIn("Not checked yet", home)
        self.assertIn("MAFFT: Check failed", page)
        self.assertIn("Pooled sequencing analysis (Sprint 3) still works", page)
        self.assertIn("every analysis that does not need geocoding still runs", page)
        self.assertIn("Claude API configuration: Disabled", page)
        self.assertIn("every genomic analysis and the rest of Sprint 4 still work", page)

    def test_pooled_analysis_runs_when_mafft_is_missing(self):
        with patch("align_sequences.shutil.which", return_value=None):
            page, token = self.review(["synthetic_pooled_snvs.csv"])
            result = self.run_reviewed(token)
        self.assertIn("Pooled samples:</strong> 9", result)
        self.assertNotIn("Analysis could not be completed", result)

    def tearDown(self):
        STAGED.clear()


if __name__ == "__main__":
    unittest.main()
