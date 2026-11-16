"""Results dashboard: the page shows what actually ran (Sprint 4 UI work, results page).

SYNTHETIC data only (the bundled integrated example files, a synthetic SARS-CoV-2-shaped SNV table, and
mocked MAFFT). Each configuration is run once; the tests check which categories, figures, and downloads
the page shows, that unsupported analyses are explained briefly instead of drawn as empty panels, that
every figure stays downloadable, and that each displayed count matches its own denominator.
"""

import io
import json
import os
import re
import sys
import unittest
from contextlib import redirect_stdout
from html import unescape
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402

from app import app  # noqa: E402
from genomics.references import load_reference_fasta  # noqa: E402
from genomics.report_exports import report_figure_names  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "integrated"
ANNOTATION = ROOT / "examples" / "annotation"
REFERENCE = (EXAMPLE / "synthetic_reference.fasta").read_text(encoding="utf-8").split("\n", 1)[1].replace("\n", "")
SARS2 = load_reference_fasta(ANNOTATION / "NC_045512.2.fasta")["Sequence"]
SITES = [(241, "T"), (3037, "T"), (14408, "T"), (23403, "G"), (28881, "A")]


def example(name, drop=()):
    data = (EXAMPLE / name).read_bytes()
    if drop:
        data = pd.read_csv(io.BytesIO(data)).drop(columns=list(drop)).to_csv(index=False).encode()
    return io.BytesIO(data), name


def sars2_counts_table():
    """SYNTHETIC: three dated samples; one identical repeated record and one record below the 20x depth limit."""

    rows = []
    for number, (date, depth) in enumerate([("2026-01-05", 200), ("2026-01-12", 400), ("2026-01-19", 300)]):
        for index, (position, alt) in enumerate(SITES):
            alt_depth = depth // (2 + index % 3) - number * 5
            rows.append({"Sample_ID": f"WW{number + 1}", "Reference_Genome": "NC_045512.2", "Position": position,
                         "Ref": SARS2[position - 1], "Alt": alt, "Alt_Frequency": round(alt_depth / depth, 6),
                         "Read_Depth": depth, "Alt_Depth": alt_depth, "Collection_Date": date,
                         "Location_Name": "Synthetic Site"})
    rows += [dict(rows[0]), dict(rows[0])]  # two identical repeats: counted once
    rows.append({**rows[1], "Position": 5000, "Ref": SARS2[4999], "Alt": "A" if SARS2[4999] != "A" else "C",
                 "Read_Depth": 12, "Alt_Depth": 6, "Alt_Frequency": 0.5})  # fails the paper's 20x QC
    return pd.DataFrame(rows).to_csv(index=False).encode()


def text_of(page):
    return unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", page)))


def section(page, identifier):
    match = re.search(rf'<section id="{identifier}" class="result-section[^"]*"[^>]*>(.*?)</section>', page, re.S)
    return match.group(1) if match else None


def category_ids(page):
    return re.findall(r'<section id="([a-z-]+)" class="result-section', page)


def images(page):
    return re.findall(r'<img src="data:image/png;base64,[^"]+" alt="([^".]+)', page)


class ResultsDashboardTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.patches = [patch.dict(os.environ, {"GEONAMES_USERNAME": "", "AI_ENRICHMENT_ENABLED": "false"}),
                       patch("align_sequences.shutil.which", return_value="mafft"),
                       patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(REFERENCE))]
        for item in cls.patches:
            item.start()
        cls.client = app.test_client()

        def post(data):
            with redirect_stdout(io.StringIO()):
                return cls.client.post("/", data=data, content_type="multipart/form-data").get_data(as_text=True)

        individual = {"fasta_file": example("synthetic_individual_sequences.fasta"),
                      "metadata_file": example("synthetic_individual_metadata.csv")}
        cls.pages = {
            "pooled": post({"snv_file": example("synthetic_pooled_snvs.csv")}),
            "individual": post(individual),
            "pooled_no_dates": post({"snv_file": example("synthetic_pooled_snvs.csv",
                                                         drop=("Collection_Date", "Location_Name"))}),
            "annotated": post({"snv_file": (io.BytesIO(sars2_counts_table()), "sars2_counts.csv"),
                               "reference_file": (io.BytesIO((ANNOTATION / "NC_045512.2.fasta").read_bytes()),
                                                  "NC_045512.2.fasta"),
                               "annotation_file": (io.BytesIO((ANNOTATION / "NC_045512.2.gb").read_bytes()),
                                                   "NC_045512.2.gb")}),
        }
        with redirect_stdout(io.StringIO()):
            cls.pages["combined"] = cls.client.post("/", data={"run_integrated_example": "1"}).get_data(as_text=True)

    @classmethod
    def tearDownClass(cls):
        for item in cls.patches:
            item.stop()

    def download(self, page, name):
        link = re.search(rf'href="(/results/[^"]+/{re.escape(name)})"', page)
        self.assertIsNotNone(link, f"{name} is not offered for download")
        response = self.client.get(link.group(1))
        self.assertEqual(response.status_code, 200, name)
        return response.data

    # --- which categories appear for each analysis -------------------------------------------------

    def test_pooled_only_shows_pooled_categories_and_the_proper_title(self):
        page = self.pages["pooled"]
        self.assertIn("Pooled Pathogen Sequencing — Sprint 3 + 5", page)
        self.assertNotIn("Sprint 2 + Sprint 3", page)
        self.assertEqual(category_ids(page)[0], "overview")
        self.assertNotIn("relationships", category_ids(page))
        self.assertNotIn("cross-source", category_ids(page))
        self.assertEqual(category_ids(page)[-1], "downloads")
        self.assertNotIn("Individual Sequence Relationships", page)
        self.assertNotIn("Cross-Source Variant Evidence", page)
        self.assertNotIn("Variants carried by individual sequences", page)
        # One short line instead of empty Sprint 2 sections.
        self.assertIn("Individual-sequence and cross-source analyses (Sprint 2)", text_of(page))
        self.assertIn("Not part of this analysis", page)

    def test_individual_only_hides_pooled_and_sprint5_panels(self):
        page = self.pages["individual"]
        self.assertIn("Individual Pathogen Sequences — Sprint 2", page)
        self.assertIn("relationships", category_ids(page))
        for absent in ["cross-source", "variation", "advanced"]:
            self.assertNotIn(absent, category_ids(page))
        for panel in ['id="sequencing-uncertainty"', 'id="mutation-annotation"', 'id="mutation-frequency-groups"']:
            self.assertNotIn(panel, page)
        self.assertIn("Potential clusters:</strong>", section(page, "relationships"))
        self.assertIn("Pooled-sequencing and cross-source analyses (Sprint 3 + 5)", text_of(page))

    def test_combined_puts_cross_source_evidence_first_after_the_overview(self):
        page = self.pages["combined"]
        self.assertIn("Combined Analysis — Sprint 2 + 3 + 5", page)
        ids = category_ids(page)
        self.assertEqual(ids[:3], ["overview", "cross-source", "relationships"])
        self.assertIn('id="shared-variant-example"', section(page, "cross-source"))
        self.assertIn("advanced", ids)
        self.assertIn("Paper method: Yue-Clayton PCoA", section(page, "advanced"))

    def test_navigation_links_match_the_rendered_categories(self):
        for name, page in self.pages.items():
            nav = re.search(r'<nav class="sections"[^>]*>(.*?)</nav>', page, re.S).group(1)
            links = re.findall(r'href="#([a-z-]+)"', nav)
            self.assertEqual(links, category_ids(page), name)
            self.assertLessEqual(len(links), 8, name)

    # --- unsupported analyses: short explanations, not empty charts ------------------------------------

    def test_missing_coverage_replaces_the_empty_coverage_chart_with_a_sentence(self):
        page = self.pages["pooled"]
        self.assertIn("Genome-wide coverage was not supplied, so coverage percentages could not be calculated.",
                      section(page, "quality"))
        self.assertNotIn("Sequencing quality by source", images(page))
        self.assertIn("Read-depth distribution", images(page))  # depth at the reported variants stays
        self.assertTrue(self.download(page, "quality.png").startswith(b"\x89PNG"))  # still downloadable

    def test_named_places_without_coordinates_get_a_summary_not_a_blank_map(self):
        page = self.pages["pooled"]
        self.assertNotIn("Geographic precision", images(page))
        self.assertIn("no map or coordinate plot is drawn", section(page, "temporal"))
        self.assertIn("Named place or site, no coordinates:", section(page, "temporal"))
        # Coordinates in the individual metadata: the geography chart is shown.
        self.assertIn("Geographic precision", images(self.pages["combined"]))

    def test_missing_dates_leave_no_time_series_and_say_why(self):
        page = self.pages["pooled_no_dates"]
        for title in ["Mutation frequencies over time with uncertainty", "Variant frequencies over time",
                      "Sampling dates and sites"]:
            self.assertNotIn(title, images(page))
        self.assertIn("Temporal frequency trends", section(page, "overview"))
        self.assertNotIn("Temporal patterns available", page)

    def test_exploratory_population_structure_keeps_its_label(self):
        advanced = section(self.pages["pooled"], "advanced")
        self.assertIn("Exploratory — not the paper's method", advanced)
        self.assertIn("Exploratory alternative (not the paper", advanced)

    def test_experimental_grouping_is_collapsed_and_off_by_default(self):
        page = self.pages["annotated"]
        panel = re.search(r'<details class="experimental-panel" id="mutation-frequency-groups"([^>]*)>', page)
        self.assertIsNotNone(panel)
        self.assertNotIn("open", panel.group(1))
        self.assertIn("Profiles only: no groups shown", page)

    def test_mutation_table_is_collapsed_and_the_full_csv_is_offered(self):
        page = self.pages["annotated"]
        self.assertRegex(page, r"<summary>Detailed mutation table \(first \d+ of \d+ mutations; the CSV has every row\)")
        annotated = pd.read_csv(io.BytesIO(self.download(page, "sprint5_annotated_mutations.csv")))
        self.assertGreater(len(annotated), 0)

    def test_hypothetical_rho_scenarios_are_behind_a_labelled_disclosure(self):
        quality = section(self.pages["pooled"], "quality")
        self.assertIn("Not calibrated", quality)
        self.assertIn("<summary>Model comparison and hypothetical overdispersion scenarios</summary>", quality)

    # --- every figure and file stays available ---------------------------------------------------

    def test_every_figure_is_shown_or_explained_and_stays_downloadable(self):
        for name, page in self.pages.items():
            downloads = section(page, "downloads")
            figures = set(re.findall(r'/results/[^"]+/([a-z_0-9]+)\.png"', downloads))
            shown = set(re.findall(r'href="/results/[^"]+/([a-z_0-9]+)\.png">Download PNG', page))
            hidden = figures - shown
            for figure in hidden:
                self.assertRegex(downloads, rf"{figure}\.png</a> <span class=\"help inline-help\">\(chart not shown",
                                 f"{name}: {figure}")

    def test_downloads_are_grouped_once_and_include_the_reports(self):
        page = self.pages["annotated"]
        downloads = section(page, "downloads")
        names = re.findall(r'<li><a href="/results/[^"]+/([^"/]+)">', downloads)
        self.assertEqual(len(names), len(set(names)))  # each file listed once
        all_files = set(re.findall(r'href="/results/[^"]+/([^"/]+\.(?:csv|png))"', page))
        self.assertEqual(all_files, set(names))  # nothing that exists is missing from the grouped list
        self.assertIn("Results workbook (XLSX)", downloads)
        self.assertIn("Analysis report (PDF)", downloads)
        self.assertNotIn('class="section-downloads"', page)  # no scattered per-section lists

    # --- counts and their denominators -------------------------------------------------------------

    def test_displayed_counts_match_their_denominators(self):
        page = self.pages["annotated"]
        table = re.search(r'<table class="counts">(.*?)</table>', page, re.S).group(1)
        shown = {unescape(label): int(value.replace(",", "")) for label, value in
                 re.findall(r"<tr><td>([^<]+)</td><td>([\d,]+)</td>", table)}
        records = pd.read_csv(io.BytesIO(self.download(page, "pooled_variants_with_qc.csv")))
        observations = pd.read_csv(io.BytesIO(self.download(page, "sprint5_vaf_observations.csv")))
        positions = pd.read_csv(io.BytesIO(self.download(page, "variant_positions.csv")))
        annotated = pd.read_csv(io.BytesIO(self.download(page, "sprint5_annotated_mutations.csv")))
        included = records["Analysis_Included"].astype(bool)
        self.assertEqual(shown["SNV observations in the upload"], len(records))
        self.assertEqual(shown["Variant records after QC"], int(included.sum()))
        self.assertEqual(shown["Unique mutation observations after duplicate handling"],
                         int(observations["Counted_Once"].astype(bool).sum()))
        self.assertEqual(shown["Distinct mutations after QC"], records.loc[included, "Variant_Key"].nunique())
        self.assertEqual(shown["Distinct mutations after QC"], len(positions))
        self.assertEqual(shown["Distinct mutations evaluated for annotation"], records["Variant_Key"].nunique())
        self.assertEqual(shown["Distinct mutations evaluated for annotation"], len(annotated))
        # The synthetic table was built so that every stage differs.
        self.assertEqual(len(records), 18)
        self.assertEqual(len(set(shown.values())), len(shown))

    def test_stat_labels_name_their_processing_stage(self):
        text = text_of(self.pages["annotated"])
        self.assertIn("Variant records after QC: 17 (of 18 SNV observations in the upload)", text)
        self.assertIn("Distinct mutations evaluated for annotation:", text)
        self.assertIn("Unique mutation observations after duplicate handling: 16", text)

    # --- Sprint 4, print layout, and the generated report ---------------------------------------

    def test_sprint4_is_a_small_link_when_enrichment_was_not_run(self):
        page = self.pages["pooled"]
        self.assertNotIn("Metadata Enrichment (Sprint 4)", page)
        self.assertNotIn('id="enrichment-report"', page)
        self.assertNotIn("Claude", section(page, "overview"))

    def test_print_styles_hide_navigation_and_keep_figures_with_captions(self):
        page = self.pages["pooled"]
        print_css = page[page.index("@media print"):]
        self.assertIn("nav.sections", print_css.split("}")[2])
        self.assertIn("figure, .summary-grid p", print_css)
        self.assertIn("break-inside: avoid", print_css)

    def test_report_summary_uses_the_same_availability_and_main_figures(self):
        page = self.pages["pooled"]
        summary = json.loads(self.download_summary(page))
        availability = summary["availability"]
        self.assertEqual(availability["title"], "Pooled Pathogen Sequencing — Sprint 3 + 5")
        figures = report_figure_names(summary["figures"], summary)
        self.assertNotIn("quality", figures)
        self.assertNotIn("geography", figures)
        self.assertIn("depth_distribution", figures)
        for metric in availability["key_metrics"]:
            self.assertIn(metric, text_of(page))  # the report states the same headline counts as the page
        self.assertTrue(self.download(page, "analysis_report.pdf").startswith(b"%PDF"))

    def download_summary(self, page):
        download_id = re.search(r"/results/([^/]+)/", page).group(1)
        response = self.client.get(f"/results/{download_id}/analysis_summary.json")
        self.assertEqual(response.status_code, 200)
        return response.data


if __name__ == "__main__":
    unittest.main()
