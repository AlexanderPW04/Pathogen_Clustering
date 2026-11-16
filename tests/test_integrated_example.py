"""End-to-end test of the SYNTHETIC paired individual + pooled demonstration.

The example in examples/integrated/ was designed so that every result below
is known in advance (see make_integrated_example.EXPECTED_CATEGORIES). MAFFT
is replaced by tests/fake_mafft.py, so this is a mocked-MAFFT test.
"""

import contextlib
import io
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402

from genomics.diversity_analysis import PAPER_PCOA  # noqa: E402
from genomics.integrated_analysis import run_integrated_analysis  # noqa: E402
from genomics.quality_control import EXPLORATORY_MODE  # noqa: E402
from genomics.variant_matching import CATEGORY_ORDER  # noqa: E402
from run_pipeline import run_pipeline  # noqa: E402

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "integrated"
sys.path.insert(0, str(EXAMPLE))
import make_integrated_example as design  # noqa: E402


def run_example(qc_mode="paper"):
    """Run the Sprint 2 pipeline, then the integrated analysis, on the example."""

    fake = make_fake_mafft(design.REFERENCE)
    with tempfile.TemporaryDirectory() as folder:
        fasta = Path(folder) / "synthetic_individual_sequences.fasta"
        metadata = Path(folder) / "synthetic_individual_metadata.csv"
        shutil.copyfile(EXAMPLE / fasta.name, fasta)
        shutil.copyfile(EXAMPLE / metadata.name, metadata)
        shared = {}
        with patch.dict(os.environ, {"GEONAMES_USERNAME": ""}):
            with patch("align_sequences.shutil.which", return_value="mafft"):
                with patch("align_sequences.subprocess.run", side_effect=fake):
                    with contextlib.redirect_stdout(io.StringIO()), \
                            contextlib.redirect_stderr(io.StringIO()):
                        status = run_pipeline(fasta, metadata, 0.05, 14, 50, analysis_data=shared)
                        result = run_integrated_analysis(
                            individual_data=shared.get("samples"),
                            mafft_path="mafft",
                            reference_path=EXAMPLE / "synthetic_reference.fasta",
                            snv_path=EXAMPLE / "synthetic_pooled_snvs.csv",
                            position_path=EXAMPLE / "synthetic_pooled_base_counts.csv",
                            signature_path=EXAMPLE / "synthetic_variant_signatures.csv",
                            settings={"qc_mode": qc_mode},
                            input_labels={"fasta": fasta.name, "snv": "synthetic_pooled_snvs.csv"},
                        )
        related = pd.read_csv(Path(folder) / (
            "synthetic_individual_metadata_pairwise_comparisons_temporal_geographic_related_pairs.csv"))
    return status, result, related


class IntegratedExampleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.status, cls.result, cls.related = run_example()

    def test_sprint_2_relationship_analysis_still_runs(self):
        self.assertEqual(self.status, 0)
        pairs = {tuple(sorted(pair)) for pair in zip(self.related["Sample_ID_1"], self.related["Sample_ID_2"])}
        # A-B are 7 days and about 1 km apart; B-C are 11 days and about 13 km apart.
        self.assertIn(("IND_A", "IND_B"), pairs)
        self.assertIn(("IND_B", "IND_C"), pairs)

    def test_individual_reference_profiles(self):
        individual = self.result["individual"]
        samples = individual["samples"].set_index("Sample_ID")
        self.assertEqual(samples.loc["IND_C", "Sequenced_Span"], "36-120")
        self.assertIn("not established", samples.loc["IND_E", "Variant_Profile_Status"])
        carried = individual["variants"].groupby("Sample_ID")["Variant_Key"].apply(set).to_dict()
        self.assertEqual(carried["IND_A"], {"SYNREF.1:15:T>C", "SYNREF.1:60:T>G"})
        self.assertEqual(carried["IND_B"], {"SYNREF.1:15:T>C", "SYNREF.1:100:T>C"})
        self.assertEqual(carried["IND_C"], {"SYNREF.1:42:T>A"})
        self.assertNotIn("IND_D", carried)  # its indels are not substitutions
        events = set(individual["indels"].loc[individual["indels"]["Sample_ID"].eq("IND_D"), "Event"])
        self.assertEqual(events, {"Insertion", "Deletion"})

    def test_every_variant_lands_in_its_designed_category(self):
        table = self.result["comparison"]["table"].set_index("Variant_Key")
        for key, expected in design.EXPECTED_CATEGORIES.items():
            with self.subTest(variant=key):
                self.assertTrue(table.loc[key, "Category"].startswith(expected),
                                f"{key}: {table.loc[key, 'Category']}")
        self.assertIn("different", table.loc["SYNREF_B.1:15:T>C", "Note"].lower())
        self.assertIn("mismatch", table.loc["SYNREF.1:110:G>A", "Note"])
        self.assertEqual(len(table), len(design.EXPECTED_CATEGORIES))

    def test_missing_coverage_stays_unresolved_not_negative(self):
        row = self.result["comparison"]["table"].set_index("Variant_Key").loc["SYNREF.1:100:T>C"]
        self.assertEqual(row["Pooled_Reported"], 0)
        self.assertEqual(row["Pooled_Adequate_Depth_Not_Reported"], 0)
        self.assertGreater(row["Pooled_Low_Or_No_Coverage"], 0)
        matrix = self.result["tables"]["pooled_observation_matrix.csv"]
        unreported = matrix[~matrix["Status"].eq("Variant observed")]
        self.assertTrue(unreported["Alt_Frequency"].isna().all())

    def test_shared_variant_dates_and_window(self):
        row = self.result["comparison"]["table"].set_index("Variant_Key").loc["SYNREF.1:15:T>C"]
        self.assertEqual(row["First_Individual_Date"], pd.Timestamp("2026-03-02"))
        self.assertEqual(row["First_Pooled_Date"], pd.Timestamp("2026-03-01"))
        self.assertEqual(row["Days_Pooled_Minus_Individual"], -1)
        self.assertGreater(row["Date_Pairs_Within_Window"], 0)
        self.assertTrue(pd.isna(row["Nearest_Km_Between_Carriers"]))  # sewersheds have no coordinates

    def test_paper_qc_exclusions(self):
        samples = self.result["pooled"]["samples"].set_index("Sample_ID")
        self.assertFalse(samples.loc["WW_S4", "Analysis_Included"])  # 78% coverage
        self.assertTrue(samples.loc["WW_N1", "Analysis_Included"])
        self.assertAlmostEqual(samples.loc["WW_S4", "Genome_Coverage_Percent"], 94 / 120 * 100)
        variants = self.result["pooled"]["variants"]
        n1_30 = variants[(variants["Sample_ID"] == "WW_N1") & (variants["Position"] == 30)].iloc[0]
        self.assertIn("fisher", n1_30["Quality_Status"])
        mismatch = variants[variants["Position"] == 110].iloc[0]
        self.assertIn("reference base", mismatch["Quality_Status"])
        # The excluded call is shown as excluded, not as a non-detection.
        matrix = self.result["tables"]["pooled_observation_matrix.csv"]
        cell = matrix[(matrix["Sample_ID"] == "WW_N1") & (matrix["Variant_Key"] == "SYNREF.1:30:G>T")]
        self.assertEqual(cell["Status"].iloc[0], "Reported but excluded by QC")

    def test_temporal_trend_is_numerically_correct(self):
        trends = self.result["tables"]["temporal_trends.csv"].set_index(["Variant_Key", "Location_Name"])
        north = trends.loc[("SYNREF.1:15:T>C", "North Synthetic Sewershed")]
        self.assertEqual(north["Reported"], 4)
        self.assertAlmostEqual(north["Change_In_Callable_Frequency"], 0.60)
        self.assertAlmostEqual(north["Slope_Per_Week"], 0.20)

    def test_paper_pcoa_is_available_and_other_reference_left_out(self):
        population = self.result["population"]
        self.assertTrue(population["available"])
        self.assertEqual(population["method"], PAPER_PCOA)
        self.assertIn("WW_X1", population["left_out"])
        self.assertNotIn("WW_S4", list(population["coordinates"]["Sample_ID"]))
        self.assertEqual(len(population["coordinates"]), 7)

    def test_signature_overlap_is_not_a_lineage_call(self):
        overlap = self.result["tables"]["variant_signature_overlap.csv"]
        row = overlap[(overlap["Signature_Name"] == "Synthetic profile 1")
                      & (overlap["Sample_ID"] == "IND_B")].iloc[0]
        self.assertEqual(row["Variants_Observed"], 2)  # positions 15 and 100
        self.assertTrue((overlap["Interpretation"] == "Profile overlap only; not a lineage assignment").all())

    def test_figures_provenance_and_summary(self):
        for name in ["variants_per_sample", "depth_vs_snvs", "genome_positions", "trajectories",
                     "heatmap", "pcoa", "cross_source", "cross_source_timeline", "quality"]:
            with self.subTest(figure=name):
                self.assertTrue(self.result["figures"][name]["png"].startswith(b"\x89PNG"))
        provenance = self.result["tables"]["provenance.csv"]
        self.assertIn("SYNREF.1", set(provenance["Value"].astype(str)))
        summary = self.result["summary"]
        self.assertEqual(summary["shared_variants"], 2)
        self.assertEqual(list(self.result["comparison"]["counts"]["Category"]), CATEGORY_ORDER)

    def test_exploratory_mode_keeps_flagged_records(self):
        _, result, _ = run_example(EXPLORATORY_MODE)
        self.assertTrue(result["pooled"]["samples"]["Analysis_Included"].all())
        self.assertEqual(len(result["population"]["coordinates"]), 8)


class CommandLineTest(unittest.TestCase):
    def test_command_line_example_writes_tables_and_figures(self):
        from run_integrated_analysis import main

        fake = make_fake_mafft(design.REFERENCE)
        with tempfile.TemporaryDirectory() as folder:
            with patch.dict(os.environ, {"GEONAMES_USERNAME": ""}):
                with patch("align_sequences.shutil.which", return_value="mafft"):
                    with patch("align_sequences.subprocess.run", side_effect=fake):
                        with contextlib.redirect_stdout(io.StringIO()), \
                                contextlib.redirect_stderr(io.StringIO()):
                            status = main(["--example", "--output-folder", folder])
            written = {path.name for path in Path(folder).iterdir()}
            comparison = pd.read_csv(Path(folder) / "cross_source_comparison.csv")
        self.assertEqual(status, 0)
        for name in ["cross_source_comparison.csv", "provenance.csv", "pcoa.png",
                     "pooled_variants_with_qc.csv",
                     "synthetic_individual_metadata_region_assignments.csv"]:
            self.assertIn(name, written)
        self.assertEqual(len(comparison), len(design.EXPECTED_CATEGORIES))

    def test_command_line_rejects_missing_inputs(self):
        from run_integrated_analysis import main

        with tempfile.TemporaryDirectory() as folder:
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(main(["--output-folder", folder]), 2)
                self.assertEqual(main(["--fasta", "missing.fasta", "--output-folder", folder]), 2)


if __name__ == "__main__":
    unittest.main()
