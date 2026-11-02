"""Tests for published-annotation import, iVar import, and paper-based QC."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from genomics.pooled_observations import pooled_sample_records, pooled_variant_records
from genomics.quality_control import (
    EXPLORATORY_MODE,
    PAPER_MODE,
    apply_pooled_quality_control,
    load_masked_sites,
)
from wastewater.processed_snv import IVAR_FORMAT, PAPER_FORMAT, load_processed_snvs


class PublishedAnnotationTest(unittest.TestCase):
    def test_published_annotations_are_kept_with_a_published_prefix(self):
        source = pd.DataFrame({
            "SAMPLE": ["253"], "COLLECTION_DATE": ["30_May_2020"], "STATE": ["AZ"],
            "LOCATION": ["Guadalupe"], "POSITION": ["23403"], "REF": ["A"], "ALT": ["G"],
            "ALT_FREQ": ["0.9"], "TOTAL_DP": ["500"], "REF_AA": ["D"], "ALT_AA": ["G"],
            "Sy/Ns": ["Ns"], "Annotation": ["S"], "GLOBAL_FREQUENCY": ["0.7"],
            "Novel SNVs based on Global frequency": ["No"],
        })
        source.attrs["fontenele_2021_table2"] = True
        with patch("wastewater.processed_snv.read_table", return_value=source):
            data, _ = load_processed_snvs("mmc2.xlsx")
        self.assertEqual(data.attrs["input_format"], PAPER_FORMAT)
        self.assertEqual(data.loc[0, "Published_Ref_AA"], "D")
        self.assertEqual(data.loc[0, "Published_Alt_AA"], "G")
        self.assertEqual(data.loc[0, "Published_Syn_NonSyn"], "Ns")
        self.assertEqual(data.loc[0, "Published_Annotation"], "S")
        self.assertEqual(data.loc[0, "Published_Novel_Global"], "No")

        variants = pooled_variant_records(data, "mmc2.xlsx")
        self.assertEqual(variants.loc[0, "Published_Global_Frequency"], "0.7")

    def test_ivar_tsv_is_read_and_indels_are_set_aside(self):
        text = (
            "SAMPLE\tREGION\tPOS\tREF\tALT\tREF_DP\tALT_DP\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\n"
            "W1\tREF.1\t10\tA\tG\t70\t30\t0.3\t100\t0.001\tTRUE\n"
            "W1\tREF.1\t20\tC\t+TT\t90\t10\t0.1\t100\t0.01\tTRUE\n"
            "W2\tREF.1\t10\tA\tG\t5\t5\t0.5\t10\t0.2\tFALSE\n"
        )
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "ivar.tsv"
            path.write_text(text, encoding="utf-8")
            data, warnings = load_processed_snvs(path)
        self.assertEqual(data.attrs["input_format"], IVAR_FORMAT)
        self.assertEqual(len(data), 2)
        self.assertIn("insertion/deletion", " ".join(warnings))
        self.assertEqual(list(data["Alt_Depth"]), [30, 5])
        self.assertEqual(list(data["Source_Fisher_P"]), ["0.001", "0.2"])
        self.assertEqual(data.loc[0, "Reference_Genome"], "REF.1")

    def test_alternate_count_cannot_exceed_depth(self):
        text = (
            "Sample_ID,Reference_Genome,Position,Ref,Alt,Alt_Frequency,Read_Depth,Alt_Depth\n"
            "W1,REF.1,10,A,G,0.5,10,12\n"
        )
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "snvs.csv"
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "larger than the total read depth"):
                load_processed_snvs(path)


class PaperQualityControlTest(unittest.TestCase):
    def build_inputs(self):
        snvs = pd.DataFrame({
            "Sample_ID": ["W1", "W1", "W2", "W3"],
            "Reference_Genome": ["REF.1"] * 4,
            "Position": [10, 20, 10, 30],
            "Ref": ["A", "C", "A", "G"],
            "Alt": ["G", "T", "G", "A"],
            "Alt_Frequency": [0.3, 0.1, 0.5, 0.2],
            "Read_Depth": pd.array([100, 10, 50, pd.NA], dtype="Int64"),
            "Collection_Date": pd.to_datetime(["2026-01-01", "2026-01-01", "2026-01-08", "2026-01-15"]),
            "Location_Name": ["Site A", "Site A", "Site B", "Site B"],
            "Source_Fisher_P": ["0.01", "0.01", "0.2", ""],
        })
        metadata = pd.DataFrame({
            "Sample_ID": ["W1", "W2", "W3"],
            "Genome_Coverage_Percent": [95.0, 50.0, pd.NA],
            "Mean_Depth": [200.0, 30.0, pd.NA],
            "Latitude": [pd.NA] * 3, "Longitude": [pd.NA] * 3,
        })
        variants = pooled_variant_records(snvs, "test")
        samples = pooled_sample_records(snvs, "test", pooled_metadata=metadata)
        return variants, samples

    def test_paper_mode_excludes_only_evaluable_failures(self):
        variants, samples = self.build_inputs()
        checked, checked_samples, report = apply_pooled_quality_control(
            variants, samples, PAPER_MODE, masked_sites={("REF.1", 20)}
        )
        rows = checked.set_index(["Sample_ID", "Position"])
        self.assertEqual(rows.loc[("W1", 10), "Quality_Status"], "Retained")
        self.assertIn("variant depth", rows.loc[("W1", 20), "Quality_Status"])
        self.assertIn("masked site", rows.loc[("W1", 20), "Quality_Status"])
        self.assertIn("fisher test", rows.loc[("W2", 10), "Quality_Status"])
        self.assertIn("sample coverage", rows.loc[("W2", 10), "Quality_Status"])
        self.assertFalse(rows.loc[("W2", 10), "Analysis_Included"])
        # W3 has no depth, p-value, or coverage: nothing failed, nothing was proven.
        self.assertEqual(rows.loc[("W3", 30), "Quality_Status"], "Retained")
        self.assertTrue(rows.loc[("W3", 30), "QC_Variant_Depth"].startswith("Not evaluable"))
        self.assertTrue(rows.loc[("W3", 30), "QC_Fisher_Test"].startswith("Not evaluable"))
        self.assertEqual(rows.loc[("W3", 30), "QC_Base_Quality"],
                         "Not verifiable from processed records")

        samples_by_id = checked_samples.set_index("Sample_ID")
        self.assertTrue(samples_by_id.loc["W1", "Analysis_Included"])
        self.assertFalse(samples_by_id.loc["W2", "Analysis_Included"])
        self.assertTrue(samples_by_id.loc["W3", "QC_Genome_Coverage"].startswith("Not evaluable"))

        steps = report.set_index("Step")["Status"]
        self.assertEqual(steps["Minimum base quality"], "Not verifiable")
        self.assertEqual(steps["Minimum variant depth"], "Partly applied (some values missing)")
        self.assertEqual(steps["Problematic-site masking"], "Applied to supplied list")
        self.assertEqual(steps["Minimum allele frequency"], "Matched")
        self.assertEqual(steps["Read mapping to the reference"], "Not reproduced here")

    def test_exploratory_mode_flags_without_excluding(self):
        variants, samples = self.build_inputs()
        checked, checked_samples, report = apply_pooled_quality_control(
            variants, samples, EXPLORATORY_MODE
        )
        self.assertTrue(checked["Analysis_Included"].all())
        self.assertTrue(checked_samples["Analysis_Included"].all())
        self.assertIn("Flagged", " ".join(checked["Quality_Status"]))
        self.assertTrue(checked["QC_Masked_Site"].str.startswith("Not reproduced").all())
        self.assertEqual(report.set_index("Step").loc["Problematic-site masking", "Status"],
                         "Not reproduced")

    def test_masked_site_list_requires_reference_and_position(self):
        with tempfile.TemporaryDirectory() as folder:
            good = Path(folder) / "mask.csv"
            good.write_text("Reference_ID,Position\nREF.1,20\n", encoding="utf-8")
            bad = Path(folder) / "bad.csv"
            bad.write_text("Position\n20\n", encoding="utf-8")
            self.assertEqual(load_masked_sites(good), {("REF.1", 20)})
            with self.assertRaisesRegex(ValueError, "Reference_ID and Position"):
                load_masked_sites(bad)


if __name__ == "__main__":
    unittest.main()
