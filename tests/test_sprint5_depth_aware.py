"""Sprint 5 Phase 1: mutation-frequency time series, website section, and upload-to-results tests.

Time-series expectations are worked out by hand from the synthetic example
files (examples/integrated). The web tests go through the real Flask routes
(/review then /run, then the result downloads). MAFFT is mocked with
tests/fake_mafft.py; GeoNames and Claude are off.
"""

import io
import math
import os
import re
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402

from app import STAGED, app  # noqa: E402
from genomics.frequency_uncertainty import Z_95  # noqa: E402
from genomics.integrated_analysis import run_integrated_analysis  # noqa: E402
from genomics.vaf_trends import NON_DETECTION, choose_series  # noqa: E402
from integrated_view import build_dashboard  # noqa: E402

EXAMPLES = Path(__file__).resolve().parents[1] / "examples" / "integrated"
REFERENCE = (EXAMPLES / "synthetic_reference.fasta").read_text(encoding="utf-8").split("\n", 1)[1].replace("\n", "")
SNVS = pd.read_csv(EXAMPLES / "synthetic_pooled_snvs.csv", dtype=str, keep_default_na=False)


def run_pooled(snvs=None, counts=True, reference=True):
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / "snvs.csv"
        (SNVS if snvs is None else snvs).to_csv(path, index=False)
        return run_integrated_analysis(
            snv_path=path,
            position_path=EXAMPLES / "synthetic_pooled_base_counts.csv" if counts else None,
            reference_path=EXAMPLES / "synthetic_reference.fasta" if reference else None)


class TimeSeriesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = run_pooled()
        cls.summary = cls.result["tables"]["sprint5_temporal_summary.csv"].set_index(
            ["Variant_Key", "Location_Name"])

    def test_supported_non_detection_enters_the_series(self):
        # South site, 15 T>C: WW_S1 (Mar 3) has 0 C reads in the base counts, WW_S2 (Mar 10)
        # 8/160 = 0.05, WW_S3 (Mar 17) 20/200 = 0.10; WW_S4 failed coverage QC.
        row = self.summary.loc[("SYNREF.1:15:T>C", "South Synthetic Sewershed")]
        self.assertEqual(row["Supported_Non_Detections"], 1)
        self.assertEqual(row["Detections"], 2)
        self.assertEqual(row["Earliest_VAF"], 0.0)
        self.assertAlmostEqual(row["Latest_VAF"], 0.10)
        self.assertAlmostEqual(row["Absolute_Change"], 0.10)
        # Days 0, 7, 14 and VAF 0, 0.05, 0.10 lie on a line of 0.05 per week.
        self.assertAlmostEqual(row["Slope_Per_Week"], 0.05)
        # Wilson upper bound for 0 of 120 reads: z^2 / (120 + z^2).
        self.assertAlmostEqual(row["Earliest_Upper"], Z_95 ** 2 / (120 + Z_95 ** 2))
        self.assertEqual(row["Earliest_Lower"], 0.0)

    def test_qc_excluded_call_is_unknown_not_zero(self):
        # North, 30 G>T: WW_N1's call fails the Fisher filter, so it is unknown; WW_N3 has 0 reads.
        row = self.summary.loc[("SYNREF.1:30:G>T", "North Synthetic Sewershed")]
        self.assertEqual(row["Unknown_Or_Uncallable"], 1)
        self.assertEqual(row["Supported_Non_Detections"], 1)
        self.assertAlmostEqual(row["Earliest_VAF"], 14 / 180)  # count-derived, reported 0.08
        self.assertIn("left out", row["Interpretation_Note"])

    def test_references_are_never_merged(self):
        self.assertIn(("SYNREF_B.1:15:T>C", "North Synthetic Sewershed"), self.summary.index)
        row = self.summary.loc[("SYNREF_B.1:15:T>C", "North Synthetic Sewershed")]
        self.assertEqual(row["Known_Dated_Time_Points"], 1)
        self.assertEqual(row["Series_Status"], "Single time point")
        self.assertTrue(math.isnan(row["Slope_Per_Week"]))

    def test_without_base_counts_unreported_samples_stay_unknown(self):
        result = run_pooled(counts=False)
        summary = result["tables"]["sprint5_temporal_summary.csv"]
        self.assertEqual(int(summary["Supported_Non_Detections"].sum()), 0)
        south = summary.set_index(["Variant_Key", "Location_Name"]).loc[
            ("SYNREF.1:15:T>C", "South Synthetic Sewershed")]
        self.assertEqual(south["Unknown_Or_Uncallable"], 2)  # WW_S1 unknown, WW_S4 excluded sample
        self.assertAlmostEqual(south["Earliest_VAF"], 0.05)  # first known value, not a fabricated 0
        self.assertTrue(any("unknown rather than a confirmed non-detection" in text
                            for text in result["depth_aware"]["messages"]))

    def test_two_dates_give_a_change_but_no_slope(self):
        snvs = SNVS[SNVS["Sample_ID"].isin(["WW_N1", "WW_N2"])]
        summary = run_pooled(snvs, counts=False, reference=False)["tables"]["sprint5_temporal_summary.csv"]
        row = summary.set_index(["Variant_Key", "Location_Name"]).loc[
            ("SYNREF.1:15:T>C", "North Synthetic Sewershed")]
        self.assertAlmostEqual(row["Absolute_Change"], 0.25 - 0.10)
        self.assertTrue(math.isnan(row["Slope_Per_Week"]))
        self.assertEqual(row["Series_Status"], "Change between two dates")

    def test_selection_prefers_series_with_more_dated_values(self):
        chosen = choose_series(self.result["tables"]["sprint5_temporal_summary.csv"], limit=2)
        self.assertEqual(len(chosen), 2)
        self.assertTrue(all(location == "North Synthetic Sewershed" for _, location in chosen))
        self.assertIn("vaf_time_series", self.result["figures"])

    def test_non_detection_points_have_intervals(self):
        self.assertEqual(NON_DETECTION, "Supported non-detection (0 alternative reads)")
        self.assertIn("depth_dependence", self.result["figures"])
        self.assertIn("vaf_depth", self.result["figures"])


class AvailabilityTest(unittest.TestCase):
    def test_no_dates(self):
        snvs = SNVS.assign(Collection_Date="")
        result = run_pooled(snvs)
        self.assertIn("Collection dates are unavailable, so temporal trends cannot be calculated.",
                      result["depth_aware"]["messages"])
        self.assertNotIn("vaf_time_series", result["figures"])
        self.assertEqual(result["notes"].get("figures", []), [])
        self.assertEqual(result["depth_aware"]["eligible_series"], 0)

    def test_no_depth(self):
        snvs = SNVS.assign(Read_Depth="", Alt_Depth="")
        result = run_pooled(snvs, counts=False)
        self.assertIn("Read depth is unavailable, so frequency uncertainty cannot be calculated.",
                      result["depth_aware"]["messages"])
        self.assertEqual(result["depth_aware"]["with_interval"], 0)
        self.assertNotIn("vaf_depth", result["figures"])

    def test_frequency_only_without_reproducible_counts(self):
        # 0.3 at 7 reads is not a whole number of reads, so nothing is reconstructed.
        snvs = SNVS.drop(columns=["Alt_Depth"]).copy()
        snvs.loc[0, ["Alt_Frequency", "Read_Depth"]] = ["0.3", "7"]
        result = run_pooled(snvs, counts=False)
        self.assertTrue(any(message.startswith("The data contain reported frequencies but no verified read counts")
                            for message in result["depth_aware"]["messages"]))
        self.assertEqual(result["depth_aware"]["with_interval"], 0)

    def test_dashboard_section(self):
        result = run_pooled()
        view = build_dashboard(result, "test")
        section = view["depth_aware"]
        self.assertEqual(section["title"], "Sprint 5 — Depth-Aware Variant Analysis")
        labels = [item["label"] for item in section["stats"]]
        for label in ["Samples evaluated", "Unique mutation observations after duplicate handling",
                      "Observations with usable depth",
                      "Observations with verified or reconstructed alternative read counts",
                      "Observations with calculated uncertainty", "Observations missing required information",
                      "Variants eligible for time-series analysis"]:
            self.assertIn(label, labels)
        self.assertIn("sprint5_vaf_observations.csv", section["downloads"])
        text = " ".join(section["limitations"]).lower()
        self.assertNotIn("tumor", text)
        self.assertNotIn("clone", text)


def token_of(page):
    match = re.search(r'name="token" value="([^"]+)"', page)
    return match.group(1) if match else None


class WebWorkflowTest(unittest.TestCase):
    """Real routes: /review, /run, and the result downloads."""

    def setUp(self):
        patches = [patch.dict(os.environ, {"GEONAMES_USERNAME": "", "AI_ENRICHMENT_ENABLED": "false"}),
                   patch("align_sequences.shutil.which", return_value="mafft"),
                   patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(REFERENCE))]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        self.addCleanup(STAGED.clear)
        self.client = app.test_client()

    def upload(self, files, mode="pooled", **form):
        uploads = [(io.BytesIO(content), name) for name, content in files]
        review = self.client.post("/review", data={"files": uploads, "intended_mode": mode, **form},
                                  content_type="multipart/form-data")
        self.assertEqual(review.status_code, 200)
        page = review.get_data(as_text=True)
        token = token_of(page)
        if token is None:
            return page, None, None
        with redirect_stdout(io.StringIO()):
            response = self.client.post("/run", data={"token": token, "intended_mode": mode, **form})
        self.assertEqual(response.status_code, 200)
        result = response.get_data(as_text=True)
        download = re.search(r"/results/([^/]+)/", result)
        return page, result, download.group(1) if download else None

    def download(self, download_id, name):
        response = self.client.get(f"/results/{download_id}/{name}")
        self.assertEqual(response.status_code, 200, name)
        return pd.read_csv(io.BytesIO(response.data)) if name.endswith(".csv") else response.data

    @staticmethod
    def example(name):
        return name, (EXAMPLES / name).read_bytes()

    def test_sprint3_alone_with_counts_and_reference(self):
        _, page, download_id = self.upload([self.example("synthetic_pooled_snvs.csv"),
                                            self.example("synthetic_pooled_base_counts.csv"),
                                            self.example("synthetic_reference.fasta")])
        self.assertIn("Sequencing Quality and Uncertainty", page)
        self.assertIn("Observations with calculated uncertainty", page)
        self.assertNotIn("could not be drawn", page)
        evidence = self.download(download_id, "sprint5_vaf_observations.csv")
        self.assertEqual(len(evidence), 22)
        self.assertTrue((evidence["Alt_Count_Source"] == "Supplied").all())
        self.assertIn("Wilson_Lower", evidence.columns)
        summary = self.download(download_id, "sprint5_temporal_summary.csv")
        self.assertGreater(int(summary["Supported_Non_Detections"].sum()), 0)
        missing = self.download(download_id, "sprint5_missing_information.csv")
        self.assertIn("Requirement", missing.columns)
        self.download(download_id, "vaf_time_series.png")
        self.download(download_id, "reference_verification.csv")
        # Sprint 3 outputs are still produced.
        for name in ["pooled_variants_with_qc.csv", "temporal_trends.csv", "pcoa_coordinates.csv",
                     "yue_clayton_pairs.csv"]:
            self.download(download_id, name)

    def test_sprint2_and_3_combined(self):
        _, page, download_id = self.upload(
            [self.example("synthetic_individual_sequences.fasta"), self.example("synthetic_individual_metadata.csv"),
             self.example("synthetic_reference.fasta"), self.example("synthetic_pooled_snvs.csv"),
             self.example("synthetic_pooled_base_counts.csv")], mode="combined")
        self.assertIn("Sequencing Uncertainty and Variant Frequency", page)
        self.assertIn("Individual Sequence Relationships", page)
        self.download(download_id, "cross_source_comparison.csv")
        self.download(download_id, "pairwise_distances.csv")
        self.download(download_id, "sprint5_vaf_observations.csv")

    def test_direct_counts_without_dates_or_coverage_table(self):
        snvs = SNVS.assign(Collection_Date="").to_csv(index=False).encode()
        _, page, download_id = self.upload([("no_dates.csv", snvs)])
        self.assertIn("Collection dates are unavailable, so temporal trends cannot be calculated.", page)
        self.assertIn("Temporal trends unavailable because collection dates were not provided", page)
        self.assertNotIn("could not be drawn", page)
        evidence = self.download(download_id, "sprint5_vaf_observations.csv")
        self.assertTrue(evidence["Wilson_Lower"].notna().all())

    def test_sparse_frequency_only_table(self):
        snvs = SNVS.drop(columns=["Alt_Depth"]).copy()
        snvs.loc[0, ["Alt_Frequency", "Read_Depth"]] = ["0.3", "7"]
        _, page, download_id = self.upload([("sparse.csv", snvs.to_csv(index=False).encode())])
        self.assertIn("no verified read counts", page)
        evidence = self.download(download_id, "sprint5_vaf_observations.csv")
        self.assertTrue(evidence["Wilson_Lower"].isna().all())
        self.assertTrue(evidence["Alt_Read_Count"].isna().all())

    def test_coverage_gaps_stay_unknown(self):
        counts = pd.read_csv(EXAMPLES / "synthetic_pooled_base_counts.csv")
        counts = counts[~((counts["Sample_ID"] == "WW_S1") & (counts["Position"] == 15))]
        _, page, download_id = self.upload([self.example("synthetic_pooled_snvs.csv"),
                                            ("counts_gap.csv", counts.to_csv(index=False).encode())])
        summary = self.download(download_id, "sprint5_temporal_summary.csv").set_index(
            ["Variant_Key", "Location_Name"])
        row = summary.loc[("SYNREF.1:15:T>C", "South Synthetic Sewershed")]
        self.assertEqual(row["Supported_Non_Detections"], 0)  # the gap is not a 0
        self.assertGreaterEqual(row["Unknown_Or_Uncallable"], 1)

    def test_ambiguous_dates_are_flagged_then_resolved_by_choice(self):
        snvs = SNVS.copy()
        snvs["Collection_Date"] = pd.to_datetime(snvs["Collection_Date"]).dt.strftime("%d/%m/%Y")
        content = snvs.to_csv(index=False).encode()
        review, page, download_id = self.upload([("dmy.csv", content)])
        # Days 01-12 (for example 01/03/2026) are ambiguous; 15/03/2026 can only be day first.
        self.assertIn("day and month order is ambiguous", review)
        self.assertIn("can only be read day first", review)
        evidence = self.download(download_id, "sprint5_vaf_observations.csv")
        ambiguous = evidence[evidence["Collection_Date_Status"].str.startswith("Ambiguous")]
        self.assertGreater(len(ambiguous), 0)
        self.assertTrue(ambiguous["Collection_Date"].isna().all())
        self.assertIn("01/03/2026", set(ambiguous["Collection_Date_Original"]))
        _, page, download_id = self.upload([("dmy.csv", content)], pooled_date_order="day_first")
        evidence = self.download(download_id, "sprint5_vaf_observations.csv")
        self.assertEqual(set(evidence["Collection_Date_Status"]), {"Parsed with the selected day/month order"})
        first = evidence[evidence["Collection_Date_Original"].eq("01/03/2026")]
        self.assertEqual(set(first["Collection_Date"]), {"2026-03-01"})

    def test_invalid_and_unsupported_inputs_give_controlled_messages(self):
        bad = SNVS.copy()
        bad.loc[0, "Read_Depth"] = "-5"
        _, page, _ = self.upload([("bad_depth.csv", bad.to_csv(index=False).encode())])
        self.assertIn("read depth", page.lower())
        self.assertNotIn("Traceback", page)
        review, page, _ = self.upload([("notes.docx", b"PK\x03\x04 not really a document")])
        self.assertNotIn("Traceback", review)
        self.assertIn("Resolve the problems listed in Step 2", page)  # /run refuses; nothing is analyzed
        self.assertNotIn('id="depth-aware"', page)

    def test_other_routes_still_work(self):
        self.assertEqual(self.client.get("/health").status_code, 200)
        self.assertEqual(self.client.get("/").status_code, 200)
        self.assertEqual(self.client.get("/examples/synthetic_pooled_snvs.csv").status_code, 200)

    @unittest.skipUnless(os.environ.get("FONTENELE_MMC2"), "set FONTENELE_MMC2 to the published mmc2.xlsx")
    def test_published_workbook(self):
        path = Path(os.environ["FONTENELE_MMC2"])
        _, page, download_id = self.upload([(path.name, path.read_bytes())])
        self.assertIn("Sprint 5 — Depth-Aware Variant Analysis", page)
        self.assertIn("reconstructed for 7973 observation(s)", page)
        self.assertNotIn("could not be drawn", page)
        evidence = self.download(download_id, "sprint5_vaf_observations.csv")
        self.assertEqual(len(evidence), 7973)
        self.assertEqual(int((evidence["Interval_Status"] == "Approximate (reconstructed counts)").sum()), 7970)
        summary = self.download(download_id, "sprint5_temporal_summary.csv")
        # Independent count from the raw workbook: mutation-site pairs reported on 2+ dates
        # (every published record passes QC except the 12 below 20x depth).
        raw = pd.read_excel(path, header=1, dtype=str, keep_default_na=False)
        raw = raw[pd.to_numeric(raw["TOTAL_DP"]) >= 20]
        raw["Site"] = raw["STATE"].str.strip() + ": " + raw["LOCATION"].str.strip()
        dates = raw.groupby(["Site", "POSITION", "REF", "ALT"])["COLLECTION_DATE"].nunique()
        self.assertEqual(int((summary["Known_Dated_Time_Points"] >= 2).sum()), int((dates >= 2).sum()))


if __name__ == "__main__":
    unittest.main()
