"""Regression tests for the Sprint 3 issues repaired at the start of Sprint 5.

Expected values are worked out by hand in the comments, not copied from the
code. Audit IDs refer to docs/SPRINT3_AUDIT_REPORT.md.
"""

import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from genomics.collection_dates import (
    AMBIGUOUS,
    INVALID,
    NOT_SUPPLIED,
    PARSED,
    PARSED_WITH_ORDER,
    PARTIAL,
    TWO_DIGIT_YEAR,
    parse_collection_dates,
    parse_one,
)
from genomics.diversity_analysis import (
    jaccard_distance_matrix,
    paper_distance_matrix,
    population_structure,
)
from genomics.integrated_analysis import position_summary, run_integrated_analysis
from genomics.pooled_observations import (
    REPEAT_CONFLICTING,
    REPEAT_IDENTICAL,
    UNIQUE_RECORD,
    add_frequency_checks,
    build_depth_lookup,
    mark_repeated_records,
)
from genomics.quality_control import (
    REFERENCE_MISMATCH,
    REFERENCE_PARTLY_MISMATCHED,
    REFERENCE_UNVERIFIED,
    REFERENCE_VERIFIED,
    reference_verification,
)
from genomics.temporal_analysis import (
    VALUE_UNKNOWN,
    ZERO_ALT_READS,
    frequency_trends,
    location_date_summary,
    observation_matrix,
)
from genomics.unified_schema import COVERED_NOT_DETECTED, VARIANT_OBSERVED
from wastewater.processed_snv import load_processed_snvs


EXAMPLES = Path(__file__).resolve().parents[1] / "examples" / "integrated"
DATES = pd.to_datetime(["2026-01-01", "2026-01-08", "2026-01-15", "2026-01-22"])


def four_samples():
    return pd.DataFrame({
        "Sample_ID": ["A", "B", "C", "E"], "Reference_ID": ["R.1"] * 4, "Collection_Date": DATES,
        "Location_Name": ["Site"] * 4, "Analysis_Included": [True] * 4,
    })


def reported_variant(sample_ids, frequencies):
    count = len(sample_ids)
    return pd.DataFrame({
        "Variant_Key": ["R.1:5:A>G"] * count, "Reference_ID": ["R.1"] * count, "Position": [5] * count,
        "Ref": ["A"] * count, "Alt": ["G"] * count, "Sample_ID": sample_ids,
        "Alt_Frequency": frequencies, "Read_Depth": [100] * count, "Analysis_Included": [True] * count,
    })


def position_table(rows, with_counts=True):
    """rows: (sample, depth, g_count); position 5 only."""

    data = pd.DataFrame([{"Sample_ID": sample, "Reference_ID": "R.1", "Position": 5, "Depth": depth,
                          "A_Count": depth - g, "C_Count": 0, "G_Count": g, "T_Count": 0}
                         for sample, depth, g in rows])
    if not with_counts:
        data = data[["Sample_ID", "Reference_ID", "Position", "Depth"]]
    data.attrs["has_base_counts"] = with_counts
    return data


class CallableTrendTest(unittest.TestCase):
    """1A / S3-02: non-detections with explicit evidence enter the trend; unknowns never become 0."""

    def matrix(self, with_counts, b_alt_reads=0):
        variants = reported_variant(["A", "C", "E"], [0.2, 0.4, 0.6])
        table = position_table([("A", 100, 20), ("B", 100, b_alt_reads), ("C", 100, 40), ("E", 100, 60)],
                               with_counts)
        return observation_matrix(variants, four_samples(), build_depth_lookup(table), 20, table)

    def test_zero_alternative_reads_enter_the_trend_as_zero(self):
        matrix = self.matrix(with_counts=True)
        b = matrix.set_index("Sample_ID").loc["B"]
        self.assertEqual(b["Status"], COVERED_NOT_DETECTED)
        self.assertEqual(b["Callable_Evidence"], ZERO_ALT_READS)
        self.assertEqual(b["Callable_Frequency"], 0.0)
        self.assertTrue(math.isnan(b["Alt_Frequency"]))  # the reported column is never filled in
        trend = frequency_trends(matrix).iloc[0]
        # x = 0, 7, 14, 21 days; y = 0.2, 0, 0.4, 0.6: mean x 10.5, mean y 0.3,
        # Sxy = 1.05 + 1.05 + 0.35 + 3.15 = 5.6, Sxx = 245 -> 0.16 per week.
        self.assertAlmostEqual(trend["Slope_Per_Week"], 0.16)
        self.assertAlmostEqual(trend["Change_In_Callable_Frequency"], 0.4)
        self.assertEqual(trend["Zero_Alt_Reads"], 1)
        self.assertEqual(trend["Callable_Dated"], 4)
        self.assertTrue(trend["Complete_Series"])

    def test_depth_without_base_counts_is_not_set_to_zero(self):
        matrix = self.matrix(with_counts=False)
        b = matrix.set_index("Sample_ID").loc["B"]
        self.assertEqual(b["Status"], COVERED_NOT_DETECTED)
        self.assertEqual(b["Callable_Evidence"], VALUE_UNKNOWN)
        self.assertTrue(math.isnan(b["Callable_Frequency"]))
        trend = frequency_trends(matrix).iloc[0]
        # Only the three reports: x = 0, 14, 21; y = 0.2, 0.4, 0.6. Mean x 35/3, mean y 0.4;
        # Sxy = 4.2, Sxx = 686/3 -> 0.0183673 per day = 0.128571 per week.
        self.assertAlmostEqual(trend["Slope_Per_Week"], 0.9 / 7, places=6)
        self.assertEqual(trend["Not_Reported_Value_Unknown"], 1)
        self.assertFalse(trend["Complete_Series"])
        self.assertIn("left out", trend["Trend_Basis"])

    def test_alternative_reads_below_the_report_give_a_count_frequency(self):
        matrix = self.matrix(with_counts=True, b_alt_reads=5)
        b = matrix.set_index("Sample_ID").loc["B"]
        self.assertAlmostEqual(b["Callable_Frequency"], 0.05)
        self.assertEqual(b["Alt_Read_Count"], 5)

    def test_low_coverage_is_missing(self):
        variants = reported_variant(["A", "C", "E"], [0.2, 0.4, 0.6])
        table = position_table([("A", 100, 20), ("B", 5, 0), ("C", 100, 40), ("E", 100, 60)])
        matrix = observation_matrix(variants, four_samples(), build_depth_lookup(table), 20, table)
        self.assertTrue(math.isnan(matrix.set_index("Sample_ID").loc["B", "Callable_Frequency"]))
        self.assertEqual(frequency_trends(matrix).iloc[0]["Low_Or_No_Coverage"], 1)

    def test_single_time_point_is_not_a_trend(self):
        matrix = pd.DataFrame({"Variant_Key": ["K"], "Location_Name": ["S"], "Sample_ID": ["A"],
                               "Collection_Date": DATES[:1], "Status": [VARIANT_OBSERVED],
                               "Alt_Frequency": [0.3]})
        trend = frequency_trends(matrix).iloc[0]
        self.assertTrue(math.isnan(trend["Change_In_Callable_Frequency"]))
        self.assertTrue(math.isnan(trend["Slope_Per_Week"]))
        self.assertIn("Fewer than two", trend["Trend_Basis"])


class SharedCoverageDistanceTest(unittest.TestCase):
    """1B / S3-03: no shared callable position means no distance, never 0."""

    @staticmethod
    def frequencies(*samples):
        """Each sample: list of per-position base frequency vectors, or None for no reads."""

        array = np.full((len(samples), len(samples[0]), 4), np.nan, dtype=np.float32)
        for index, positions in enumerate(samples):
            for position, vector in enumerate(positions):
                if vector is not None:
                    array[index, position] = vector
        return array

    def test_pair_without_shared_positions_is_undefined(self):
        a, g = [1, 0, 0, 0], [0, 0, 1, 0]
        data = self.frequencies([a, None], [None, a], [a, g])
        distances, compared = paper_distance_matrix(data)
        self.assertTrue(np.isnan(distances[0, 1]))  # old code returned 0 here
        self.assertEqual(compared[0, 1], 0)
        # Sample 1 vs 3 share position 1 only: identical A -> 0. Sample 2 vs 3 share
        # position 2: A vs G -> 1 - 0/(1 + 1 - 0) = 1.
        self.assertAlmostEqual(distances[0, 2], 0.0)
        self.assertAlmostEqual(distances[1, 2], 1.0)

    def test_minimum_shared_positions(self):
        a = [1, 0, 0, 0]
        data = self.frequencies([a, a], [a, None])
        distances, compared = paper_distance_matrix(data, minimum_shared_positions=2)
        self.assertEqual(compared[0, 1], 1)
        self.assertTrue(np.isnan(distances[0, 1]))

    def test_jaccard_of_two_empty_sets_is_undefined(self):
        distances = jaccard_distance_matrix([set(), set(), {"x"}])
        self.assertTrue(np.isnan(distances[0, 1]))
        self.assertAlmostEqual(distances[0, 2], 1.0)

    def samples(self, ids):
        return pd.DataFrame({"Sample_ID": ids, "Reference_ID": ["R.1"] * len(ids),
                             "Collection_Date": pd.date_range("2026-01-01", periods=len(ids)),
                             "Location_Name": ["S"] * len(ids), "Analysis_Included": [True] * len(ids)})

    @staticmethod
    def counts(rows):
        table = pd.DataFrame([{"Sample_ID": s, "Reference_ID": "R.1", "Position": p, "A_Count": a,
                               "C_Count": c, "G_Count": 0, "T_Count": 0, "Depth": a + c}
                              for s, p, a, c in rows])
        table.attrs["has_base_counts"] = True
        return table

    def test_pcoa_skipped_when_a_pair_has_no_shared_coverage(self):
        table = self.counts([("A", 1, 10, 0), ("B", 2, 10, 0), ("C", 1, 10, 0), ("C", 2, 0, 10)])
        variants = pd.DataFrame({"Sample_ID": ["A"], "Variant_Key": ["k"], "Analysis_Included": [True]})
        result = population_structure(self.samples(["A", "B", "C"]), variants, table)
        self.assertFalse(result["available"])
        self.assertIn("A/B", result["message"])
        self.assertIn("not replaced by 0", result["message"])
        pairs = result["pairs"].set_index(["Sample_ID_1", "Sample_ID_2"])
        self.assertFalse(pairs.loc[("A", "B"), "Supported"])
        self.assertEqual(pairs.loc[("A", "C"), "Shared_Callable_Positions"], 1)

    def test_sample_without_reads_is_left_out_and_pcoa_runs(self):
        table = self.counts([(s, p, 10, 0) for s in "ABC" for p in (1, 2)] + [("D", 1, 0, 0)])
        table.loc[(table["Sample_ID"] == "C") & (table["Position"] == 2), ["A_Count", "C_Count"]] = [0, 10]
        variants = pd.DataFrame({"Sample_ID": ["A"], "Variant_Key": ["k"], "Analysis_Included": [True]})
        result = population_structure(self.samples(["A", "B", "C", "D"]), variants, table)
        self.assertTrue(result["available"])
        self.assertEqual(result["left_out_no_data"], ["D"])
        self.assertEqual(list(result["coordinates"]["Sample_ID"]), ["A", "B", "C"])
        # A and B are identical at both positions; C differs completely at position 2.
        self.assertAlmostEqual(result["distances"].loc["A", "B"], 0.0)
        self.assertAlmostEqual(result["distances"].loc["A", "C"], 1.0)


class CollectionDateTest(unittest.TestCase):
    """1C / S3-05: day/month order is never guessed."""

    def test_parsing_rules(self):
        cases = {
            "2021-04-03": (pd.Timestamp("2021-04-03"), PARSED),
            "2020-04-23 00:00:00": (pd.Timestamp("2020-04-23"), PARSED),
            "13/04/2021": (pd.Timestamp("2021-04-13"), PARSED),
            "04/13/2021": (pd.Timestamp("2021-04-13"), PARSED),
            "03/03/2021": (pd.Timestamp("2021-03-03"), PARSED),
            "3 April 2021": (pd.Timestamp("2021-04-03"), PARSED),
            "Apr 3, 2021": (pd.Timestamp("2021-04-03"), PARSED),
            "03/04/2021": (None, AMBIGUOUS),
            "03/04/21": (None, TWO_DIGIT_YEAR),
            "2021-04": (None, PARTIAL),
            "31/02/2021": (None, INVALID),
            "13/14/2021": (None, INVALID),
            "soon": (None, INVALID),
            "": (None, NOT_SUPPLIED),
        }
        for text, expected in cases.items():
            with self.subTest(text=text):
                self.assertEqual(parse_one(text), expected)

    def test_selected_order_resolves_and_contradictions_are_invalid(self):
        self.assertEqual(parse_one("03/04/2021", "day_first"), (pd.Timestamp("2021-04-03"), PARSED_WITH_ORDER))
        self.assertEqual(parse_one("03/04/2021", "month_first"), (pd.Timestamp("2021-03-04"), PARSED_WITH_ORDER))
        self.assertEqual(parse_one("04/13/2021", "day_first"), (None, INVALID))
        self.assertEqual(parse_one("23_Apr_2020", fixed_format="%d_%b_%Y"), (pd.Timestamp("2020-04-23"), PARSED))

    def test_original_text_is_preserved(self):
        parsed = parse_collection_dates(pd.Series(["03/04/2021", "2021-04-05", ""]))
        self.assertEqual(list(parsed["Original"]), ["03/04/2021", "2021-04-05", ""])
        self.assertTrue(pd.isna(parsed["Date"].iloc[0]))
        self.assertEqual(parsed["Date"].iloc[1], pd.Timestamp("2021-04-05"))

    def snv_file(self, folder, dates):
        rows = [{"Sample_ID": f"S{index}", "Reference_Genome": "R.1", "Position": 10, "Ref": "A",
                 "Alt": "G", "Alt_Frequency": 0.5, "Read_Depth": 100, "Collection_Date": date}
                for index, date in enumerate(dates)]
        path = Path(folder) / "snvs.csv"
        pd.DataFrame(rows).to_csv(path, index=False)
        return path

    def test_loader_flags_ambiguous_dates_and_uses_a_selected_order(self):
        with tempfile.TemporaryDirectory() as folder:
            path = self.snv_file(folder, ["03/04/2021", "2021-04-20"])
            data, warnings = load_processed_snvs(path)
            first = data.set_index("Sample_ID").loc["S0"]
            self.assertTrue(pd.isna(first["Collection_Date"]))
            self.assertEqual(first["Collection_Date_Original"], "03/04/2021")
            self.assertEqual(first["Collection_Date_Status"], AMBIGUOUS)
            self.assertTrue(any("ambiguous" in warning for warning in warnings))
            data, _ = load_processed_snvs(path, date_order="day_first")
            self.assertEqual(data.set_index("Sample_ID").loc["S0", "Collection_Date"], pd.Timestamp("2021-04-03"))

    def test_loader_rejects_dates_that_are_not_dates(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(ValueError, "collection date"):
                load_processed_snvs(self.snv_file(folder, ["2021-02-30"]))


class MissingDateFigureTest(unittest.TestCase):
    """1D / S3-04: no dates -> a plain message, no figure error, other analyses still run."""

    def run_with_dates(self, transform):
        snvs = pd.read_csv(EXAMPLES / "synthetic_pooled_snvs.csv", dtype=str, keep_default_na=False)
        snvs["Collection_Date"] = transform(snvs["Collection_Date"])
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "snvs.csv"
            snvs.to_csv(path, index=False)
            return run_integrated_analysis(snv_path=path)

    def test_no_dates(self):
        result = self.run_with_dates(lambda dates: "")
        self.assertEqual(result["notes"].get("figures", []), [])
        self.assertNotIn("trajectories", result["figures"])
        self.assertTrue(any("Temporal trends unavailable because collection dates were not provided" in text
                            for text in result["notes"]["temporal"]))
        for name in ["variants_per_sample", "heatmap", "observation_status", "pcoa"]:
            self.assertIn(name, result["figures"])
        status = {row["Analysis"]: row["Status"] for row in result["analyses"]}
        self.assertEqual(status["Temporal frequency trends"], "Skipped")

    def test_partial_and_ambiguous_dates(self):
        result = self.run_with_dates(lambda dates: [value if index % 2 else "03/04/2026"
                                                    for index, value in enumerate(dates)])
        self.assertEqual(result["notes"].get("figures", []), [])
        self.assertTrue(any("ambiguous" in text for text in result["notes"]["quality"]))


class ReferenceVerificationTest(unittest.TestCase):
    """1E / S3-07: verified, mismatched, unverified, and missing reference labels."""

    reference = {"Reference_ID": "R.1", "Length": 10, "Sequence": "ACGTACGTAC", "Checksum": "x"}

    def variants(self, checks, label="R.1"):
        return pd.DataFrame({"Reference_ID": [label] * len(checks), "QC_Reference_Base": checks})

    def status(self, checks, reference=reference, label="R.1"):
        return reference_verification(self.variants(checks, label), reference, {}).iloc[0]["Status"]

    def test_statuses(self):
        self.assertEqual(self.status(["Pass", "Pass"]), REFERENCE_VERIFIED)
        self.assertEqual(self.status(["Pass", "Pass", "Fail (Ref base does not match the reference)"]),
                         REFERENCE_PARTLY_MISMATCHED)
        self.assertEqual(self.status(["Pass", "Fail (Ref base does not match the reference)",
                                      "Fail (position beyond the reference)"]), REFERENCE_MISMATCH)
        self.assertEqual(self.status(["Not checked (no reference supplied)"], reference=None),
                         REFERENCE_UNVERIFIED)
        self.assertEqual(self.status(["Not checked (record uses a different reference)"], label="Other.1"),
                         REFERENCE_UNVERIFIED)

    def test_pipeline_labels_pooled_only_runs_unverified(self):
        result = run_integrated_analysis(snv_path=EXAMPLES / "synthetic_pooled_snvs.csv")
        table = result["tables"]["reference_verification.csv"].set_index("Reference_ID")
        self.assertEqual(table.loc["SYNREF.1", "Status"], REFERENCE_UNVERIFIED)
        with_reference = run_integrated_analysis(snv_path=EXAMPLES / "synthetic_pooled_snvs.csv",
                                                 reference_path=EXAMPLES / "synthetic_reference.fasta")
        table = with_reference["tables"]["reference_verification.csv"].set_index("Reference_ID")
        # The example has one deliberate Ref mismatch at position 110 among 21 SYNREF.1 records.
        self.assertEqual(table.loc["SYNREF.1", "Status"], REFERENCE_PARTLY_MISMATCHED)
        self.assertEqual(table.loc["SYNREF.1", "Ref_Base_Matches"], 20)
        self.assertEqual(table.loc["SYNREF_B.1", "Status"], REFERENCE_UNVERIFIED)


class FrequencyAndDuplicateTest(unittest.TestCase):
    """1F: reported vs count-derived frequency (S3-06) and duplicate-aware medians (S3-13)."""

    def test_reported_frequency_is_kept_and_disagreement_flagged(self):
        table = pd.DataFrame({"Alt_Frequency": [0.30, 0.333, 0.5],
                              "Alt_Depth": pd.array([40, 1, pd.NA], dtype="Int64"),
                              "Read_Depth": pd.array([100, 3, 50], dtype="Int64")})
        checked = add_frequency_checks(table)
        self.assertEqual(list(checked["Alt_Frequency"]), [0.30, 0.333, 0.5])  # never overwritten
        self.assertAlmostEqual(checked["Count_Derived_Frequency"][0], 0.40)
        self.assertAlmostEqual(checked["Frequency_Difference"][0], -0.10)
        self.assertIn("Differs", checked["Frequency_Check"][0])
        self.assertEqual(checked["Frequency_Check"][1], "Agrees with alternate reads / depth")  # 1/3 vs 0.333
        self.assertEqual(checked["Frequency_Check"][2], "No alternate read count supplied")

    def duplicated_variants(self):
        return pd.DataFrame({
            "Variant_Key": ["R.1:5:A>G"] * 4, "Reference_ID": ["R.1"] * 4, "Position": [5] * 4,
            "Ref": ["A"] * 4, "Alt": ["G"] * 4, "Sample_ID": ["S1", "S1", "S2", "S2"],
            "Alt_Frequency": [0.2, 0.6, 0.4, 0.4], "Read_Depth": [100, 100, 80, 80],
            "Alt_Depth": pd.array([pd.NA] * 4, dtype="Int64"), "Analysis_Included": [True] * 4,
            "Collection_Date": pd.to_datetime(["2026-01-01"] * 4), "Location_Name": ["Site"] * 4,
        })

    def test_repeated_records_are_classified(self):
        marked = mark_repeated_records(self.duplicated_variants())
        self.assertEqual(list(marked["Counted_Once"]), [True, False, True, False])
        self.assertEqual(marked["Duplicate_Status"][0], REPEAT_CONFLICTING)
        self.assertEqual(marked["Duplicate_Status"][2], REPEAT_IDENTICAL)
        self.assertTrue(marked["Duplicate_Status"][3].endswith("(not counted again)"))
        unique = mark_repeated_records(self.duplicated_variants().iloc[[0, 2]].copy())
        self.assertEqual(set(unique["Duplicate_Status"]), {UNIQUE_RECORD})

    def test_medians_count_each_sample_once(self):
        variants = self.duplicated_variants()
        # S1 first row 0.2, S2 0.4 -> median 0.3 (the old median over all rows was 0.4).
        summary = position_summary(variants, None, "R.1")
        self.assertAlmostEqual(summary.iloc[0]["Median_Reported_Frequency"], 0.3)
        samples = pd.DataFrame({"Sample_ID": ["S1", "S2"], "Location_Name": ["Site"] * 2,
                                "Collection_Date": pd.to_datetime(["2026-01-01"] * 2),
                                "Analysis_Included": [True, True]})
        site = location_date_summary(variants, samples).iloc[0]
        self.assertAlmostEqual(site["Median_Reported_Frequency"], 0.3)
        self.assertEqual(site["Variant_Records"], 4)  # records are still counted as records


if __name__ == "__main__":
    unittest.main()
