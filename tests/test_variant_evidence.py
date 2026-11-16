"""Sprint 5 data model and depth-aware VAF uncertainty.

Expected values come from published tables (Newcombe 1998, Table I, Wilson
score method) or from an independent bisection that inverts the binomial
score test, not from the closed-form formula used in the code.
"""

import math
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from genomics.frequency_uncertainty import Z_95, depth_dependence_table, wilson_interval
from genomics.pooled_observations import pooled_variant_records
from genomics.quality_control import apply_pooled_quality_control, reference_verification
from genomics.variant_evidence import (
    BASE_COUNTS,
    FREQUENCY_ONLY,
    INTERVAL_RECONSTRUCTED,
    INTERVAL_SUPPLIED,
    NO_DEPTH,
    NO_INTERVAL_CONFLICT,
    NO_INTERVAL_DEPTH,
    NO_INTERVAL_DISAGREE,
    NO_INTERVAL_REPEAT,
    NO_INTERVAL_SITE_SUM,
    RECONSTRUCTED,
    SUPPLIED,
    build_variant_evidence,
    reconstruction_check,
    validate_evidence,
)
from wastewater.processed_snv import load_processed_snvs


def score_test_bound(x, n, upper, z=Z_95):
    """Independent Wilson bound: solve |x/n - p| = z sqrt(p(1-p)/n) by bisection."""

    p_hat = x / n
    if (upper and x == n) or (not upper and x == 0):
        return p_hat
    low, high = (p_hat, 1.0) if upper else (0.0, p_hat)
    for _ in range(200):
        middle = (low + high) / 2
        inside = abs(p_hat - middle) <= z * math.sqrt(middle * (1 - middle) / n)
        if upper:
            low, high = (middle, high) if inside else (low, middle)
        else:
            low, high = (low, middle) if inside else (middle, high)
    return (low + high) / 2


class WilsonIntervalTest(unittest.TestCase):
    def test_published_values(self):
        # Newcombe (1998) Stat Med 17:857-872, Table I, Wilson score interval (method 3).
        published = {(81, 263): (0.2553, 0.3662), (15, 148): (0.0624, 0.1605),
                     (0, 20): (0.0000, 0.1611), (1, 29): (0.0061, 0.1718)}
        for (x, n), (low, high) in published.items():
            with self.subTest(x=x, n=n):
                lower, upper = wilson_interval([x], [n])
                self.assertAlmostEqual(lower[0], low, places=4)
                self.assertAlmostEqual(upper[0], high, places=4)

    def test_matches_score_test_inversion(self):
        for x, n in [(4, 20), (40, 200), (3, 7), (999, 1000), (1, 12317)]:
            with self.subTest(x=x, n=n):
                lower, upper = wilson_interval([x], [n])
                self.assertAlmostEqual(lower[0], score_test_bound(x, n, False), places=9)
                self.assertAlmostEqual(upper[0], score_test_bound(x, n, True), places=9)

    def test_same_vaf_narrows_with_depth(self):
        # 4/20 and 40/200 are both 0.20; the deeper one must be narrower.
        lower, upper = wilson_interval([4, 40], [20, 200])
        self.assertAlmostEqual(lower[0], 0.080658, places=6)
        self.assertAlmostEqual(upper[0], 0.416017, places=6)
        self.assertAlmostEqual(lower[1], 0.150452, places=6)
        self.assertAlmostEqual(upper[1], 0.260855, places=6)
        self.assertLess(upper[1] - lower[1], upper[0] - lower[0])
        table = depth_dependence_table()
        self.assertEqual(list(table["Depth"]), [20, 50, 100, 300, 1000])
        self.assertTrue((table["VAF"] == 0.2).all())
        self.assertTrue(table["Interval_Width"].is_monotonic_decreasing)

    def test_boundaries(self):
        # x = 0: lower 0 and upper z^2 / (n + z^2); x = n: upper 1 and lower n / (n + z^2).
        lower, upper = wilson_interval([0, 50], [50, 50])
        z2 = Z_95 ** 2
        self.assertEqual(lower[0], 0.0)
        self.assertAlmostEqual(upper[0], z2 / (50 + z2))
        self.assertEqual(upper[1], 1.0)
        self.assertAlmostEqual(lower[1], 50 / (50 + z2))
        self.assertGreater(upper[0], 0.0)  # x = 0 still has a positive upper bound

    def test_invalid_counts_give_no_interval(self):
        lower, upper = wilson_interval([5, 3, np.nan, 2.5, -1, 1], [4, 0, 10, 10, 10, np.nan])
        self.assertTrue(np.isnan(lower).all() and np.isnan(upper).all())


class ReconstructionTest(unittest.TestCase):
    def test_unique_reconstruction_at_table_precision(self):
        # The table writes up to 6 significant digits, so '0.14' means 0.140000:
        # 7/50 = 0.14 exactly; 6/50 and 8/50 do not round to it. 5/103 = 0.0485437 (6 digits).
        counts, details = reconstruction_check([0.14, 0.0485437], [50, 103], [2, 6])
        self.assertEqual(list(counts), [7, 5])
        self.assertTrue(details["used"])
        self.assertEqual(details["precision"], 6)

    def test_ambiguous_count_is_not_reconstructed(self):
        # At 2-digit precision 139/1000 and 141/1000 also round to 0.14.
        counts, details = reconstruction_check([0.14, 0.5], [1000, 10], [2, 1])
        self.assertTrue(np.isnan(counts[0]))
        self.assertEqual(counts[1], 5)
        self.assertEqual(details["ambiguous"], 1)

    def test_inconsistent_table_reconstructs_nothing(self):
        # 0.3 x 7 = 2.1 reads: no whole count reproduces 0.3 at depth 7, so the
        # frequency and depth do not share a denominator for this table.
        counts, details = reconstruction_check([0.3, 0.25], [7, 40], [1, 2])
        self.assertTrue(np.isnan(counts).all())
        self.assertFalse(details["used"])
        self.assertEqual(details["not_reproduced"], 1)


def evidence_from_rows(rows, position_table=None, reference=None):
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / "snvs.csv"
        pd.DataFrame(rows).to_csv(path, index=False)
        data, _ = load_processed_snvs(path)
    variants = pooled_variant_records(data, "test table")
    samples = pd.DataFrame({"Sample_ID": sorted(variants["Sample_ID"].unique())})
    for column in ["Genome_Coverage_Percent"]:
        samples[column] = np.nan
    variants, samples, _ = apply_pooled_quality_control(variants, samples, reference=reference)
    status = reference_verification(variants, reference, {})
    return build_variant_evidence(variants, position_table, status, reference)


def row(sample, position, ref, alt, frequency, depth, alt_depth="", **extra):
    return {"Sample_ID": sample, "Reference_Genome": "R.1", "Position": position, "Ref": ref, "Alt": alt,
            "Alt_Frequency": frequency, "Read_Depth": depth, "Alt_Depth": alt_depth,
            "Collection_Date": "2026-01-01", **extra}


class VariantEvidenceTest(unittest.TestCase):
    def test_direct_counts(self):
        evidence, _ = evidence_from_rows([row("S1", 10, "A", "G", "0.2", 200, 40)])
        record = evidence.iloc[0]
        self.assertEqual(record["Alt_Count_Source"], SUPPLIED)
        self.assertEqual(record["Alt_Read_Count"], 40)
        self.assertEqual(record["Total_Depth"], 200)
        self.assertAlmostEqual(record["VAF"], 0.2)
        self.assertEqual(record["Interval_Status"], INTERVAL_SUPPLIED)
        self.assertAlmostEqual(record["Wilson_Lower"], score_test_bound(40, 200, False), places=9)
        self.assertAlmostEqual(record["Wilson_Upper"], score_test_bound(40, 200, True), places=9)
        self.assertTrue(pd.isna(record["Ref_Read_Count"]))  # never derived as 200 - 40
        self.assertEqual(validate_evidence(evidence), [])

    def test_supplied_reference_count_is_kept(self):
        evidence, _ = evidence_from_rows([row("S1", 10, "A", "G", "0.2", 200, 40, Ref_Depth=150)])
        self.assertEqual(evidence.iloc[0]["Ref_Read_Count"], 150)

    def test_reported_and_count_frequency_disagree(self):
        evidence, _ = evidence_from_rows([row("S1", 10, "A", "G", "0.3", 100, 40)])
        record = evidence.iloc[0]
        self.assertEqual(record["Reported_Frequency"], 0.3)
        self.assertAlmostEqual(record["Count_Derived_Frequency"], 0.4)
        self.assertEqual(record["Interval_Status"], NO_INTERVAL_DISAGREE)
        self.assertTrue(math.isnan(record["Wilson_Lower"]))
        self.assertEqual(record["VAF"], 0.3)  # the measured source value, not overwritten

    def test_reconstructed_counts_are_labeled_approximate(self):
        evidence, details = evidence_from_rows([row("S1", 10, "A", "G", "0.25", 40),
                                                row("S2", 10, "A", "G", "0.0485437", 103)])
        self.assertTrue(details["used"])
        self.assertEqual(list(evidence["Alt_Count_Source"]), [RECONSTRUCTED, RECONSTRUCTED])
        self.assertEqual(list(evidence["Alt_Read_Count"]), [10, 5])
        self.assertTrue((evidence["Interval_Status"] == INTERVAL_RECONSTRUCTED).all())
        self.assertIn("counts reconstructed", evidence.iloc[0]["Validation_Warnings"])

    def test_frequency_only_when_reconstruction_fails(self):
        evidence, details = evidence_from_rows([row("S1", 10, "A", "G", "0.3", 7),
                                                row("S2", 10, "A", "G", "0.25", 40)])
        self.assertFalse(details["used"])
        self.assertTrue((evidence["Alt_Count_Source"] == FREQUENCY_ONLY).all())
        self.assertTrue(evidence["Wilson_Lower"].isna().all())
        self.assertTrue(evidence["Alt_Read_Count"].isna().all())

    def test_missing_depth_is_not_zero_depth(self):
        evidence, _ = evidence_from_rows([row("S1", 10, "A", "G", "0.25", "")])
        record = evidence.iloc[0]
        self.assertEqual(record["Alt_Count_Source"], NO_DEPTH)
        self.assertEqual(record["Interval_Status"], NO_INTERVAL_DEPTH)
        self.assertTrue(pd.isna(record["Total_Depth"]))
        self.assertIn("depth not supplied", record["Coverage_Status"])

    def test_zero_depth_and_invalid_counts_are_rejected_at_import(self):
        for bad in [row("S1", 10, "A", "G", "0.25", 0), row("S1", 10, "A", "G", "0.25", 10, 11),
                    row("S1", 10, "A", "G", "0.25", 10, 2.5), row("S1", 10, "A", "G", "0.25", 10, -1)]:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                evidence_from_rows([bad])

    def test_multiallelic_sites(self):
        evidence, _ = evidence_from_rows([row("S1", 10, "A", "G", "0.3", 100, 30),
                                          row("S1", 10, "A", "T", "0.2", 100, 20),
                                          row("S2", 10, "A", "G", "0.6", 100, 60),
                                          row("S2", 10, "A", "T", "0.5", 100, 50)])
        by_sample = evidence.groupby("Sample_ID")
        self.assertEqual(list(evidence["Alleles_At_Site"]), [2, 2, 2, 2])
        self.assertTrue((by_sample.get_group("S1")["Interval_Status"] == INTERVAL_SUPPLIED).all())
        # 0.6 + 0.5 > 1: the allele model is inconsistent, so no interval is made.
        self.assertTrue((by_sample.get_group("S2")["Interval_Status"] == NO_INTERVAL_SITE_SUM).all())

    def test_duplicate_records(self):
        evidence, _ = evidence_from_rows([row("S1", 10, "A", "G", "0.2", 100, 20),
                                          row("S1", 10, "A", "G", "0.2", 100, 20),
                                          row("S2", 10, "A", "G", "0.2", 100, 20),
                                          row("S2", 10, "A", "G", "0.4", 100, 40)])
        status = list(evidence["Interval_Status"])
        self.assertEqual(status[0], INTERVAL_SUPPLIED)
        self.assertEqual(status[1], NO_INTERVAL_REPEAT)
        self.assertEqual(status[2], NO_INTERVAL_CONFLICT)
        self.assertEqual(status[3], NO_INTERVAL_REPEAT)

    def test_base_counts(self):
        table = pd.DataFrame([{"Sample_ID": "S1", "Reference_ID": "R.1", "Position": 10, "Depth": 120,
                               "A_Count": 90, "C_Count": 0, "G_Count": 30, "T_Count": 0}])
        table.attrs["has_base_counts"] = True
        evidence, _ = evidence_from_rows([row("S1", 10, "A", "G", "0.25", 120)], position_table=table)
        record = evidence.iloc[0]
        self.assertEqual(record["Alt_Count_Source"], BASE_COUNTS)
        self.assertEqual(record["Alt_Read_Count"], 30)
        self.assertEqual(record["Ref_Read_Count"], 90)
        self.assertAlmostEqual(record["VAF"], 0.25)

    def test_unverified_reference_is_labeled(self):
        evidence, _ = evidence_from_rows([row("S1", 10, "A", "G", "0.2", 200, 40)])
        self.assertEqual(evidence.iloc[0]["Reference_Status"], "Supplied but unverified")
        self.assertIn("reference not verified", evidence.iloc[0]["Validation_Warnings"])

    @unittest.skipUnless(os.environ.get("FONTENELE_MMC2"), "set FONTENELE_MMC2 to the published mmc2.xlsx")
    def test_published_workbook_counts_are_reconstructed(self):
        data, _ = load_processed_snvs(os.environ["FONTENELE_MMC2"])
        variants = pooled_variant_records(data, "mmc2")
        samples = pd.DataFrame({"Sample_ID": sorted(variants["Sample_ID"].unique()),
                                "Genome_Coverage_Percent": np.nan})
        variants, _, _ = apply_pooled_quality_control(variants, samples)
        evidence, details = build_variant_evidence(variants, None, reference_verification(variants))
        self.assertTrue(details["used"])
        self.assertEqual(details["precision"], 6)
        self.assertEqual(int((evidence["Alt_Count_Source"] == RECONSTRUCTED).sum()), 7973)
        self.assertEqual(validate_evidence(evidence), [])
        product = evidence["Reported_Frequency"] * evidence["Total_Depth"].astype(float)
        self.assertLess(float((product - evidence["Alt_Read_Count"].astype(float)).abs().max()), 0.5)


if __name__ == "__main__":
    unittest.main()
