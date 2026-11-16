"""Sprint 5 Phase 3: numerical checks of mutation profiles and exploratory grouping.

Small deterministic inputs with known answers. Optional cross-checks against
SciPy run only where SciPy is installed (it is not a project dependency).
"""

import math
import unittest

import numpy as np
import pandas as pd

from genomics import chi_square as cs
from genomics import mutation_grouping as mg
from genomics.frequency_uncertainty import wilson_interval
from genomics.mutation_profiles import (
    FREQUENCY_ONLY,
    MEASURED,
    REPLICATE_UNIT,
    build_mutation_profiles,
    effective_score_interval,
    variance_factor,
)
from genomics.temporal_analysis import ZERO_ALT_READS
from genomics.unified_schema import COVERAGE_UNKNOWN, COVERED_NOT_DETECTED, VARIANT_OBSERVED

try:
    from scipy import stats
except ImportError:  # SciPy is optional
    stats = None


def arrays_from(F, N):
    F, N = np.asarray(F, float), np.asarray(N, float)
    return {"F": F, "X": np.round(F * N), "N": N, "Q": N ** 2}


class ChiSquareTest(unittest.TestCase):
    def test_closed_forms(self):
        for x in [0.1, 1.0, 5.0, 30.0]:
            self.assertAlmostEqual(cs.chi2_sf(x, 2), math.exp(-x / 2), places=12)
            self.assertAlmostEqual(cs.chi2_sf(x, 1), math.erfc(math.sqrt(x / 2)), places=12)
        self.assertAlmostEqual(cs.chi2_isf(0.05, 1), 3.841458820694124, places=8)
        self.assertAlmostEqual(cs.ncx2_sf(4.0, 3, 0.0), cs.chi2_sf(4.0, 3), places=14)

    def test_power_noncentrality_gives_the_requested_power(self):
        for k in [1, 4, 12]:
            lam = cs.noncentrality_for_power(k, 0.01, 0.8)
            self.assertAlmostEqual(cs.ncx2_sf(cs.chi2_isf(0.01, k), k, lam), 0.8, places=6)

    @unittest.skipIf(stats is None, "SciPy not installed (optional cross-check)")
    def test_against_scipy(self):
        for x, k in [(0.5, 1), (100, 80), (250, 200), (0.01, 0.5)]:
            self.assertAlmostEqual(cs.chi2_sf(x, k), stats.chi2.sf(x, k), places=11)
        for x, k, lam in [(6.63, 1, 5), (20, 5, 10), (50, 30, 3)]:
            self.assertAlmostEqual(cs.ncx2_sf(x, k, lam), stats.ncx2.sf(x, k, lam), places=10)


class DistanceTest(unittest.TestCase):
    def test_ordinary_distance_uses_shared_units_only_and_never_fills_zeros(self):
        F = np.array([[0.1, 0.5, np.nan, 0.3],
                      [0.2, 0.5, 0.9, np.nan],
                      [np.nan, np.nan, 0.9, 0.3]])
        result = mg.ordinary_distances(F, ["R", "R", "R"], minimum_shared=2)
        self.assertAlmostEqual(result["distance"][0, 1], math.sqrt((0.1 ** 2 + 0) / 2))
        # Rows 0 and 2 share one unit only: no distance (inf), not 0 or a filled-in value.
        self.assertEqual(result["shared"][0, 2], 1)
        self.assertTrue(np.isinf(result["distance"][0, 2]))

    def test_different_references_are_never_compared(self):
        F = np.tile([0.1, 0.4, 0.7], (2, 1))
        result = mg.ordinary_distances(F, ["SEG1", "SEG2"], minimum_shared=2)
        self.assertTrue(np.isinf(result["distance"][0, 1]))

    def test_pair_statistic_is_pearson_chi_square_at_rho_zero(self):
        x = np.array([[30, 120], [45, 100]], float)
        n = np.array([[200, 400], [210, 380]], float)
        arrays = {"F": x / n, "X": x, "N": n, "Q": n ** 2}
        result = mg.depth_aware_distances(arrays, ["R", "R"], np.zeros(2), minimum_shared=2)
        expected = 0.0
        for s in range(2):
            table = np.array([[x[0, s], n[0, s] - x[0, s]], [x[1, s], n[1, s] - x[1, s]]])
            expect = table.sum(axis=1, keepdims=True) * table.sum(axis=0, keepdims=True) / table.sum()
            expected += ((table - expect) ** 2 / expect).sum()
        self.assertAlmostEqual(result["chi2"][0, 1], expected, places=9)
        self.assertEqual(result["df"][0, 1], 2)

    def test_power_gate_blocks_low_depth_pairs(self):
        F = np.array([[0.2, 0.5, 0.8], [0.2, 0.5, 0.8]])
        low = mg.depth_aware_distances(arrays_from(F, np.full_like(F, 10)), ["R", "R"], np.zeros(3), 2)
        high = mg.depth_aware_distances(arrays_from(F, np.full_like(F, 1000)), ["R", "R"], np.zeros(3), 2)
        self.assertFalse(low["comparable"][0, 1])
        self.assertTrue(np.isinf(low["distance"][0, 1]))
        self.assertTrue(high["comparable"][0, 1])

    def test_variance_factor_matches_phase2_and_combines_libraries(self):
        self.assertAlmostEqual(float(variance_factor(300, 300 ** 2, 0.02)), (1 + 299 * 0.02) / 300)
        self.assertAlmostEqual(float(variance_factor(300, 300 ** 2, 0.0)), 1 / 300)
        # Two libraries of 100 and 200 reads: Var = p(1-p)[100(1+99r) + 200(1+199r)] / 300^2.
        r = 0.01
        self.assertAlmostEqual(float(variance_factor(300, 100 ** 2 + 200 ** 2, r)),
                               (100 * (1 + 99 * r) + 200 * (1 + 199 * r)) / 300 ** 2)


class LinkageTest(unittest.TestCase):
    def test_complete_linkage_heights_and_forest(self):
        D = np.array([[0, 1, 4, np.inf],
                      [1, 0, 3, np.inf],
                      [4, 3, 0, np.inf],
                      [np.inf, np.inf, np.inf, 0]], float)
        merges = mg.complete_linkage(D)
        self.assertEqual([(a, b, h) for a, b, h, _ in merges], [(0, 1, 1.0), (4, 2, 4.0)])  # max(4, 3)
        _, _, _, roots = mg.tree_members(merges, 4)
        self.assertEqual(sorted(roots), [3, 5])  # mutation 3 is never joined across an infinite distance

    @unittest.skipIf(stats is None, "SciPy not installed (optional cross-check)")
    def test_matches_scipy_complete_linkage(self):
        from scipy.cluster.hierarchy import linkage
        from scipy.spatial.distance import squareform
        rng = np.random.default_rng(3)
        points = rng.random((12, 3))
        D = np.sqrt(((points[:, None] - points[None]) ** 2).sum(-1))
        ours = sorted(h for _, _, h, _ in mg.complete_linkage(D))
        theirs = sorted(linkage(squareform(D, checks=False), "complete")[:, 2])
        np.testing.assert_allclose(ours, theirs)

    def test_homogeneity_test_is_pearson_k_by_2(self):
        x = np.array([[30.0], [45.0], [50.0]])
        n = np.array([[200.0], [210.0], [260.0]])
        statistic, df, p_value = mg.homogeneity_test([0, 1, 2], {"F": x / n, "X": x, "N": n, "Q": n ** 2},
                                                     np.zeros(1))
        table = np.hstack([x, n - x])
        expect = table.sum(axis=1, keepdims=True) * table.sum(axis=0, keepdims=True) / table.sum()
        self.assertAlmostEqual(statistic, ((table - expect) ** 2 / expect).sum(), places=9)
        self.assertEqual(df, 2)
        self.assertAlmostEqual(p_value, cs.chi2_sf(statistic, 2))


class GroupingBehaviourTest(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(11)
        units = 10
        self.patterns = [np.linspace(0.1, 0.8, units), np.linspace(0.7, 0.15, units),
                         0.4 + 0.3 * np.sin(np.arange(units))]
        rows, truth = [], []
        for label, pattern in enumerate(self.patterns):
            for _ in range(4):
                rows.append(rng.binomial(500, pattern) / 500)
                truth.append(label)
        self.F, self.truth = np.array(rows), np.array(truth)
        self.N = np.full_like(self.F, 500)

    def test_identical_profiles_group_and_distinct_profiles_do_not(self):
        arrays = arrays_from(self.F, self.N)
        for method in ["ordinary", "depth"]:
            result = mg.run_method(method, arrays, ["R"] * len(self.F), np.zeros(10), 2)
            labels = mg.label_vector(result["groups"], len(self.F))
            self.assertAlmostEqual(mg.adjusted_rand_index(self.truth, labels), 1.0, msg=method)

    def test_flat_unrelated_profiles_are_not_grouped_by_the_depth_aware_method(self):
        rng = np.random.default_rng(5)
        level = np.clip(0.3 + rng.normal(0, 0.03, (6, 10)), 0.01, 0.99)
        arrays = arrays_from(np.round(level * 500) / 500, np.full((6, 10), 500))
        result = mg.run_method("depth", arrays, ["R"] * 6, np.zeros(10), 2)
        self.assertEqual(result["groups"], [])
        ordinary = mg.run_method("ordinary", arrays, ["R"] * 6, np.zeros(10), 2)
        self.assertTrue(ordinary["groups"])  # the depth-blind baseline joins them (a documented weakness)

    def test_single_unit_gives_no_comparable_pair(self):
        arrays = arrays_from(self.F[:, :1], self.N[:, :1])
        result = mg.run_method("depth", arrays, ["R"] * len(self.F), np.zeros(1), 2)
        self.assertEqual(result["groups"], [])
        self.assertFalse(np.isfinite(result["pairs"]["distance"][~np.eye(len(self.F), dtype=bool)]).any())

    def test_assumed_rho_changes_what_counts_as_similar_or_comparable(self):
        # Two patterns 0.03 apart at 2000x. Binomial: the difference is detectable, so not grouped.
        F = np.vstack([np.linspace(0.2, 0.6, 8), np.linspace(0.2, 0.6, 8) + 0.03])
        arrays = arrays_from(F, np.full_like(F, 2000))
        binomial = mg.run_method("depth", arrays, ["R", "R"], np.zeros(8), 2)
        self.assertEqual(binomial["groups"], [])
        # rho 0.002: the measurement is noisier, the 0.03 gap is within noise, so they are grouped.
        moderate = mg.run_method("depth", arrays, ["R", "R"], np.full(8, 0.002), 2)
        self.assertEqual(moderate["groups"], [[0, 1]])
        # rho 0.05: too noisy to resolve a 0.10 difference -> not comparable, abstain (never grouped).
        large = mg.run_method("depth", arrays, ["R", "R"], np.full(8, 0.05), 2)
        self.assertEqual(large["groups"], [])
        self.assertFalse(large["pairs"]["comparable"][0, 1])


class MetricTest(unittest.TestCase):
    def test_known_values(self):
        self.assertAlmostEqual(mg.adjusted_rand_index([0, 0, 1, 1], [1, 1, 0, 0]), 1.0)
        # scikit-learn: adjusted_rand_score([0,0,1,1],[0,0,1,2]) = 0.5714285714
        self.assertAlmostEqual(mg.adjusted_rand_index([0, 0, 1, 1], [0, 0, 1, 2]), 4 / 7)
        # scikit-learn: normalized_mutual_info_score([0,0,1,1],[0,0,1,2]) = 0.8
        self.assertAlmostEqual(mg.normalized_mutual_information([0, 0, 1, 1], [0, 0, 1, 2]), 0.8)
        precision, recall = mg.pairwise_precision_recall([0, 0, 1, 1], [0, 0, 0, 1])
        self.assertAlmostEqual(precision, 1 / 3)
        self.assertAlmostEqual(recall, 1 / 2)
        self.assertTrue(math.isnan(mg.pairwise_precision_recall([0, 1], [0, 1])[0]))


def tiny_inputs():
    """Two samples x two mutations, plus a verified replicate pair, with all matrix states."""

    samples = pd.DataFrame({"Sample_ID": ["S1", "S2", "L1", "L2"],
                            "Collection_Date": pd.to_datetime(["2026-01-01", "2026-01-08", "2026-01-15", "2026-01-15"]),
                            "Location_Name": ["Site"] * 4, "Analysis_Included": True})
    key_a, key_b = "R:10:A>G", "R:20:C>T"
    matrix = pd.DataFrame([
        (key_a, "S1", VARIANT_OBSERVED, "Reported frequency", 0.25, 40, 160),
        (key_a, "S2", COVERED_NOT_DETECTED, ZERO_ALT_READS, 0.0, 0, 150),
        (key_a, "L1", VARIANT_OBSERVED, "Reported frequency", 0.30, 30, 100),
        (key_a, "L2", VARIANT_OBSERVED, "Reported frequency", 0.35, 70, 200),
        (key_b, "S1", COVERAGE_UNKNOWN, COVERAGE_UNKNOWN, np.nan, pd.NA, None),
        (key_b, "S2", VARIANT_OBSERVED, "Reported frequency", 0.5, 50, 100),
        (key_b, "L1", VARIANT_OBSERVED, "Reported frequency", 0.4, np.nan, np.nan),
        (key_b, "L2", COVERAGE_UNKNOWN, COVERAGE_UNKNOWN, np.nan, pd.NA, None),
    ], columns=["Variant_Key", "Sample_ID", "Status", "Callable_Evidence", "Callable_Frequency", "Alt_Read_Count",
                "Depth"])
    matrix = matrix.merge(samples[["Sample_ID", "Collection_Date", "Location_Name"]], on="Sample_ID")
    parts = matrix["Variant_Key"].str.extract(r"(?P<Reference_ID>[^:]+):(?P<Position>\d+):(?P<Ref>\w)>(?P<Alt>\w)")
    matrix = pd.concat([matrix, parts], axis=1)
    matrix["Position"] = matrix["Position"].astype(int)
    rows = []
    for sample, key, alt, depth, vaf, status in [
            ("S1", key_a, 40, 160, 0.25, "Calculated (supplied counts)"),
            ("L1", key_a, 30, 100, 0.30, "Calculated (supplied counts)"),
            ("L2", key_a, 70, 200, 0.35, "Calculated (supplied counts)"),
            ("S2", key_b, 50, 100, 0.5, "Calculated (supplied counts)"),
            ("L1", key_b, pd.NA, pd.NA, 0.4, "Not calculated: no read counts")]:
        rows.append({"Sample_ID": sample, "Variant_Key": key, "VAF": vaf, "Alt_Read_Count": alt, "Total_Depth": depth,
                     "Alt_Count_Source": "Supplied" if alt is not pd.NA else "Frequency only",
                     "Interval_Status": status, "Counted_Once": True, "Analysis_Included": True,
                     "Collection_Date": samples.set_index("Sample_ID").loc[sample, "Collection_Date"],
                     "Location_Name": "Site", "Reference_ID": "R", "Reference_Status": "Verified compatible",
                     "Annotated_Region": "", "Specimen_ID": "SPEC" if sample.startswith("L") else "",
                     "Library_ID": "LIB_" + sample})
    evidence = pd.DataFrame(rows)
    evidence["Alt_Read_Count"] = evidence["Alt_Read_Count"].astype("Int64")
    evidence["Total_Depth"] = evidence["Total_Depth"].astype("Int64")
    return matrix, evidence, samples


class ProfileTest(unittest.TestCase):
    def setUp(self):
        self.matrix, self.evidence, self.samples = tiny_inputs()
        self.built = build_mutation_profiles(self.matrix, self.evidence, None, self.samples, minimum_shared=2)
        self.profiles = self.built["profiles"].set_index(["Variant_Key", "Unit_ID"])

    def test_verified_replicates_become_one_unit_with_summed_counts(self):
        row = self.profiles.loc[("R:10:A>G", "SPEC [2 libraries]")]
        self.assertEqual(row["Unit_Type"], REPLICATE_UNIT)
        self.assertEqual(int(row["Alt_Read_Count"]), 100)
        self.assertEqual(int(row["Total_Depth"]), 300)
        self.assertEqual(row["Depth_Square_Sum"], 100 ** 2 + 200 ** 2)
        self.assertAlmostEqual(row["VAF"], 100 / 300)

    def test_unknown_is_never_zero_and_supported_zero_is_kept(self):
        self.assertNotIn(("R:20:C>T", "S1"), self.profiles.index)  # coverage unknown: no row, no value
        zero = self.profiles.loc[("R:10:A>G", "S2")]
        self.assertEqual(zero["VAF"], 0.0)
        self.assertEqual(zero["Count_Basis"], MEASURED)

    def test_frequency_only_values_are_labeled_and_excluded_from_depth_aware(self):
        row = self.profiles.loc[("R:20:C>T", "SPEC [2 libraries]")]
        self.assertEqual(row["Count_Basis"], FREQUENCY_ONLY)
        self.assertTrue(pd.isna(row["Alt_Read_Count"]))
        eligibility = self.built["eligibility"].set_index("Variant_Key")
        self.assertTrue(eligibility.loc["R:20:C>T", "Eligible_Ordinary"])
        self.assertFalse(eligibility.loc["R:20:C>T", "Eligible_Depth_Aware"])
        self.assertIn("read counts in 1 unit", eligibility.loc["R:20:C>T", "Exclusion_Reason"])
        self.assertEqual(eligibility.loc["R:20:C>T", "Unknown_Coverage"], 2)

    def test_single_library_interval_is_the_phase1_wilson_interval(self):
        row = self.profiles.loc[("R:10:A>G", "S1")]
        lower, upper = wilson_interval(np.array([40.0]), np.array([160.0]))
        self.assertAlmostEqual(row["Wilson_Lower"], lower[0])
        self.assertAlmostEqual(row["Wilson_Upper"], upper[0])
        np.testing.assert_allclose(effective_score_interval([0.25], [160]), (lower, upper))

    def test_unverified_specimen_keeps_separate_units(self):
        evidence = self.evidence.drop(columns=["Library_ID"])
        built = build_mutation_profiles(self.matrix, evidence, None, self.samples, minimum_shared=2)
        self.assertEqual(set(built["units"]["Unit_ID"]), {"S1", "S2", "L1", "L2"})
        self.assertEqual(int(built["units"].set_index("Unit_ID").loc["L1", "Same_Site_Date_Units"]), 1)


if __name__ == "__main__":
    unittest.main()
