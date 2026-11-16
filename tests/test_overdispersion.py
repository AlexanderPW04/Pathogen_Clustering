"""Sprint 5 Phase 2: overdispersion calibration from technical replicates and its integration.

Independent checks: the moment estimator's two expectation identities are
verified by exact enumeration of beta-binomial probabilities; the plug-in
likelihood optimum is compared with a brute-force grid using the gamma-form
pmf; statistical behaviour is checked on seeded synthetic replicates.
"""

import math
import unittest

import numpy as np
import pandas as pd

from genomics import beta_binomial as bb
from genomics import overdispersion as od
from genomics.overdispersion_analysis import (
    CALIBRATED,
    DEFAULT_SENSITIVITY_RHOS,
    NO_REPLICATES,
    NO_REPLICATE_FOR_SAMPLE,
    add_overdispersion_columns,
    parse_rho_values,
    run_overdispersion_analysis,
)
from genomics.frequency_uncertainty import wilson_interval


def gamma_log_pmf(k, n, p, rho):
    alpha, beta = p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + math.lgamma(k + alpha) + math.lgamma(n - k + beta) - math.lgamma(n + alpha + beta)
            - math.lgamma(alpha) - math.lgamma(beta) + math.lgamma(alpha + beta))


def simulate_pairs(generator, sites, rho, depth=(200, 3000)):
    p = generator.uniform(0.02, 0.8, sites)
    n1, n2 = generator.integers(*depth, sites), generator.integers(*depth, sites)

    def draw(n):
        if rho == 0:
            return generator.binomial(n, p)
        return generator.binomial(n, generator.beta(p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho))
    return draw(n1), n1, draw(n2), n2


def site_table(x1, n1, x2, n2):
    return pd.DataFrame({"Alt_Reads_1": x1, "Depth_1": n1, "Alt_Reads_2": x2, "Depth_2": n2,
                         "Used": True})


class MomentIdentityTest(unittest.TestCase):
    """The estimating equation rests on two expectations; check both by exact enumeration."""

    def test_expectations_hold_exactly(self):
        for p, rho, n1, n2 in [(0.2, 0.02, 30, 45), (0.5, 0.1, 12, 20), (0.05, 0.0, 40, 25)]:
            f1 = bb.pmf_all(n1, p, rho)
            f2 = bb.pmf_all(n2, p, rho)
            k1, k2 = np.arange(n1 + 1), np.arange(n2 + 1)
            joint = np.outer(f1, f2)  # replicates are independent given the specimen mean p
            discordance = (k1[:, None] / n1 - k2[None, :] / n2) ** 2
            total = k1[:, None] + k2[None, :]
            depth = n1 + n2
            expected_discordance = float(np.sum(joint * discordance))
            expected_spread = float(np.sum(joint * total * (depth - total)))
            c = (1 + (n1 - 1) * rho) / n1 + (1 + (n2 - 1) * rho) / n2
            s = depth + rho * (n1 * (n1 - 1) + n2 * (n2 - 1))
            self.assertAlmostEqual(expected_discordance, p * (1 - p) * c, places=10)
            self.assertAlmostEqual(expected_spread, p * (1 - p) * (depth ** 2 - s), places=7)


class MomentEstimatorTest(unittest.TestCase):
    def test_identical_replicates_give_zero(self):
        x = np.array([20, 50, 100, 7])
        n = np.array([200, 500, 1000, 70])
        self.assertEqual(od.moment_estimate(x, n, x, n), 0.0)

    def test_recovers_rho_on_many_sites(self):
        generator = np.random.default_rng(1)
        for rho in [0.005, 0.02, 0.05]:
            estimate = od.moment_estimate(*simulate_pairs(generator, 3000, rho))
            self.assertAlmostEqual(estimate, rho, delta=0.15 * rho)

    def test_binomial_data_give_rho_near_zero(self):
        estimate = od.moment_estimate(*simulate_pairs(np.random.default_rng(2), 3000, 0.0))
        self.assertLess(estimate, 0.0005)

    def test_weighted_form_matches_plain_form(self):
        data = simulate_pairs(np.random.default_rng(3), 40, 0.02)
        weights = np.ones((3, 40))
        np.testing.assert_allclose(od.moment_estimate(*data, weights=weights), od.moment_estimate(*data), rtol=1e-9)

    def test_bootstrap_is_reproducible_and_brackets_the_estimate(self):
        data = simulate_pairs(np.random.default_rng(4), 40, 0.02)
        first, second = od.bootstrap_estimates(*data, replicates=300), od.bootstrap_estimates(*data, replicates=300)
        np.testing.assert_array_equal(first, second)
        estimate = od.moment_estimate(*data)
        self.assertLess(np.quantile(first, 0.025), estimate)
        self.assertGreater(np.quantile(first, 0.975), estimate)


class PluginLikelihoodTest(unittest.TestCase):
    def test_refined_optimum_matches_brute_force(self):
        x1, n1, x2, n2 = simulate_pairs(np.random.default_rng(5), 15, 0.02, depth=(150, 400))
        pooled = (x1 + x2) / (n1 + n2)

        def independent(rho):  # gamma-form pmf, not the module's rising-factorial form
            return sum(gamma_log_pmf(int(k), int(n), float(p), rho)
                       for xs, ns in ((x1, n1), (x2, n2)) for k, n, p in zip(xs, ns, pooled))
        grid = np.linspace(0.0005, 0.08, 400)
        brute = grid[int(np.argmax([independent(rho) for rho in grid]))]
        estimate, statistic = od.plugin_likelihood_estimate(x1, n1, x2, n2)
        self.assertAlmostEqual(estimate, brute, delta=0.0005)
        self.assertGreater(statistic, 0)

    def test_plugin_underestimates_compared_with_moment(self):
        # Fitting p to the same two counts halves the apparent discordance (documented bias).
        generator = np.random.default_rng(6)
        plugin, moment = [], []
        for _ in range(12):
            data = simulate_pairs(generator, 20, 0.02)
            plugin.append(od.plugin_likelihood_estimate(*data)[0])
            moment.append(od.moment_estimate(*data))
        self.assertLess(np.mean(plugin), 0.75 * 0.02)
        self.assertAlmostEqual(np.mean(moment), 0.02, delta=0.005)


class EstimateStatusTest(unittest.TestCase):
    def test_too_few_sites_gives_no_estimate(self):
        data = simulate_pairs(np.random.default_rng(7), od.MINIMUM_SITES - 1, 0.02)
        result = od.estimate_overdispersion(site_table(*data))
        self.assertEqual(result["Status"], od.TOO_FEW_SITES)
        self.assertTrue(np.isnan(result["Rho_Estimate"]) and np.isnan(result["Rho_Used_For_Intervals"]))

    def test_identical_replicates_are_estimated_zero_with_positive_upper_bound(self):
        x = np.arange(20, 220, 10)
        n = np.full(len(x), 1000)
        result = od.estimate_overdispersion(site_table(x, n, x, n))
        self.assertEqual(result["Status"], od.ESTIMATED_ZERO)
        self.assertEqual(result["Rho_Estimate"], 0.0)
        self.assertGreaterEqual(result["Rho_Used_For_Intervals"], 0.0)

    def test_wildly_discordant_replicates_are_unstable(self):
        x1 = np.array([10, 990] * 10)
        x2 = np.array([990, 10] * 10)
        n = np.full(20, 1000)
        result = od.estimate_overdispersion(site_table(x1, n, x2, n), include_plugin=False)
        self.assertEqual(result["Status"], od.UNSTABLE)
        self.assertTrue(np.isnan(result["Rho_Used_For_Intervals"]))

    def test_estimated_rho_uses_upper_bound(self):
        data = simulate_pairs(np.random.default_rng(8), 40, 0.01)
        result = od.estimate_overdispersion(site_table(*data), include_plugin=False)
        self.assertEqual(result["Status"], od.ESTIMATED)
        self.assertEqual(result["Rho_Used_For_Intervals"], result["Rho_Upper_95"])
        self.assertLess(result["Rho_Lower_95"], result["Rho_Estimate"])
        self.assertGreater(result["Share_In_Prediction_Range"], 0.85)


def evidence_rows(rows):
    """Minimal evidence records as made by build_variant_evidence (the columns calibration reads)."""

    table = pd.DataFrame(rows)
    defaults = {"Reference_ID": "R1", "Ref": "A", "Collection_Date": pd.Timestamp("2026-03-01"),
                "Location_Name": "Site", "Counted_Once": True, "Analysis_Included": True,
                "Alt_Count_Source": "Supplied", "Interval_Status": "Calculated (supplied counts)"}
    for column, value in defaults.items():
        if column not in table:
            table[column] = [value] * len(table)
        else:  # rows that did not set the column take the default
            table[column] = [value if pd.isna(cell) else cell for cell in table[column]]
    if "Library_ID" not in table:  # Phase 3: each sample documents its own library
        table["Library_ID"] = "LIB_" + table["Sample_ID"].astype(str)
    table["Variant_Key"] = table["Reference_ID"] + ":" + table["Position"].astype(str) + ":" + table["Ref"] + ">" + table["Alt"]
    return table


class ReplicatePairTest(unittest.TestCase):
    def test_pairs_need_two_samples_with_the_same_date_and_location(self):
        rows = [
            {"Sample_ID": "A1", "Specimen_ID": "A", "Position": 5, "Alt": "G"},
            {"Sample_ID": "A2", "Specimen_ID": "A", "Position": 5, "Alt": "G"},
            {"Sample_ID": "B1", "Specimen_ID": "B", "Position": 5, "Alt": "G"},
            {"Sample_ID": "B2", "Specimen_ID": "B", "Position": 5, "Alt": "G",
             "Collection_Date": pd.Timestamp("2026-03-08")},
            {"Sample_ID": "C1", "Specimen_ID": "C", "Position": 5, "Alt": "G"},
            {"Sample_ID": "D1", "Specimen_ID": "D", "Position": 5, "Alt": "G"},
            {"Sample_ID": "D2", "Specimen_ID": "D", "Position": 5, "Alt": "G", "Location_Name": "Other"},
            {"Sample_ID": "E1", "Specimen_ID": "", "Position": 5, "Alt": "G"},
        ]
        table = evidence_rows(rows)
        table.loc[3, "Collection_Date"] = pd.Timestamp("2026-03-08")
        table.loc[6, "Location_Name"] = "Other"
        pairs = {specimen: problem for specimen, _, problem in od.replicate_pairs(table)}
        self.assertEqual(pairs["A"], "")
        self.assertIn("different collection dates", pairs["B"])
        self.assertIn("exactly two", pairs["C"])
        self.assertIn("different locations", pairs["D"])
        self.assertNotIn("", pairs)


class CalibrationSiteTest(unittest.TestCase):
    def setUp(self):
        positions = []
        for sample in ["A1", "A2"]:
            for position in range(1, 8):
                positions.append({"Sample_ID": sample, "Reference_ID": "R1", "Position": position,
                                  "Depth": 500, "A_Count": 480, "C_Count": 0, "G_Count": 20, "T_Count": 0})
        self.positions = pd.DataFrame(positions)
        self.positions.attrs["has_base_counts"] = True

    def records(self):
        rows = [
            # 1: reported in both replicates
            {"Sample_ID": "A1", "Position": 1, "Alt": "G", "Alt_Read_Count": 100, "Total_Depth": 500},
            {"Sample_ID": "A2", "Position": 1, "Alt": "G", "Alt_Read_Count": 110, "Total_Depth": 500},
            # 2: reported in A1 only; A2 count comes from A/C/G/T counts (20 of 500)
            {"Sample_ID": "A1", "Position": 2, "Alt": "G", "Alt_Read_Count": 40, "Total_Depth": 500},
            # 3: reconstructed counts are not used
            {"Sample_ID": "A1", "Position": 3, "Alt": "G", "Alt_Read_Count": 60, "Total_Depth": 500,
             "Alt_Count_Source": "Reconstructed", "Interval_Status": "Approximate (reconstructed counts)"},
            {"Sample_ID": "A2", "Position": 3, "Alt": "G", "Alt_Read_Count": 61, "Total_Depth": 500},
            # 4: two alternative alleles at one position
            {"Sample_ID": "A1", "Position": 4, "Alt": "G", "Alt_Read_Count": 60, "Total_Depth": 500},
            {"Sample_ID": "A2", "Position": 4, "Alt": "T", "Alt_Read_Count": 30, "Total_Depth": 500},
            # 5: depth 100 is not above 100
            {"Sample_ID": "A1", "Position": 5, "Alt": "G", "Alt_Read_Count": 20, "Total_Depth": 100},
            {"Sample_ID": "A2", "Position": 5, "Alt": "G", "Alt_Read_Count": 22, "Total_Depth": 100},
            # 6: mean frequency above 90%
            {"Sample_ID": "A1", "Position": 6, "Alt": "G", "Alt_Read_Count": 480, "Total_Depth": 500},
            {"Sample_ID": "A2", "Position": 6, "Alt": "G", "Alt_Read_Count": 470, "Total_Depth": 500},
        ]
        table = evidence_rows(rows)
        table["Specimen_ID"] = "A"
        return table

    def test_each_site_is_used_or_explained(self):
        sites = od.calibration_sites(self.records(), self.positions, "A", ["A1", "A2"]).set_index("Position")
        self.assertTrue(sites.loc[1, "Used"])
        self.assertTrue(sites.loc[2, "Used"])
        self.assertEqual(int(sites.loc[2, "Alt_Reads_2"]), 20)
        self.assertEqual(sites.loc[2, "Count_Source_2"], "Base counts (allele not reported)")
        self.assertIn("not measured", sites.loc[3, "Exclusion_Reason"])
        self.assertTrue(sites.loc[4, "Exclusion_Reason"].str.contains("more than one alternative").all())
        self.assertIn("depth not above 100", sites.loc[5, "Exclusion_Reason"])
        self.assertIn("outside 1%", sites.loc[6, "Exclusion_Reason"])

    def test_unreported_allele_without_base_counts_is_unknown_not_zero(self):
        sites = od.calibration_sites(self.records(), None, "A", ["A1", "A2"]).set_index("Position")
        self.assertFalse(sites.loc[2, "Used"])
        self.assertTrue(pd.isna(sites.loc[2, "Alt_Reads_2"]))
        self.assertIn("count unknown", sites.loc[2, "Exclusion_Reason"])


class SensitivitySettingTest(unittest.TestCase):
    def test_parsing(self):
        self.assertEqual(parse_rho_values("0.02, 0.002"), [0.002, 0.02])
        self.assertEqual(parse_rho_values(""), [])
        self.assertEqual(parse_rho_values("0, 0.01; 0.01"), [0.01])
        self.assertEqual(parse_rho_values(None), DEFAULT_SENSITIVITY_RHOS)
        for bad in ["abc", "1", "-0.1", "0.1, 0.2, 0.3, 0.4, 0.5, 0.6"]:
            with self.assertRaises(ValueError):
                parse_rho_values(bad)


def phase1_like_evidence(samples, specimen=None, seed=9, sites=30, rho=0.01):
    """Evidence rows with supplied counts and Phase 1 Wilson intervals, for integration checks."""

    generator = np.random.default_rng(seed)
    rows = []
    p = generator.uniform(0.05, 0.6, sites)
    for sample in samples:
        n = generator.integers(300, 2000, sites)
        x = generator.binomial(n, generator.beta(p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho))
        for position, (k, depth) in enumerate(zip(x, n), start=1):
            rows.append({"Sample_ID": sample, "Position": position, "Alt": "G", "Alt_Read_Count": int(k),
                         "Total_Depth": int(depth), "VAF": k / depth})
    table = evidence_rows(rows)
    table["Specimen_ID"] = specimen or ""
    lower, upper = wilson_interval(table["Alt_Read_Count"], table["Total_Depth"])
    table["Wilson_Lower"], table["Wilson_Upper"] = lower, upper
    table["Alt_Read_Count"] = table["Alt_Read_Count"].astype("Int64")
    table["Total_Depth"] = table["Total_Depth"].astype("Int64")
    return table


class IntegrationTest(unittest.TestCase):
    def test_no_replicates_never_invents_rho(self):
        evidence = phase1_like_evidence(["S1", "S2"])
        out = run_overdispersion_analysis(evidence, None, [0.02])
        table = out["evidence"]
        self.assertTrue(table["Rho"].isna().all())
        self.assertTrue(table["BetaBinomial_Lower"].isna().all())
        self.assertEqual(set(table["Calibration_Status"]), {NO_REPLICATES})
        self.assertEqual(out["summary"]["calibrated_specimens"], 0)
        self.assertIn("could not be estimated from this dataset", out["summary"]["headline"])
        self.assertEqual(set(out["sensitivity"]["Rho_Source"]), {"Hypothetical; not calibrated"})
        # Wilson baseline columns are untouched.
        pd.testing.assert_series_equal(table["Wilson_Lower"], evidence["Wilson_Lower"])

    def test_replicated_specimen_gets_calibrated_intervals_only_for_its_samples(self):
        replicated = phase1_like_evidence(["A1", "A2"], specimen="SPEC")
        single = phase1_like_evidence(["S1"], seed=10)
        evidence = pd.concat([replicated, single], ignore_index=True)
        out = run_overdispersion_analysis(evidence, None, [])
        table = out["evidence"]
        mine = table[table["Sample_ID"].isin(["A1", "A2"])]
        self.assertTrue(mine["Calibration_Status"].eq(CALIBRATED).all())
        self.assertTrue((mine["BetaBinomial_Lower"] <= mine["VAF"] + 1e-12).all())
        self.assertTrue((mine["BetaBinomial_Upper"] >= mine["VAF"] - 1e-12).all())
        self.assertTrue((mine["BetaBinomial_Width"] >= mine["Wilson_Upper"] - mine["Wilson_Lower"] - 1e-12).all())
        other = table[table["Sample_ID"].eq("S1")]
        self.assertTrue(other["Calibration_Status"].eq(NO_REPLICATE_FOR_SAMPLE).all())
        self.assertTrue(other["Rho"].isna().all())
        self.assertEqual(out["summary"]["calibrated_specimens"], 1)
        # With no hypothetical values, only the borrowed value is offered, and it says so.
        self.assertEqual(len(out["scenarios"]), 1)
        self.assertIn("Borrowed", out["scenarios"][0][1])

    def test_samples_without_counts_get_no_interval(self):
        evidence = phase1_like_evidence(["S1"])
        evidence.loc[0, ["Wilson_Lower", "Wilson_Upper"]] = np.nan
        evidence.loc[0, "Interval_Status"] = "Not calculated: no read counts"
        table = add_overdispersion_columns(evidence, pd.DataFrame(columns=od.CALIBRATION_COLUMNS))
        self.assertTrue(table.loc[0, "Calibration_Status"].startswith("No interval"))
        self.assertEqual(table.loc[0, "Baseline_Model"], "None (no valid counts)")


if __name__ == "__main__":
    unittest.main()
