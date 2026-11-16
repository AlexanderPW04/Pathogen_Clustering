"""SYNTHETIC checks of the Phase 2 interval procedures under their own beta-binomial model.

Fast versions of validation/simulate_overdispersion.py (the full study writes
the tables quoted in the report). Passing them does not show that real
wastewater libraries follow the model.
"""

import unittest

import numpy as np

from validation.simulate_overdispersion import (
    estimated_rho_table,
    known_rho_table,
    phase1_scenario_table,
    plugin_bias_table,
    sensitivity_table,
)


class KnownRhoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = known_rho_table()

    def test_score_interval_meets_criterion_c1_for_every_rho(self):
        for rho, group in self.table.groupby("True_Rho"):
            values = group["Beta-binomial score coverage"]
            self.assertAlmostEqual(float(values.mean()), 0.95, delta=0.01, msg=f"rho={rho}")
            self.assertGreaterEqual(float(values.min()), 0.90, msg=f"rho={rho}")

    def test_wilson_fails_when_rho_is_positive(self):
        deep = self.table[(self.table["True_Rho"] == 0.02) & (self.table["Depth"] == 300)]
        self.assertLess(float(deep["Wilson coverage"].max()), 0.90)  # 0.85 at p = 0.01 and 0.99
        interior = deep[deep["True_VAF"].between(0.02, 0.98)]
        self.assertLess(float(interior["Wilson coverage"].max()), 0.60)  # about 0.52-0.56

    def test_rho_zero_score_interval_is_the_wilson_interval(self):
        zero = self.table[self.table["True_Rho"] == 0]
        np.testing.assert_allclose(zero["Beta-binomial score coverage"], zero["Wilson coverage"])

    def test_likelihood_ratio_interval_is_not_promoted(self):
        # It fails the minimum-coverage part of C1 near the boundaries, so it is not used.
        self.assertLess(float(self.table["Beta-binomial likelihood ratio coverage"].min()), 0.90)


class Phase1ScenarioTest(unittest.TestCase):
    def test_phase1_undercoverage_is_reproduced_and_corrected(self):
        table = phase1_scenario_table(replicates=4000)
        deep = table[table["Depth"] == 300]
        self.assertTrue((deep["Wilson_Coverage"] < 0.60).all())
        gap = (table["BetaBinomial_Simulated_Coverage"] - table["BetaBinomial_Exact_Coverage"]).abs()
        standard_error = np.sqrt(table["BetaBinomial_Exact_Coverage"] * (1 - table["BetaBinomial_Exact_Coverage"]) / 4000)
        self.assertTrue((gap <= 4 * standard_error).all())


class EstimatedRhoTest(unittest.TestCase):
    def test_upper_bound_rule_meets_c3_at_the_minimum_site_count(self):
        coverage, bias = estimated_rho_table(datasets=60, sites_list=[10], rhos=[0.0, 0.02], bootstrap=300)
        upper = coverage[coverage["Rule"] == "Upper 95% bound"]
        self.assertTrue((upper["Mean_Coverage"] >= 0.93).all())
        self.assertTrue((upper["Min_Cell_Coverage"] >= 0.90).all())
        self.assertLess(abs(float(bias.loc[bias["True_Rho"] == 0.02, "Relative_Bias"].iloc[0])), 0.15)

    def test_published_plugin_estimator_is_biased_low(self):
        table = plugin_bias_table(datasets=10, sites=20).set_index("True_Rho")
        self.assertLess(table.loc[0.02, "Plugin_Relative_Bias"], -0.25)
        self.assertLess(abs(table.loc[0.02, "Moment_Relative_Bias"]), 0.25)


class SensitivityTest(unittest.TestCase):
    def test_assumed_rho_only_matches_nominal_when_it_equals_the_truth(self):
        table = sensitivity_table()
        deep = table[table["Depth"] == 300].set_index(["True_Rho", "Assumed_Rho"])["Exact_Coverage"]
        self.assertAlmostEqual(deep[(0.02, 0.02)], 0.95, delta=0.01)
        self.assertLess(deep[(0.05, 0.02)], 0.90)      # assumed too small: still too narrow
        self.assertGreater(deep[(0.0, 0.05)], 0.999)   # assumed too large: needlessly wide


if __name__ == "__main__":
    unittest.main()
