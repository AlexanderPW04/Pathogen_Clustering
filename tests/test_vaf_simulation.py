"""SYNTHETIC checks of the Wilson interval's behaviour under its own binomial model.

These do not validate performance on real wastewater samples or lineage
detection; see validation/simulate_vaf_intervals.py for the full tables.
"""

import unittest

import numpy as np

from validation.simulate_vaf_intervals import (
    DEPTHS,
    TRUE_VAFS,
    coverage_table,
    exact_coverage,
    missing_coverage_table,
    overdispersion_table,
)


class SimulationTest(unittest.TestCase):
    def test_exact_coverage_is_near_nominal(self):
        # Exact binomial sums over every possible count (no random numbers). The Wilson
        # interval's coverage oscillates around 95%; it stays above 90% on this grid.
        values = np.array([[exact_coverage(p, n) for p in TRUE_VAFS] for n in DEPTHS])
        self.assertGreater(values.min(), 0.90)
        self.assertAlmostEqual(float(values.mean()), 0.95, delta=0.01)
        self.assertAlmostEqual(exact_coverage(0.5, 300), exact_coverage(0.5, 300))  # deterministic

    def test_simulation_is_reproducible_and_unbiased(self):
        first, second = coverage_table(replicates=2000), coverage_table(replicates=2000)
        self.assertTrue(first.equals(second))
        self.assertLess(first["Bias"].abs().max(), 0.005)  # x/n is unbiased for p under the model
        # Simulated and exact coverage agree within Monte Carlo error (sd about 0.005 at 2000 runs).
        self.assertLess((first["Simulated_Coverage"] - first["Exact_Coverage"]).abs().max(), 0.025)
        widths = first[first["True_VAF"] == 0.2].sort_values("Depth")["Mean_Width"]
        self.assertTrue(widths.is_monotonic_decreasing)

    def test_overdispersion_lowers_coverage(self):
        table = overdispersion_table(replicates=2000)
        deep = table[table["Depth"] == 300]["Coverage"]
        self.assertLess(float(deep.max()), 0.80)  # the binomial interval is too narrow for overdispersed reads

    def test_supported_non_detections_remove_slope_bias(self):
        table = missing_coverage_table(replicates=400)
        for row in table.to_dict("records"):
            truth = row["True_Slope_Per_Week"]
            self.assertAlmostEqual(row["Mean_Slope_With_Supported_Non_Detections"], truth, delta=0.0015)
        shallow = table[table["Depth"] == 50]
        # Dropping unreported low-frequency samples flattens the slope at 50x.
        self.assertTrue((shallow["Mean_Slope_Detections_Only"].abs()
                         < shallow["True_Slope_Per_Week"].abs()).all())


if __name__ == "__main__":
    unittest.main()
