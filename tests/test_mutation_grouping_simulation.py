"""SYNTHETIC checks of exploratory mutation grouping (fast versions of validation/simulate_mutation_groups.py).

They pin the study's main findings, including the negative ones, so a later
change cannot silently alter them. Passing them does not show that real
pooled samples follow the simulation model, and a group is not a lineage.
"""

import unittest

import numpy as np
import pandas as pd

from validation import simulate_mutation_groups as sim

DEPTH_AWARE = ("Depth-aware (rho=0)", "depth", 0.0)
ORDINARY = ("Ordinary (h=0.10)", "ordinary", 0.0)


def mean_metric(structure, depth, rho, missing, method, metric, seeds=10, detection_only=False, assumed=None):
    values = []
    for seed in range(seeds):
        rng = np.random.default_rng([99, seed, sim.STRUCTURES.index(structure), depth, int(rho * 1000),
                                     int(missing * 10), detection_only])
        data = sim.simulate(rng, structure, depth, rho, missing, detection_only=detection_only)
        name, kind, _ = method
        row = sim.run_methods(data, 2, [(name, kind, rho if assumed is None else assumed)])[0]
        values.append(row[metric])
    return float(np.nanmean(values))


class NullCalibrationTest(unittest.TestCase):
    def test_homogeneity_test_does_not_exceed_its_level_much(self):
        table = sim.null_calibration_study(seeds=600)
        self.assertLess(float(table["Type_I_Error"].max()), 0.03)  # nominal 0.01; full study max 0.0155


class FalseGroupingTest(unittest.TestCase):
    def test_attempt1_failed_and_attempt2_fixes_look_alike_controls(self):
        attempt1 = mean_metric("lookalike", 300, 0.0, 0.3, sim.ATTEMPT1, "False_Grouping_Share")
        attempt2 = mean_metric("lookalike", 300, 0.0, 0.3, DEPTH_AWARE, "False_Grouping_Share")
        self.assertGreater(attempt1, 0.5)   # study: 0.78-0.81
        self.assertLess(attempt2, 0.10)

    def test_ordinary_clustering_groups_unrelated_mutations_far_more_often(self):
        ordinary = mean_metric("negative", 100, 0.0, 0.0, ORDINARY, "False_Grouping_Share")
        depth_aware = mean_metric("negative", 100, 0.0, 0.0, DEPTH_AWARE, "False_Grouping_Share")
        self.assertGreater(ordinary, 0.2)
        self.assertLess(depth_aware, ordinary / 2)

    def test_detection_only_profiles_mislead_ordinary_clustering(self):
        value = mean_metric("negative", 300, 0.0, 0.0, ORDINARY, "False_Grouping_Share", detection_only=True)
        self.assertGreater(value, 0.25)  # study: 0.43; detection-only input never receives groups


class RecoveryTest(unittest.TestCase):
    def test_planted_groups_are_recovered_at_high_depth_without_overdispersion(self):
        self.assertGreater(mean_metric("planted", 300, 0.0, 0.0, DEPTH_AWARE, "ARI"), 0.75)
        self.assertGreater(mean_metric("planted", 300, 0.0, 0.0, DEPTH_AWARE, "Precision"), 0.9)

    def test_overdispersion_with_the_true_rho_mostly_abstains(self):
        # A documented failure of criterion P2: correct but uninformative.
        recall = mean_metric("planted", 300, 0.02, 0.0, ("Depth-aware (rho=0.02)", "depth", 0.02), "Recall")
        self.assertLess(recall, 0.4)

    def test_low_depth_abstains_instead_of_grouping(self):
        self.assertEqual(mean_metric("planted", 20, 0.0, 0.0, DEPTH_AWARE, "Grouped_Mutations"), 0.0)


class RuleTest(unittest.TestCase):
    def test_minimum_overlap_rule_reports_failure_instead_of_choosing_silently(self):
        summary = pd.DataFrame({"Method": ["Depth-aware (rho=0)"] * 2, "Structure": ["negative"] * 2,
                                "Minimum_Shared": [2, 3], "False_Grouping_Share": [0.2, 0.15]})
        self.assertIsNone(sim.choose_minimum_shared(summary, candidates=(2, 3)))
        summary.loc[1, "False_Grouping_Share"] = 0.05
        self.assertEqual(sim.choose_minimum_shared(summary, candidates=(2, 3)), 3)


if __name__ == "__main__":
    unittest.main()
