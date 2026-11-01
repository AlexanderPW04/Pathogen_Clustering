"""Known-answer tests for temporal analysis, Yue-Clayton distances, and PCoA."""

import math
import unittest

import numpy as np
import pandas as pd

from genomics.diversity_analysis import (
    EXPLORATORY_ORDINATION,
    PAPER_PCOA,
    frequency_arrays,
    jaccard_distance_matrix,
    paper_distance_matrix,
    population_structure,
    principal_coordinates,
    yue_clayton_dissimilarity,
)
from genomics.pooled_observations import build_depth_lookup
from genomics.temporal_analysis import (
    depth_regression,
    frequency_trends,
    observation_matrix,
    variants_per_sample,
)
from genomics.unified_schema import (
    COVERAGE_UNKNOWN,
    COVERED_NOT_DETECTED,
    LOW_COVERAGE,
    VARIANT_OBSERVED,
)


def position_rows(sample_id, counts_by_position):
    rows = []
    for position, (a, c, g, t) in counts_by_position.items():
        rows.append({"Sample_ID": sample_id, "Reference_ID": "REF.1", "Position": position,
                     "Depth": a + c + g + t, "A_Count": a, "C_Count": c, "G_Count": g,
                     "T_Count": t})
    return rows


def base_count_table():
    rows = (
        position_rows("A", {1: (10, 0, 0, 0), 2: (10, 0, 0, 0), 3: (0, 0, 10, 0)})
        + position_rows("B", {1: (10, 0, 0, 0), 2: (5, 5, 0, 0), 3: (0, 0, 10, 0)})
        + position_rows("C", {1: (0, 10, 0, 0), 2: (0, 0, 0, 0), 3: (0, 0, 10, 0)})
    )
    table = pd.DataFrame(rows)
    table.attrs["has_base_counts"] = True
    return table


class YueClaytonTest(unittest.TestCase):
    def test_known_values(self):
        self.assertAlmostEqual(yue_clayton_dissimilarity([1, 0, 0, 0], [1, 0, 0, 0]), 0.0)
        self.assertAlmostEqual(yue_clayton_dissimilarity([1, 0, 0, 0], [0, 1, 0, 0]), 1.0)
        # shared 0.5; 0.5 + 1 - 0.5 = 1, so similarity 0.5.
        self.assertAlmostEqual(yue_clayton_dissimilarity([0.5, 0.5, 0, 0], [1, 0, 0, 0]), 0.5)
        self.assertIsNone(yue_clayton_dissimilarity([0, 0, 0, 0], [1, 0, 0, 0]))

    def test_paper_sum_sets_zero_depth_positions_to_zero(self):
        frequencies = frequency_arrays(base_count_table(), ["A", "B", "C"], "REF.1")
        self.assertTrue(np.isnan(frequencies[2, 1, 0]))  # C has no reads at position 2
        distances, compared = paper_distance_matrix(frequencies)
        # A-B: 0 + 0.5 + 0; A-C: 1 + (zero depth -> 0) + 0; B-C: 1 + 0 + 0.
        self.assertAlmostEqual(distances[0, 1], 0.5)
        self.assertAlmostEqual(distances[0, 2], 1.0)
        self.assertAlmostEqual(distances[1, 2], 1.0)
        self.assertEqual(compared[0, 2], 2)
        self.assertTrue(np.allclose(distances, distances.T))

    def test_jaccard_on_reported_sets(self):
        distances = jaccard_distance_matrix([{"x", "y"}, {"y"}, set()])
        self.assertAlmostEqual(distances[0, 1], 0.5)
        self.assertAlmostEqual(distances[0, 2], 1.0)


class PrincipalCoordinatesTest(unittest.TestCase):
    def test_euclidean_distances_are_recovered(self):
        points = np.array([[0.0, 0.0], [3.0, 0.0], [0.0, 4.0], [1.0, 1.0]])
        distances = np.sqrt(((points[:, None, :] - points[None, :, :]) ** 2).sum(axis=2))
        pcoa = principal_coordinates(distances)
        coordinates = pcoa["coordinates"]
        rebuilt = np.sqrt(((coordinates[:, None, :] - coordinates[None, :, :]) ** 2).sum(axis=2))
        self.assertTrue(np.allclose(rebuilt, distances))
        self.assertEqual(coordinates.shape[1], 2)
        self.assertAlmostEqual(float(np.sum(pcoa["explained"])), 1.0)
        self.assertEqual(pcoa["negative_eigenvalues"], 0)

    def test_non_euclidean_distances_report_negative_eigenvalues(self):
        # 1 + 1 < 3 breaks the triangle inequality, so no exact embedding exists.
        distances = np.array([[0, 1, 1], [1, 0, 3], [1, 3, 0]], dtype=float)
        pcoa = principal_coordinates(distances)
        self.assertGreaterEqual(pcoa["negative_eigenvalues"], 1)
        self.assertGreater(pcoa["negative_magnitude_ratio"], 0)


class PopulationStructureTest(unittest.TestCase):
    def samples(self, ids):
        return pd.DataFrame({
            "Sample_ID": ids, "Reference_ID": ["REF.1"] * len(ids),
            "Collection_Date": pd.date_range("2026-01-01", periods=len(ids), freq="7D"),
            "Location_Name": ["Site"] * len(ids), "Analysis_Included": [True] * len(ids),
        })

    def variants(self):
        return pd.DataFrame({
            "Sample_ID": ["A", "B", "C"], "Variant_Key": ["REF.1:2:A>C"] * 2 + ["REF.1:1:A>C"],
            "Analysis_Included": [True, True, True],
        })

    def test_paper_method_runs_only_with_base_counts(self):
        result = population_structure(self.samples(["A", "B", "C"]), self.variants(),
                                      base_count_table())
        self.assertTrue(result["available"])
        self.assertEqual(result["method"], PAPER_PCOA)
        self.assertAlmostEqual(result["distances"].loc["A", "B"], 0.5)
        self.assertIn("PCoA_1", result["coordinates"].columns)
        diversity = result["diversity"].set_index("Sample_ID")
        self.assertAlmostEqual(diversity.loc["B", "Mean_Nucleotide_Diversity"], 0.5 / 3)

        fallback = population_structure(self.samples(["A", "B", "C"]), self.variants(), None)
        self.assertEqual(fallback["method"], EXPLORATORY_ORDINATION)
        self.assertIn("cannot be computed", fallback["message"])

        too_few = population_structure(self.samples(["A", "B"]), self.variants(), base_count_table())
        self.assertFalse(too_few["available"])
        self.assertIn("at least 3", too_few["message"])


class TemporalAnalysisTest(unittest.TestCase):
    def test_unreported_variants_are_described_not_set_to_zero(self):
        variants = pd.DataFrame({
            "Variant_Key": ["REF.1:2:A>C"], "Reference_ID": ["REF.1"], "Position": [2],
            "Ref": ["A"], "Alt": ["C"], "Sample_ID": ["B"], "Alt_Frequency": [0.5],
            "Read_Depth": [10], "Analysis_Included": [True],
        })
        samples = pd.DataFrame({
            "Sample_ID": ["A", "B", "C", "D"], "Reference_ID": ["REF.1"] * 4,
            "Collection_Date": pd.to_datetime(["2026-01-01", "2026-01-08", "2026-01-15", "2026-01-22"]),
            "Location_Name": ["Site"] * 4, "Analysis_Included": [True] * 4,
        })
        table = base_count_table()
        table.loc[(table["Sample_ID"] == "C") & (table["Position"] == 2), "Depth"] = 5
        matrix = observation_matrix(variants, samples, build_depth_lookup(table), minimum_depth=10)
        status = dict(zip(matrix["Sample_ID"], matrix["Status"]))
        self.assertEqual(status["A"], COVERED_NOT_DETECTED)
        self.assertEqual(status["B"], VARIANT_OBSERVED)
        self.assertEqual(status["C"], LOW_COVERAGE)
        self.assertEqual(status["D"], COVERAGE_UNKNOWN)  # D is not in the depth table
        frequency = dict(zip(matrix["Sample_ID"], matrix["Alt_Frequency"]))
        self.assertEqual(frequency["B"], 0.5)
        self.assertTrue(math.isnan(frequency["A"]))  # not 0%
        self.assertTrue(math.isnan(frequency["D"]))

    def test_reported_frequency_trend_per_week(self):
        matrix = pd.DataFrame({
            "Variant_Key": ["K"] * 4, "Location_Name": ["Site"] * 4,
            "Collection_Date": pd.to_datetime(["2026-01-01", "2026-01-08", "2026-01-15", "2026-01-22"]),
            "Status": [VARIANT_OBSERVED] * 3 + [COVERAGE_UNKNOWN],
            "Alt_Frequency": [0.1, 0.3, 0.5, np.nan],
        })
        trend = frequency_trends(matrix).iloc[0]
        self.assertEqual(trend["Reported"], 3)
        self.assertEqual(trend["Coverage_Unknown"], 1)
        self.assertAlmostEqual(trend["Change_In_Reported_Frequency"], 0.4)
        self.assertAlmostEqual(trend["Slope_Per_Week"], 0.2)
        self.assertEqual(trend["Last_Reported_Date"], pd.Timestamp("2026-01-15"))

    def test_depth_regression_known_line_and_small_sample_flag(self):
        per_sample = pd.DataFrame({
            "Mean_Depth": [10.0, 20.0, 30.0, 40.0],
            "Distinct_SNVs_After_QC": [21, 41, 61, 81],
            "Analysis_Included": [True] * 4,
        })
        fit = depth_regression(per_sample, bootstrap_rounds=50, permutations=50)
        self.assertAlmostEqual(fit["slope"], 2.0)
        self.assertAlmostEqual(fit["intercept"], 1.0)
        self.assertAlmostEqual(fit["r_squared"], 1.0)
        self.assertTrue(fit["small_sample"])
        self.assertIsNone(depth_regression(per_sample.head(2)))

    def test_variants_per_sample_counts_records_and_distinct_snvs(self):
        variants = pd.DataFrame({
            "Sample_ID": ["A", "A", "A"], "Variant_Key": ["k1", "k1", "k2"],
            "Analysis_Included": [True, True, False],
        })
        samples = pd.DataFrame({
            "Sample_ID": ["A", "B"], "Collection_Date": pd.to_datetime(["2026-01-01"] * 2),
            "Location_Name": ["Site"] * 2, "Genome_Coverage_Percent": [95, np.nan],
            "Mean_Depth": [100, np.nan],
        })
        counts = variants_per_sample(variants, samples).set_index("Sample_ID")
        self.assertEqual(counts.loc["A", "SNV_Records"], 3)
        self.assertEqual(counts.loc["A", "Distinct_SNVs"], 2)
        self.assertEqual(counts.loc["A", "Distinct_SNVs_After_QC"], 1)
        self.assertEqual(counts.loc["B", "Distinct_SNVs"], 0)


if __name__ == "__main__":
    unittest.main()
