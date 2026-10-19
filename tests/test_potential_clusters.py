import unittest

import pandas as pd

from detect_potential_clusters import build_clusters, find_related_pairs


class PotentialClusterTests(unittest.TestCase):
    def make_pair(
        self,
        p_distance: float = 0.01,
        days_apart: float = 5,
        km_apart: float | None = 10,
        same_location: bool = False,
    ) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "Sample_ID_1": ["A"],
                "Sample_ID_2": ["B"],
                "Pathogen": ["Pathogen A"],
                "Gene_Segment": ["Segment 1"],
                "P_Distance": [p_distance],
                "Differences": [1],
                "Days_Apart": [days_apart],
                "Km_Apart": [km_apart],
                "Same_Location": [same_location],
            }
        )

    def apply_thresholds(self, pairs: pd.DataFrame) -> pd.DataFrame:
        return find_related_pairs(pairs, 0.02, 14, 50)

    def test_pair_passes_all_three_criteria(self) -> None:
        related = self.apply_thresholds(self.make_pair())

        self.assertEqual(len(related), 1)

    def test_pair_fails_genetic_threshold(self) -> None:
        related = self.apply_thresholds(self.make_pair(p_distance=0.03))

        self.assertTrue(related.empty)

    def test_pair_fails_date_threshold(self) -> None:
        related = self.apply_thresholds(self.make_pair(days_apart=15))

        self.assertTrue(related.empty)

    def test_pair_fails_geographic_threshold(self) -> None:
        related = self.apply_thresholds(self.make_pair(km_apart=51))

        self.assertTrue(related.empty)

    def test_missing_coordinates_with_matching_location_names(self) -> None:
        pairs = self.make_pair(km_apart=None, same_location=True)
        related = self.apply_thresholds(pairs)

        self.assertEqual(len(related), 1)

    def test_transitive_pairs_form_one_three_sample_cluster(self) -> None:
        related_pairs = pd.DataFrame(
            {
                "Sample_ID_1": ["A", "B"],
                "Sample_ID_2": ["B", "C"],
                "Pathogen": ["Pathogen A", "Pathogen A"],
                "Gene_Segment": ["Segment 1", "Segment 1"],
            }
        )

        membership = build_clusters(related_pairs)

        self.assertEqual(membership["Cluster_ID"].nunique(), 1)
        self.assertEqual(set(membership["Sample_ID"]), {"A", "B", "C"})


if __name__ == "__main__":
    unittest.main()
