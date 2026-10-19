import unittest

import pandas as pd

from compare_genetic_distances import compare_samples, compare_sequences


class GeneticComparisonTests(unittest.TestCase):
    def test_identical_sequences(self) -> None:
        result = compare_sequences("ACGT", "ACGT")

        self.assertEqual(result["Differences"], 0)
        self.assertEqual(result["Comparable_Positions"], 4)
        self.assertEqual(result["P_Distance"], 0.0)
        self.assertEqual(result["Percent_Identity"], 100.0)

    def test_one_mismatch(self) -> None:
        result = compare_sequences("ACGT", "ACGA")

        self.assertEqual(result["Differences"], 1)
        self.assertEqual(result["Comparable_Positions"], 4)
        self.assertEqual(result["P_Distance"], 0.25)
        self.assertEqual(result["Percent_Identity"], 75.0)

    def test_n_values_are_ignored(self) -> None:
        result = compare_sequences("ACNT", "ACGT")

        self.assertEqual(result["Differences"], 0)
        self.assertEqual(result["Comparable_Positions"], 3)
        self.assertEqual(result["P_Distance"], 0.0)
        self.assertEqual(result["Percent_Identity"], 100.0)

    def test_different_sequence_lengths_are_reported(self) -> None:
        data = pd.DataFrame(
            {
                "Sample_ID": ["S1", "S2"],
                "Pathogen": ["Pathogen A", "Pathogen A"],
                "Gene_Segment": ["Segment 1", "Segment 1"],
                "Sequence": ["ACGT", "ACGTA"],
            }
        )

        results, errors = compare_samples(data)

        self.assertTrue(results.empty)
        self.assertEqual(len(errors), 1)
        self.assertIn("S1 and S2", errors[0])
        self.assertIn("lengths are 4 and 5", errors[0])

    def test_different_pathogen_or_gene_groups_are_not_compared(self) -> None:
        data = pd.DataFrame(
            {
                "Sample_ID": ["S1", "S2", "S3"],
                "Pathogen": ["Pathogen A", "Pathogen B", "Pathogen A"],
                "Gene_Segment": ["Segment 1", "Segment 1", "Segment 2"],
                "Sequence": ["ACGT", "ACGT", "ACGT"],
            }
        )

        results, errors = compare_samples(data)

        self.assertTrue(results.empty)
        self.assertEqual(errors, [])


if __name__ == "__main__":
    unittest.main()
