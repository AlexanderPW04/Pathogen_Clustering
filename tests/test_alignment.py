import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from align_sequences import align_samples, find_mafft
from compare_genetic_distances import compare_samples


class AlignmentTests(unittest.TestCase):
    def test_linux_mafft_executable_is_found_on_path(self):
        with patch("align_sequences.shutil.which", return_value="/usr/bin/mafft"):
            self.assertEqual(find_mafft(), "/usr/bin/mafft")

    def test_missing_mafft_gives_clear_error(self):
        with patch("align_sequences.shutil.which", return_value=None):
            with self.assertRaises(RuntimeError) as error:
                find_mafft()
        self.assertIn("Install MAFFT and add it to PATH", str(error.exception))

    def test_groups_are_aligned_separately_and_ids_are_preserved(self):
        samples = pd.DataFrame(
            {
                "Sample_ID": ["A1.1", "A2.1", "B1.1", "B2.1"],
                "Pathogen": ["Virus A", "Virus A", "Virus B", "Virus B"],
                "Gene_Segment": ["Region 1"] * 4,
                "Sequence": ["ACGT", "ACGTAA", "TTT", "TTAT"],
            }
        )
        calls = []

        def fake_mafft(command, **kwargs):
            fasta_text = Path(command[-1]).read_text(encoding="utf-8")
            calls.append(fasta_text)
            if ">A1.1" in fasta_text:
                output = ">A1.1\nACGT--\n>A2.1\nACGTAA\n"
            else:
                output = ">B1.1\nTT-T\n>B2.1\nTTAT\n"
            return subprocess.CompletedProcess(command, 0, stdout=output)

        with patch("align_sequences.subprocess.run", side_effect=fake_mafft):
            aligned = align_samples(samples, "mafft")

        self.assertEqual(len(calls), 2)
        self.assertTrue(all(not (">A1.1" in call and ">B1.1" in call) for call in calls))
        self.assertEqual(aligned["Sample_ID"].tolist(), samples["Sample_ID"].tolist())
        self.assertEqual(aligned["Sequence"].tolist(), ["ACGT--", "ACGTAA", "TT-T", "TTAT"])
        self.assertEqual(samples.loc[0, "Sequence"], "ACGT")

        pairs, errors = compare_samples(aligned)
        self.assertEqual(errors, [])
        self.assertEqual(len(pairs), 2)
        self.assertEqual(pairs.loc[0, "Comparable_Positions"], 4)
        self.assertEqual(pairs.loc[0, "Differences"], 0)

    def test_changed_mafft_ids_are_rejected(self):
        samples = pd.DataFrame(
            {
                "Sample_ID": ["A", "B"],
                "Pathogen": ["Virus", "Virus"],
                "Gene_Segment": ["Region", "Region"],
                "Sequence": ["ACGT", "ACGTA"],
            }
        )
        result = subprocess.CompletedProcess(
            ["mafft"], 0, stdout=">A\nACGT-\n>WRONG\nACGTA\n"
        )

        with patch("align_sequences.subprocess.run", return_value=result):
            with self.assertRaises(RuntimeError) as error:
                align_samples(samples, "mafft")
        self.assertIn("changed or omitted sample IDs", str(error.exception))

    def test_same_pathogen_different_segments_are_aligned_separately(self):
        samples = pd.DataFrame(
            {
                "Sample_ID": ["R1A", "R1B", "R2A", "R2B"],
                "Pathogen": ["Virus"] * 4,
                "Gene_Segment": ["Region 1", "Region 1", "Region 2", "Region 2"],
                "Sequence": ["ACGT", "ACGTA", "TTT", "TTTT"],
            }
        )
        calls = []

        def fake_mafft(command, **kwargs):
            fasta_text = Path(command[-1]).read_text(encoding="utf-8")
            calls.append(fasta_text)
            if ">R1A" in fasta_text:
                output = ">R1A\nACGT-\n>R1B\nACGTA\n"
            else:
                output = ">R2A\nTTT-\n>R2B\nTTTT\n"
            return subprocess.CompletedProcess(command, 0, stdout=output)

        with patch("align_sequences.subprocess.run", side_effect=fake_mafft):
            aligned = align_samples(samples, "mafft")

        self.assertEqual(len(calls), 2)
        self.assertTrue(all(not (">R1A" in call and ">R2A" in call) for call in calls))
        self.assertEqual(aligned["Sequence"].tolist(), ["ACGT-", "ACGTA", "TTT-", "TTTT"])

    def test_windows_batch_file_is_called_through_cmd(self):
        samples = pd.DataFrame(
            {
                "Sample_ID": ["A", "B"],
                "Pathogen": ["Virus", "Virus"],
                "Gene_Segment": ["Region", "Region"],
                "Sequence": ["ACGT", "ACGT"],
            }
        )
        result = subprocess.CompletedProcess(
            ["mafft"], 0, stdout=">A\nACGT\n>B\nACGT\n"
        )

        with patch("align_sequences.subprocess.run", return_value=result) as run:
            align_samples(samples, r"C:\MAFFT\mafft.bat")

        self.assertEqual(run.call_args.args[0][:2], ["cmd.exe", "/c"])


if __name__ == "__main__":
    unittest.main()
