import contextlib
import io
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from run_pipeline import main, run_pipeline


class PipelineIntegrationTest(unittest.TestCase):
    def test_command_line_passes_region_options_to_pipeline(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            fasta_path = Path(folder) / "sequences.fasta"
            metadata_path = Path(folder) / "metadata.csv"
            fasta_path.write_text(">A\nACGT\n", encoding="utf-8")
            metadata_path.write_text("Accession\nA\n", encoding="utf-8")
            arguments = [
                "run_pipeline.py", str(fasta_path), str(metadata_path),
                "--max-p-distance", "0.05", "--max-days-apart", "14",
                "--max-km-apart", "50", "--default-gene-segment", "Segment 4",
                "--confirm-same-region",
            ]
            with patch("sys.argv", arguments):
                with patch("run_pipeline.run_pipeline", return_value=0) as runner:
                    status = main()

        self.assertEqual(status, 0)
        runner.assert_called_once_with(
            fasta_path, metadata_path, 0.05, 14.0, 50.0,
            default_gene_segment="Segment 4", confirm_same_region=True,
        )

    def test_command_line_accepts_complete_genome_confirmation(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            fasta_path = Path(folder) / "sequences.fasta"
            metadata_path = Path(folder) / "metadata.csv"
            fasta_path.write_text(">A\nACGT\n", encoding="utf-8")
            metadata_path.write_text("Accession\nA\n", encoding="utf-8")
            arguments = [
                "run_pipeline.py", str(fasta_path), str(metadata_path),
                "--max-p-distance", "0.05", "--max-days-apart", "14",
                "--max-km-apart", "50", "--complete-genomes",
            ]
            with patch("sys.argv", arguments):
                with patch("run_pipeline.run_pipeline", return_value=0) as runner:
                    status = main()

        self.assertEqual(status, 0)
        self.assertTrue(runner.call_args.kwargs["complete_genomes"])

    def test_complete_pipeline_creates_final_outputs(self) -> None:
        data = pd.DataFrame(
            {
                "Accession": ["A", "B", "C"],
                "Organism_Name": ["Pathogen A"] * 3,
                "Segment": ["Segment 1"] * 3,
                "Collection_Date": ["2026-01-01", "2026-01-03", "2026-01-05"],
                "Geo_Location": ["Location A"] * 3,
                "Host": ["Human"] * 3,
            }
        )

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = Path(folder) / "sequences.fasta"
            fasta_path.write_text(
                ">A.1\nACGT\n>B.1\nACGA\n>C.1\nACGA\n",
                encoding="utf-8",
            )
            metadata_path = Path(folder) / "samples.xlsx"
            data.to_excel(metadata_path, index=False, engine="openpyxl")

            def fake_mafft(command, **kwargs):
                fasta_text = Path(command[-1]).read_text(encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, stdout=fasta_text)

            with patch.dict(os.environ, {"GEONAMES_USERNAME": ""}):
                with patch("align_sequences.shutil.which", return_value="mafft"):
                    with patch("align_sequences.subprocess.run", side_effect=fake_mafft):
                        with contextlib.redirect_stdout(io.StringIO()) as output:
                            status = run_pipeline(fasta_path, metadata_path, 0.25, 10, 1)

            related_path = Path(folder) / (
                "samples_pairwise_comparisons_temporal_geographic_related_pairs.csv"
            )
            clusters_path = Path(folder) / (
                "samples_pairwise_comparisons_temporal_geographic_"
                "cluster_membership.csv"
            )

            self.assertEqual(status, 0)
            self.assertTrue(related_path.is_file())
            self.assertTrue(clusters_path.is_file())

            membership = pd.read_csv(clusters_path)
            self.assertEqual(
                set(membership["Sample_ID"]), {"A.1", "B.1", "C.1"}
            )
            self.assertEqual(membership["Cluster_ID"].nunique(), 1)
            self.assertIn("1/6 Validating and matching FASTA", output.getvalue())
            self.assertIn("3/6 Aligning each pathogen/gene group", output.getvalue())
            self.assertIn("6/6 Detecting potential clusters", output.getvalue())

    def test_missing_mafft_stops_before_creating_results(self) -> None:
        metadata = pd.DataFrame(
            {
                "Accession": ["A", "B"],
                "Pathogen": ["Virus", "Virus"],
                "Gene_Segment": ["Region", "Region"],
                "Collection_Date": ["2026-01-01", "2026-01-02"],
                "Location_Name": ["Phoenix", "Phoenix"],
                "Host_Source": ["Human", "Human"],
            }
        )

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = Path(folder) / "sequences.fasta"
            fasta_path.write_text(">A\nACGT\n>B\nACGT\n", encoding="utf-8")
            metadata_path = Path(folder) / "metadata.csv"
            metadata.to_csv(metadata_path, index=False)

            with patch("align_sequences.shutil.which", return_value=None):
                with contextlib.redirect_stderr(io.StringIO()) as error_output:
                    status = run_pipeline(fasta_path, metadata_path, 0.1, 7, 20)

            self.assertEqual(status, 2)
            self.assertIn("Install MAFFT and add it to PATH", error_output.getvalue())
            self.assertFalse((Path(folder) / "metadata_pairwise_comparisons.csv").exists())

    def test_different_fasta_segments_are_aligned_separately(self) -> None:
        metadata = pd.DataFrame(
            {
                "Accession": ["A", "B", "C", "D"],
                "Pathogen": ["Virus"] * 4,
                "Collection_Date": ["2026-01-01"] * 4,
                "Location_Name": ["Phoenix"] * 4,
                "Latitude": [33.4] * 4,
                "Longitude": [-112.0] * 4,
                "Host_Source": ["Human"] * 4,
            }
        )
        mafft_inputs = []

        def fake_mafft(command, **kwargs):
            fasta_text = Path(command[-1]).read_text(encoding="utf-8")
            mafft_inputs.append(fasta_text)
            if ">A" in fasta_text:
                output = ">A\nACGT-\n>B\nACGTA\n"
            else:
                output = ">C\nTTT-\n>D\nTTTA\n"
            return subprocess.CompletedProcess(command, 0, stdout=output)

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = Path(folder) / "sequences.fasta"
            fasta_path.write_text(
                ">A segment 4\nACGT\n>B segment 4\nACGTA\n"
                ">C segment 5\nTTT\n>D segment 5\nTTTA\n",
                encoding="utf-8",
            )
            metadata_path = Path(folder) / "metadata.csv"
            metadata.to_csv(metadata_path, index=False)

            with patch("align_sequences.shutil.which", return_value="mafft"):
                with patch("align_sequences.subprocess.run", side_effect=fake_mafft):
                    with contextlib.redirect_stdout(io.StringIO()):
                        status = run_pipeline(fasta_path, metadata_path, 1, 7, 10)

            pairs = pd.read_csv(Path(folder) / "metadata_pairwise_comparisons.csv")
            assignments = pd.read_csv(Path(folder) / "metadata_region_assignments.csv")

        self.assertEqual(status, 0)
        self.assertEqual(len(mafft_inputs), 2)
        self.assertTrue(all(not (">A" in text and ">C" in text) for text in mafft_inputs))
        self.assertEqual(set(pairs["Gene_Segment"]), {"Segment 4", "Segment 5"})
        self.assertEqual(len(pairs), 2)
        self.assertTrue((assignments["Gene_Segment_Source"] == "FASTA header").all())

    def test_both_existing_examples_still_produce_clusters(self) -> None:
        example_folder = (
            Path(__file__).resolve().parents[1]
            / "outputs"
            / "01a0e5b9-470b-7032-a8c4-c8465c89040b"
        )
        examples = [
            ("example_sequences.fasta", "example_metadata.csv", 36),
            ("ncbi_virus_style_sequences.fasta", "ncbi_virus_style_metadata.csv", 45),
        ]

        def fake_mafft(command, **kwargs):
            fasta_text = Path(command[-1]).read_text(encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout=fasta_text)

        for fasta_name, metadata_name, expected_pairs in examples:
            with self.subTest(metadata=metadata_name):
                with tempfile.TemporaryDirectory() as folder:
                    fasta_path = Path(folder) / fasta_name
                    metadata_path = Path(folder) / metadata_name
                    shutil.copyfile(example_folder / fasta_name, fasta_path)
                    shutil.copyfile(example_folder / metadata_name, metadata_path)

                    with patch.dict(os.environ, {"GEONAMES_USERNAME": ""}):
                        with patch("align_sequences.shutil.which", return_value="mafft"):
                            with patch(
                                "align_sequences.subprocess.run",
                                side_effect=fake_mafft,
                            ):
                                with contextlib.redirect_stdout(io.StringIO()):
                                    status = run_pipeline(
                                        fasta_path, metadata_path, 0.05, 7, 10
                                    )

                    stem = metadata_path.stem
                    genomic = Path(folder) / f"{stem}_pairwise_comparisons.csv"
                    related = Path(folder) / (
                        f"{stem}_pairwise_comparisons_temporal_geographic_"
                        "related_pairs.csv"
                    )
                    membership = Path(folder) / (
                        f"{stem}_pairwise_comparisons_temporal_geographic_"
                        "cluster_membership.csv"
                    )

                    self.assertEqual(status, 0)
                    self.assertEqual(len(pd.read_csv(genomic)), expected_pairs)
                    self.assertEqual(len(pd.read_csv(related)), 4)
                    members = pd.read_csv(membership)
                    self.assertEqual(members["Cluster_ID"].nunique(), 2)
                    self.assertEqual(len(members), 5)


if __name__ == "__main__":
    unittest.main()
