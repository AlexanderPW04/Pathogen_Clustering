import io
import os
import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from app import app, find_geography_unresolved_candidates


class FlaskInterfaceTest(unittest.TestCase):
    def test_health_check(self) -> None:
        with app.test_client() as client:
            response = client.get("/health")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.get_data(as_text=True), "OK")

    def test_upload_size_limit_has_readable_message(self) -> None:
        with patch.dict(app.config, {"MAX_CONTENT_LENGTH": 128}):
            with app.test_client() as client:
                response = client.post(
                    "/",
                    data={
                        "fasta_file": (io.BytesIO(b"A" * 256), "large.fasta"),
                        "metadata_file": (io.BytesIO(b"Accession\nA\n"), "small.csv"),
                    },
                    content_type="multipart/form-data",
                )
        self.assertEqual(response.status_code, 413)
        self.assertIn("Upload too large", response.get_data(as_text=True))

    def test_unexpected_error_does_not_show_internal_details(self) -> None:
        with patch("app.run_pipeline", side_effect=RuntimeError("private server path")):
            with app.test_client() as client:
                response = client.post(
                    "/",
                    data={
                        "fasta_file": (io.BytesIO(b">A\nACGT\n"), "a.fasta"),
                        "metadata_file": (io.BytesIO(b"Accession\nA\n"), "m.csv"),
                    },
                    content_type="multipart/form-data",
                )
        self.assertEqual(response.status_code, 500)
        self.assertIn("analysis could not be completed", response.get_data(as_text=True))
        self.assertNotIn("private server path", response.get_data(as_text=True))

    def test_dual_upload_runs_pipeline_and_displays_cluster(self) -> None:
        metadata = pd.DataFrame(
            {
                "Accession": ["A", "B"],
                "Pathogen": ["Pathogen A", "Pathogen A"],
                "Gene_Segment": ["Segment 1", "Segment 1"],
                "Collection_Date": ["2026-01-01", "2026-01-02"],
                "Location_Name": ["Location A", "Location A"],
                "Latitude": [33.45, 33.45],
                "Longitude": [-112.07, -112.07],
                "Host_Source": ["Human", "Human"],
            }
        )

        form_data = {
            "fasta_file": (io.BytesIO(b">A\nACGT\n>B\nACGA\n"), "sequences.fasta"),
            "metadata_file": (
                io.BytesIO(metadata.to_csv(index=False).encode("utf-8")),
                "metadata.csv",
            ),
            "max_p_distance": "0.25",
            "max_days_apart": "10",
            "max_km_apart": "1",
        }

        def fake_mafft(command, **kwargs):
            fasta_text = Path(command[-1]).read_text(encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout=fasta_text)

        with patch("align_sequences.shutil.which", return_value="mafft"):
            with patch("align_sequences.subprocess.run", side_effect=fake_mafft):
                with app.test_client() as client:
                    response = client.post(
                        "/",
                        data=form_data,
                        content_type="multipart/form-data",
                    )

        page = response.get_data(as_text=True)
        self.assertEqual(response.status_code, 200)
        self.assertIn("Potential clusters:</strong> 1", page)
        self.assertIn("<h3>C1</h3>", page)
        self.assertIn("A, B", page)
        self.assertIn("A: Segment 1 (Metadata)", page)
        self.assertIn("sequences.fasta", page)
        self.assertIn("metadata.csv", page)
        self.assertIn("Sequences aligned by MAFFT:</strong> 2", page)
        self.assertIn("Unique samples analyzed:</strong> 2", page)
        self.assertIn("Pairwise comparisons:</strong> 1", page)
        self.assertIn("Strict related pairs:</strong> 1", page)
        self.assertNotIn("Synthetic example analysis", page)
        self.assertIn("MAFFT aligns sequences", page)
        self.assertIn("GeoNames may estimate them", page)

    def test_missing_geonames_username_is_shown_on_page(self) -> None:
        metadata = pd.DataFrame(
            {
                "Accession": ["A", "B"],
                "Pathogen": ["Virus", "Virus"],
                "Collection_Date": ["2026-01-01", "2026-01-02"],
                "Location_Name": ["Phoenix", "Phoenix"],
                "Host_Source": ["Human", "Human"],
            }
        )
        form_data = {
            "fasta_file": (io.BytesIO(b">A\nACGT\n>B\nACGT\n"), "sequences.fasta"),
            "metadata_file": (
                io.BytesIO(metadata.to_csv(index=False).encode("utf-8")),
                "metadata.csv",
            ),
            "max_p_distance": "0.05",
            "max_days_apart": "7",
            "max_km_apart": "10",
            "genomic_region_mode": "complete_genomes",
            "confirm_complete_genomes": "on",
        }

        def fake_mafft(command, **kwargs):
            fasta_text = Path(command[-1]).read_text(encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout=fasta_text)

        with patch.dict(os.environ, {"GEONAMES_USERNAME": ""}):
            with patch("align_sequences.shutil.which", return_value="mafft"):
                with patch("align_sequences.subprocess.run", side_effect=fake_mafft):
                    with app.test_client() as client:
                        response = client.post(
                            "/", data=form_data, content_type="multipart/form-data"
                        )

        page = response.get_data(as_text=True)
        self.assertIn("GeoNames lookup skipped", page)
        self.assertIn("Potential clusters:</strong> 1", page)
        self.assertIn("Complete genome (User-confirmed complete genome)", page)

    def test_complete_genomes_needs_explicit_confirmation(self) -> None:
        with app.test_client() as client:
            response = client.post(
                "/",
                data={"genomic_region_mode": "complete_genomes"},
            )
        self.assertIn(
            "Please confirm that all sequences are complete genomes",
            response.get_data(as_text=True),
        )

    def test_automatic_detection_reports_unknown_region(self) -> None:
        metadata = pd.DataFrame(
            {
                "Accession": ["A", "B"],
                "Pathogen": ["Virus", "Virus"],
                "Collection_Date": ["2026-01-01", "2026-01-02"],
                "Location_Name": ["Tempe", "Tempe"],
                "Host_Source": ["Synthetic", "Synthetic"],
            }
        )
        with app.test_client() as client:
            response = client.post(
                "/",
                data={
                    "fasta_file": (io.BytesIO(b">A\nACGT\n>B\nACGT\n"), "a.fasta"),
                    "metadata_file": (
                        io.BytesIO(metadata.to_csv(index=False).encode("utf-8")),
                        "m.csv",
                    ),
                },
                content_type="multipart/form-data",
            )
        page = response.get_data(as_text=True)
        self.assertIn("Cannot determine genomic region", page)
        self.assertIn("Select Complete genomes only if", page)

    def test_invalid_upload_and_missing_mafft_are_visible(self) -> None:
        with app.test_client() as client:
            invalid = client.post(
                "/",
                data={
                    "fasta_file": (io.BytesIO(b"not fasta"), "wrong.txt"),
                    "metadata_file": (io.BytesIO(b"Accession\nA\n"), "data.csv"),
                },
                content_type="multipart/form-data",
            )
        self.assertIn("must use .fasta, .fa, or .fna", invalid.get_data(as_text=True))

        metadata = pd.DataFrame(
            {
                "Accession": ["A", "B"],
                "Pathogen": ["Virus", "Virus"],
                "Gene_Segment": ["Segment 1", "Segment 1"],
                "Collection_Date": ["2026-01-01", "2026-01-02"],
                "Location_Name": ["Tempe", "Tempe"],
                "Host_Source": ["Synthetic", "Synthetic"],
            }
        )
        with patch("align_sequences.shutil.which", return_value=None):
            with app.test_client() as client:
                missing_mafft = client.post(
                    "/",
                    data={
                        "fasta_file": (io.BytesIO(b">A\nACGT\n>B\nACGT\n"), "a.fasta"),
                        "metadata_file": (
                            io.BytesIO(metadata.to_csv(index=False).encode("utf-8")),
                            "m.csv",
                        ),
                    },
                    content_type="multipart/form-data",
                )
        self.assertIn(
            "Install MAFFT and add it to PATH", missing_mafft.get_data(as_text=True)
        )

    def test_state_only_pair_is_a_candidate_not_a_strict_cluster(self) -> None:
        metadata = pd.DataFrame(
            {
                "Accession": ["A", "B"],
                "Pathogen": ["Example pathogen"] * 2,
                "Gene_Segment": ["Segment 1"] * 2,
                "Collection_Date": ["2026-01-01", "2026-01-02"],
                "Location_Name": ["USA: ND", "USA: ND"],
                "Host_Source": ["Synthetic"] * 2,
            }
        )
        form_data = {
            "fasta_file": (io.BytesIO(b">A\nACGT\n>B\nACGT\n"), "states.fasta"),
            "metadata_file": (
                io.BytesIO(metadata.to_csv(index=False).encode("utf-8")),
                "states.csv",
            ),
            "max_p_distance": "0.05",
            "max_days_apart": "14",
            "max_km_apart": "50",
        }

        def fake_mafft(command, **kwargs):
            fasta_text = Path(command[-1]).read_text(encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout=fasta_text)

        with patch.dict(os.environ, {"GEONAMES_USERNAME": ""}):
            with patch("align_sequences.shutil.which", return_value="mafft"):
                with patch("align_sequences.subprocess.run", side_effect=fake_mafft):
                    with app.test_client() as client:
                        response = client.post(
                            "/", data=form_data, content_type="multipart/form-data"
                        )

        page = response.get_data(as_text=True)
        self.assertIn("Potential clusters:</strong> 0", page)
        self.assertIn("candidate pairs with unresolved geography:</strong> 1", page)
        self.assertIn("Distinct locations too coarse for distance comparison:</strong> 1", page)
        self.assertIn("Samples lacking sufficiently precise geographic metadata:</strong> 2", page)
        self.assertIn("Unique samples represented in candidate pairs:</strong> 2", page)
        self.assertIn("Show candidate pairs (1)", page)
        self.assertIn('<details class="candidate-details">', page)
        self.assertIn("A – B", page)

    def test_known_distance_over_cutoff_is_not_unresolved(self) -> None:
        pairs = pd.DataFrame(
            {
                "Sample_ID_1": ["A", "C"],
                "Sample_ID_2": ["B", "D"],
                "Pathogen": ["Virus", "Virus"],
                "Gene_Segment": ["Segment 1", "Segment 1"],
                "P_Distance": [0.0, 0.0],
                "Days_Apart": [1, 1],
                "Km_Apart": [100.0, None],
                "Same_Location": [False, None],
            }
        )
        candidates = find_geography_unresolved_candidates(pairs, [0.05, 14, 50])
        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0]["sample_1"], "C")

    def test_example_button_runs_synthetic_dataset(self) -> None:
        def fake_mafft(command, **kwargs):
            lines = Path(command[-1]).read_text(encoding="utf-8").splitlines()
            for index, line in enumerate(lines):
                if line == ">SYN_C":
                    lines[index + 1] += "-"
            return subprocess.CompletedProcess(
                command, 0, stdout="\n".join(lines) + "\n"
            )

        with patch.dict(os.environ, {"GEONAMES_USERNAME": ""}):
            with patch("align_sequences.shutil.which", return_value="mafft"):
                with patch("align_sequences.subprocess.run", side_effect=fake_mafft):
                    with app.test_client() as client:
                        response = client.post(
                            "/",
                            data={
                                "run_example": "1",
                                "max_p_distance": "0.05",
                                "max_days_apart": "14",
                                "max_km_apart": "50",
                            },
                        )

        page = response.get_data(as_text=True)
        self.assertIn("synthetic_sequences.fasta", page)
        self.assertIn("synthetic_metadata.csv", page)
        self.assertIn("Sequences aligned by MAFFT:</strong> 8", page)
        self.assertIn("Pairwise comparisons:</strong> 28", page)
        self.assertIn("Potential clusters:</strong> 2", page)
        self.assertIn("Synthetic example analysis", page)

    def test_state_only_ncbi_upload_keeps_strict_and_candidate_results_separate(self) -> None:
        project = Path(__file__).resolve().parents[1]
        form_data = {
            "fasta_file": (
                io.BytesIO((project / "sequences.fasta").read_bytes()),
                "sequences.fasta",
            ),
            "metadata_file": (
                io.BytesIO((project / "sequences.csv").read_bytes()),
                "sequences.csv",
            ),
            "max_p_distance": "0.05",
            "max_days_apart": "14",
            "max_km_apart": "50",
        }

        def fake_mafft(command, **kwargs):
            fasta_text = Path(command[-1]).read_text(encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout=fasta_text)

        with patch.dict(os.environ, {"GEONAMES_USERNAME": ""}):
            with patch("align_sequences.shutil.which", return_value="mafft"):
                with patch("align_sequences.subprocess.run", side_effect=fake_mafft):
                    with app.test_client() as client:
                        response = client.post(
                            "/", data=form_data, content_type="multipart/form-data"
                        )

        page = response.get_data(as_text=True)
        self.assertIn("Potential clusters:</strong> 0", page)
        self.assertIn("their geographic evidence is insufficient", page)
        self.assertIn("Unique samples analyzed:</strong> 10", page)
        self.assertIn("Pairwise comparisons:</strong> 45", page)
        self.assertIn("Strict related pairs:</strong> 0", page)
        self.assertIn("Genomic-temporal candidate pairs with unresolved geography:</strong> 22", page)
        self.assertIn("Unique samples represented in candidate pairs:</strong> 9", page)
        self.assertIn("Candidate pairs are relationships among 9 of the 10 original samples", page)
        self.assertIn("Show candidate pairs (22)", page)
        self.assertIn('<details class="candidate-details">', page)
        self.assertIn("P-distance: 0.001174", page)
        self.assertNotIn("Synthetic example analysis", page)
        self.assertIn("Distinct locations too coarse for distance comparison:</strong>", page)
        self.assertNotIn(
            "candidate pairs with unresolved geography:</strong> 0", page
        )


if __name__ == "__main__":
    unittest.main()
