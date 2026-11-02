"""Flask tests for the unified dashboard (mocked MAFFT, GeoNames disabled)."""

import io
import os
import re
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402

from app import app  # noqa: E402

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "integrated"
REFERENCE_TEXT = (EXAMPLE / "synthetic_reference.fasta").read_text(encoding="utf-8")
REFERENCE = REFERENCE_TEXT.split("\n", 1)[1].replace("\n", "")


def upload(name):
    return (io.BytesIO((EXAMPLE / name).read_bytes()), name)


def post(data, multipart=False):
    with patch.dict(os.environ, {"GEONAMES_USERNAME": ""}):
        with patch("align_sequences.shutil.which", return_value="mafft"):
            with patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(REFERENCE)):
                with app.test_client() as client:
                    if multipart:
                        response = client.post("/", data=data, content_type="multipart/form-data")
                    else:
                        response = client.post("/", data=data)
                    page = response.get_data(as_text=True)
                    links = re.findall(r'href="(/results/[^"]+)"', page)
                    downloads = {link: client.get(link) for link in links[:3]}
    return response, page, links, downloads


class DashboardTest(unittest.TestCase):
    def test_integrated_example_shows_one_dashboard_with_both_sources(self):
        response, page, links, downloads = post({"run_integrated_example": "1"})
        self.assertEqual(response.status_code, 200)
        self.assertIn("SYNTHETIC integrated example", page)
        self.assertNotIn("Synthetic example analysis", page)
        # Integrated (both-source) sections come before source-specific ones.
        headings = ["1. Overview", "2. Cross-source evidence", "3. Genomic variation",
                    "4. Temporal patterns", "5. Geographic context", "6. Sequencing quality",
                    "7. Population structure", "8. Individual sequence relationships"]
        positions = [page.index(f"<h2>{heading}</h2>") for heading in headings]
        self.assertEqual(positions, sorted(positions))
        self.assertIn("Shared: observed in both sources", page)
        self.assertIn("Incomparable: reference not verified as the same", page)
        self.assertIn("Variants observed in both sources:</strong> 2", page)
        self.assertIn("Paper method: Yue-Clayton PCoA", page)
        self.assertIn("Which Fontenele et al. (2021) methods were reproduced", page)
        self.assertIn("Profile overlap only; not a lineage assignment", page)
        self.assertIn('src="data:image/png;base64,', page)
        self.assertIn("Potential clusters:</strong>", page)  # Sprint 2 results are kept
        self.assertTrue(any(link.endswith("cross_source_comparison.csv") for link in links))
        for link, download in downloads.items():
            self.assertEqual(download.status_code, 200, link)

    def test_result_download_returns_csv_and_unknown_ids_404(self):
        _, page, links, _ = post({"run_integrated_example": "1"})
        link = next(link for link in links if link.endswith("cross_source_comparison.csv"))
        with app.test_client() as client:
            csv_response = client.get(link)
            missing = client.get("/results/not-a-real-id/cross_source_comparison.csv")
        self.assertEqual(csv_response.status_code, 200)
        self.assertIn("Variant_Key", csv_response.get_data(as_text=True))
        self.assertEqual(missing.status_code, 404)

    def test_pooled_only_upload_runs_without_individual_sequences(self):
        _, page, _, _ = post({"snv_file": upload("synthetic_pooled_snvs.csv"), "qc_mode": "paper"},
                             multipart=True)
        self.assertNotIn("Analysis could not be completed", page)
        self.assertIn("Pooled samples:</strong> 9", page)
        self.assertIn("EXPLORATORY", page.upper())  # no base counts: not the paper's PCoA
        self.assertIn("cannot be computed", page)
        self.assertIn("Cross-source comparison needs individual sequences", page)
        self.assertIn("Upload individual FASTA sequences and metadata", page)
        self.assertNotIn("Potential clusters:</strong>", page)

    def test_individual_sequences_with_reference_get_variant_profiles(self):
        _, page, links, _ = post({
            "fasta_file": upload("synthetic_individual_sequences.fasta"),
            "metadata_file": upload("synthetic_individual_metadata.csv"),
            "reference_file": upload("synthetic_reference.fasta"),
        }, multipart=True)
        self.assertIn("Individual sequences:</strong> 5 <span", page)
        self.assertIn("(4 mapped to the reference)", page)
        self.assertIn("Reference mapping not established", page)
        self.assertIn("Cross-source comparison needs individual sequences", page)
        self.assertTrue(any(link.endswith("individual_variants.csv") for link in links))

    def test_pooled_side_files_need_an_snv_table_and_bad_reference_is_reported(self):
        _, page, _, _ = post({"position_file": upload("synthetic_pooled_base_counts.csv")},
                             multipart=True)
        self.assertIn("Please choose a FASTA sequence file with its metadata", page)

        _, page, _, _ = post({
            "snv_file": upload("synthetic_pooled_snvs.csv"),
            "position_file": upload("synthetic_pooled_base_counts.csv"),
        }, multipart=True)
        self.assertIn("Paper method: Yue-Clayton PCoA", page)

        two_records = (io.BytesIO(b">A\nACGT\n>B\nACGT\n"), "reference.fasta")
        _, page, _, _ = post({"snv_file": upload("synthetic_pooled_snvs.csv"),
                              "reference_file": two_records}, multipart=True)
        self.assertIn("exactly one sequence", page)

    def test_invalid_date_window_and_qc_mode_are_rejected(self):
        _, page, _, _ = post({"snv_file": upload("synthetic_pooled_snvs.csv"),
                              "date_window_days": "-3"}, multipart=True)
        self.assertIn("date window must be a non-negative number", page)
        _, page, _, _ = post({"snv_file": upload("synthetic_pooled_snvs.csv"),
                              "qc_mode": "strict"}, multipart=True)
        self.assertIn("paper-reproduction or exploratory", page)


if __name__ == "__main__":
    unittest.main()
