"""Optional integration test that uses a real, untracked NCBI Virus download.

The ordinary suite does not depend on this dataset. It runs only when the
10-record H5N1 segment 4 files `sequences.fasta` and `sequences.csv` are placed
in the project root (they are excluded from Git by .gitignore). Otherwise it is
reported as skipped with the reason below, rather than failing.
"""

import io
import os
import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

from app import app


PROJECT = Path(__file__).resolve().parents[1]
REAL_FASTA = PROJECT / "sequences.fasta"
REAL_METADATA = PROJECT / "sequences.csv"
REAL_FILES_PRESENT = REAL_FASTA.is_file() and REAL_METADATA.is_file()


@unittest.skipUnless(
    REAL_FILES_PRESENT,
    "Optional real-data integration test: put the NCBI Virus sequences.fasta "
    "and sequences.csv in the project root to run it.",
)
class OptionalRealNcbiDataIntegrationTest(unittest.TestCase):
    def test_state_only_ncbi_upload_keeps_strict_and_candidate_results_separate(self) -> None:
        form_data = {
            "fasta_file": (io.BytesIO(REAL_FASTA.read_bytes()), "sequences.fasta"),
            "metadata_file": (io.BytesIO(REAL_METADATA.read_bytes()), "sequences.csv"),
            "max_p_distance": "0.05",
            "max_days_apart": "14",
            "max_km_apart": "50",
        }

        def fake_mafft(command, **kwargs):
            # MAFFT is mocked as a pass-through; the real files are equal length.
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
