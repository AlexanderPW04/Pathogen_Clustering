"""Mode-specific upload fields: declared roles, conflicts, and each analysis mode.

MAFFT is mocked with tests/fake_mafft.py and GeoNames is disabled, so these
tests check the upload and review behavior, not the external tools. Mode
switching in the browser (JavaScript) was checked separately in a browser.
"""

import io
import os
import re
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402

from app import STAGED, UPLOAD_FIELDS, app  # noqa: E402
from genomics.input_inspection import apply_declared_role, empty_info  # noqa: E402

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "integrated"
REFERENCE = (EXAMPLE / "synthetic_reference.fasta").read_text(encoding="utf-8").split("\n", 1)[1].replace("\n", "")
_SAMPLES = pd.read_csv(EXAMPLE / "synthetic_pooled_snvs.csv").drop_duplicates("Sample_ID")
# Pooled metadata must list every pooled sample (an existing review rule).
POOLED_METADATA = _SAMPLES.assign(Genome_Coverage_Percent=95)[
    ["Sample_ID", "Collection_Date", "Location_Name", "Genome_Coverage_Percent"]].to_csv(index=False).encode()


def example(name):
    return (io.BytesIO((EXAMPLE / name).read_bytes()), name)


def token_of(page):
    match = re.search(r'name="token" value="([^"]+)"', page)
    return match.group(1) if match else None


class ModeUploadTest(unittest.TestCase):
    def setUp(self):
        patches = [patch.dict(os.environ, {"GEONAMES_USERNAME": ""}),
                   patch("align_sequences.shutil.which", return_value="mafft"),
                   patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(REFERENCE))]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        self.addCleanup(STAGED.clear)
        self.client = app.test_client()

    def review(self, mode, **fields):
        data = {"intended_mode": mode}
        for key, value in fields.items():
            data[key] = value if isinstance(value, list) else [value]
        page = self.client.post("/review", data=data, content_type="multipart/form-data").get_data(as_text=True)
        return page, token_of(page)

    def run_reviewed(self, token):
        with redirect_stdout(io.StringIO()):
            return self.client.post("/run", data={"token": token}).get_data(as_text=True)

    # Page --------------------------------------------------------------------------
    def test_home_shows_the_automatic_box_and_hides_mode_fields(self):
        page = self.client.get("/").get_data(as_text=True)
        self.assertIn('data-modes="auto" >', page)
        self.assertIn('id="files" type="file" name="files" multiple required', page)
        for field in UPLOAD_FIELDS:
            self.assertIn(f'name="file_{field["role"]}"', page)
        # Mode fields start hidden and are not required while Automatic detection is chosen.
        fasta_block = page[page.index('data-modes="individual combined"'):]
        self.assertIn("hidden", fasta_block[:fasta_block.index(">")])
        self.assertNotIn('name="file_fasta" accept=".fasta,.fa,.fna,.gz" multiple required', page)
        self.assertIn("<summary>Advanced inputs (optional)</summary>", page)
        self.assertIn("it does not accept raw FASTQ reads or call variants", page)
        self.assertIn("same verified reference", page)
        self.assertIn("Kept but not uploaded in this mode", page)  # the switching script

    def test_required_and_optional_status_per_mode(self):
        status = {field["role"]: field["modes"] for field in UPLOAD_FIELDS}
        self.assertEqual(status["fasta"], {"individual": "required", "combined": "required"})
        self.assertEqual(status["metadata"], {"individual": "required", "combined": "required"})
        self.assertEqual(status["reference"], {"individual": "optional", "pooled": "optional", "combined": "required"})
        self.assertEqual(status["snv"], {"pooled": "required", "combined": "required"})
        for role in ["pooled_metadata", "positions", "masked_sites", "signatures"]:
            self.assertEqual(status[role], {"pooled": "optional", "combined": "optional"})
        advanced = {field["role"] for field in UPLOAD_FIELDS if field.get("advanced")}
        self.assertEqual(advanced, {"masked_sites", "signatures"})

    def test_a_chosen_mode_is_shown_again_after_going_back(self):
        page = self.client.post("/back", data={"target": "upload", "intended_mode": "pooled"}).get_data(as_text=True)
        self.assertIn('value="pooled" checked', page)
        self.assertIn('name="file_snv" accept=".csv,.tsv,.txt,.xlsx,.vcf,.gz" multiple required', page)
        self.assertNotIn('id="files" type="file" name="files" multiple required', page)

    # No files ----------------------------------------------------------------------
    def test_no_files_in_any_mode_is_a_clear_message(self):
        for mode in ["auto", "individual", "pooled", "combined"]:
            page, token = self.review(mode)
            self.assertIn("Choose one or more files to upload", page)
            self.assertIsNone(token)

    # Individual (Sprint 2) -----------------------------------------------------------
    def test_individual_required_files_only(self):
        page, token = self.review("individual", file_fasta=example("synthetic_individual_sequences.fasta"),
                                  file_metadata=example("synthetic_individual_metadata.csv"))
        self.assertEqual(page.count("Your choice"), 2)
        self.assertIn("uploaded as: Individual sequences (FASTA)", page)
        self.assertIn("Detected analysis:</strong> Individual sequences (Sprint 2)", page)
        self.assertIn("the files support it", page)
        self.assertIn("5 sample IDs", page)
        self.assertNotIn("disabled>Run analysis", page)
        result = self.run_reviewed(token)
        self.assertIn("Pairwise comparisons:</strong> 6", result)
        self.assertIn("Individual sequence analysis (Sprint 2)", result)

    def test_individual_with_the_optional_reference(self):
        page, token = self.review("individual", file_fasta=example("synthetic_individual_sequences.fasta"),
                                  file_metadata=example("synthetic_individual_metadata.csv"),
                                  file_reference=example("synthetic_reference.fasta"))
        self.assertEqual(page.count("Your choice"), 3)
        result = self.run_reviewed(token)
        self.assertIn("(4 mapped to the reference)", result)

    # Pooled (Sprint 3) ---------------------------------------------------------------
    def test_pooled_required_only_and_with_optional_files(self):
        page, token = self.review("pooled", file_snv=example("synthetic_pooled_snvs.csv"))
        self.assertIn("Detected analysis:</strong> Pooled sequencing results (Sprint 3 + 5)", page)
        result = self.run_reviewed(token)
        self.assertIn("Pooled samples:</strong> 9", result)
        self.assertIn("Exploratory alternative", result)  # no A/C/G/T counts: not the paper's PCoA
        page, token = self.review("pooled", file_snv=example("synthetic_pooled_snvs.csv"),
                                  file_reference=example("synthetic_reference.fasta"),
                                  file_pooled_metadata=(io.BytesIO(POOLED_METADATA), "pooled_meta.csv"),
                                  file_positions=example("synthetic_pooled_base_counts.csv"),
                                  file_signatures=example("synthetic_variant_signatures.csv"))
        self.assertEqual(page.count("Your choice"), 5)
        result = self.run_reviewed(token)
        self.assertIn("Paper method: Yue-Clayton PCoA", result)

    # Combined (Sprint 2 + 3) -----------------------------------------------------------
    def test_combined_mode_reports_missing_files_and_runs_with_all(self):
        page, _ = self.review("combined", file_fasta=example("synthetic_individual_sequences.fasta"),
                              file_metadata=example("synthetic_individual_metadata.csv"),
                              file_snv=example("synthetic_pooled_snvs.csv"))
        self.assertIn("these files were not found: reference genome FASTA", page)
        page, token = self.review("combined", file_fasta=example("synthetic_individual_sequences.fasta"),
                                  file_metadata=example("synthetic_individual_metadata.csv"),
                                  file_reference=example("synthetic_reference.fasta"),
                                  file_snv=example("synthetic_pooled_snvs.csv"),
                                  file_positions=example("synthetic_pooled_base_counts.csv"))
        self.assertIn("Detected analysis:</strong> Combined (Sprint 2 + 3 + 5)", page)
        result = self.run_reviewed(token)
        self.assertIn("Variants observed in both sources:</strong> 2", result)
        self.assertIn('id="shared-variant-example"', result)

    # Wrong or ambiguous roles ----------------------------------------------------------
    def test_a_file_in_the_wrong_field_is_a_conflict_not_a_guess(self):
        page, _ = self.review("individual", file_fasta=example("synthetic_individual_metadata.csv"),
                              file_metadata=example("synthetic_individual_sequences.fasta"))
        self.assertEqual(page.count("Conflict — please choose"), 2)
        self.assertIn("Uploaded as &#39;Individual sequences (FASTA)&#39;, but the content looks like", page)
        self.assertIn("Choose a role for this file", page)
        self.assertIn("disabled>Run analysis", page)

    def test_declared_role_is_kept_when_the_content_fits(self):
        # A one-sequence FASTA described as a reference would be detected as the
        # reference; uploaded in the reference field it is simply kept.
        info = apply_declared_role({**empty_info("x.fasta"), "kind": "FASTA", "candidates": ["reference", "fasta"]},
                                   "fasta")
        self.assertEqual((info["suggested_role"], info["role_status"]), ("fasta", "declared"))
        info = apply_declared_role({**empty_info("x.csv"), "kind": "table", "candidates": ["snv"]}, "metadata")
        self.assertEqual((info["suggested_role"], info["role_status"]), (None, "conflict"))
        page, _ = self.review("pooled", file_reference=example("synthetic_reference.fasta"),
                              file_snv=example("synthetic_pooled_snvs.csv"))
        self.assertEqual(page.count("Your choice"), 2)

    def test_ambiguous_files_in_the_automatic_box_still_need_a_choice(self):
        lone = (io.BytesIO(b">S1 sample\nACGTACGTACGTACGTACGT\n"), "one.fasta")
        page, _ = self.review("auto", files=[lone])
        self.assertIn("Please choose", page)
        self.assertIn("could be the reference genome or a single individual sequence", page)
        # The same file in a role field is not ambiguous.
        page, _ = self.review("individual", file_reference=(io.BytesIO(b">S1 sample\nACGTACGTACGTACGTACGT\n"), "one.fasta"))
        self.assertIn("Your choice", page)

    def test_automatic_box_and_role_fields_are_combined(self):
        page, token = self.review("auto", files=[example("synthetic_individual_sequences.fasta")],
                                  file_metadata=example("synthetic_individual_metadata.csv"))
        self.assertIn("Detected</span>", page)
        self.assertIn("Your choice", page)
        self.assertIn("Detected analysis:</strong> Individual sequences (Sprint 2)", page)


if __name__ == "__main__":
    unittest.main()
