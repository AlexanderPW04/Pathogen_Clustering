"""Optional real-data checks of the Phase 4A importer on published wastewater datasets.

The published files are not in Git. Download them (see the retrieval notes in
validation/check_published_datasets.py) and set any of:

    FONTENELE_MMC2=<path to Fontenele et al. 2021 mmc2.xlsx>
    CRITS_CHRISTOPH_DIR=<folder with mBio.02703-20-st001.xlsx, -st004.xlsx, reference.fasta,
                         and optionally the GitHub inStrain *_SNVs.tsv and Wastewater_SNVs_v3.csv>
    SCHUMANN_DIR=<folder with Schumann et al. 2022 mmc1.xlsx (and optionally mmc2.xlsx, mmc3.xlsx)>

Each test is skipped when its data are absent; nothing is substituted.
"""

import os
import tempfile
import unittest
from pathlib import Path

from validation.check_published_datasets import crits_christoph, fontenele, schumann

MMC2 = os.environ.get("FONTENELE_MMC2")
CRITS = os.environ.get("CRITS_CHRISTOPH_DIR")
SCHUMANN = os.environ.get("SCHUMANN_DIR")


class PublishedDatasetTest(unittest.TestCase):
    def assert_all_pass(self, report):
        failed = [row for row in report["checks"] if row["Result"] != "PASS"]
        self.assertFalse(failed, failed)
        self.assertTrue(report["checks"])

    @unittest.skipUnless(MMC2 and Path(MMC2).exists(), "Set FONTENELE_MMC2")
    def test_fontenele_regression_dataset(self):
        self.assert_all_pass(fontenele(MMC2, tempfile.mkdtemp()))

    @unittest.skipUnless(CRITS and Path(CRITS).exists(), "Set CRITS_CHRISTOPH_DIR")
    def test_crits_christoph_base_counts_and_unmatched_metadata(self):
        self.assert_all_pass(crits_christoph(CRITS, tempfile.mkdtemp()))

    @unittest.skipUnless(SCHUMANN and Path(SCHUMANN).exists(), "Set SCHUMANN_DIR")
    def test_schumann_wide_matrix(self):
        self.assert_all_pass(schumann(SCHUMANN, tempfile.mkdtemp(), include_large=False))


if __name__ == "__main__":
    unittest.main()
