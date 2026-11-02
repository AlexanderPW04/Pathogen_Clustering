"""Optional check against the real Fontenele et al. (2021) Supplementary Table 2.

The workbook (mmc2.xlsx) is not stored in Git. To run this test, download it
from the article's supplementary material and set FONTENELE_MMC2 to its path:

    set FONTENELE_MMC2=C:\\path\\to\\mmc2.xlsx      (Windows)
    export FONTENELE_MMC2=/path/to/mmc2.xlsx        (Linux/macOS)
"""

import os
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from validation.check_fontenele_2021 import compare  # noqa: E402

WORKBOOK = os.environ.get("FONTENELE_MMC2", "")


@unittest.skipUnless(WORKBOOK and Path(WORKBOOK).is_file(),
                     "Optional real-data test: set FONTENELE_MMC2 to the path of mmc2.xlsx.")
class PublishedSupplementTest(unittest.TestCase):
    def test_import_reproduces_the_counts_the_paper_reports(self):
        table, _ = compare(WORKBOOK)
        values = dict(zip(table["Quantity"], table["Application calculates"]))
        self.assertEqual(values["SNV records (all samples)"], 7973)
        self.assertEqual(values["Samples with >90% coverage (analyzed)"], 52)
        self.assertEqual(values["States represented"], 10)
        self.assertEqual(values["SNVs per sample (range)"], "24 to 793")
        self.assertEqual(values["Samples with spike D614G (A23403G)"], 52)
        self.assertEqual(values["Novel SNVs in more than 8 samples"], 4)
        # Without per-position A/C/G/T counts the paper's PCoA is not claimed.
        self.assertTrue(values["Population structure available from this table"].startswith("Exploratory"))


if __name__ == "__main__":
    unittest.main()
