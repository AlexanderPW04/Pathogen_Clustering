"""Tests for the optional offline-output converters (synthetic text only).

These do not run BWA, SAMtools, or iVar; they check the Python converters
that read those tools' text outputs.
"""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "offline_workflow"))
from convert_offline_outputs import (  # noqa: E402
    combine_ivar_tables,
    count_pileup_bases,
    pileup_to_base_counts,
)

from genomics.pooled_observations import load_position_table  # noqa: E402
from wastewater.processed_snv import IVAR_FORMAT, load_processed_snvs  # noqa: E402


class PileupTest(unittest.TestCase):
    def test_markers_indels_and_low_quality_bases(self):
        # '.' and ',' are the reference base A; '^I' starts a read; '$' ends one;
        # 'g' has quality '#' (Q2) and is skipped; '+2AG' is an insertion; '*' a deletion.
        counts = count_pileup_bases("a", ".$,^ICg+2AGt*", "III#II", minimum_quality=20)
        self.assertEqual(counts, {"A": 2, "C": 1, "G": 0, "T": 1})

    def test_pileup_file_becomes_a_valid_position_table(self):
        text = (
            "SYNREF.1\t1\tG\t3\t..A\tIII\n"
            "SYNREF.1\t2\tT\t0\t*\t*\n"
        )
        with tempfile.TemporaryDirectory() as folder:
            pileup = Path(folder) / "s1.mpileup"
            pileup.write_text(text, encoding="utf-8")
            counts = pileup_to_base_counts(pileup, "S1")
            output = Path(folder) / "counts.csv"
            counts.to_csv(output, index=False)
            table = load_position_table(output)
        self.assertEqual(list(table.loc[0, ["A_Count", "C_Count", "G_Count", "T_Count"]]), [1, 0, 2, 0])
        self.assertEqual(table.loc[1, "Depth"], 0)
        self.assertTrue(table.attrs["has_base_counts"])


class IvarCombineTest(unittest.TestCase):
    def test_per_sample_ivar_files_are_combined_and_imported(self):
        header = "REGION\tPOS\tREF\tALT\tREF_DP\tALT_DP\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\n"
        with tempfile.TemporaryDirectory() as folder:
            tsv_folder = Path(folder) / "ivar_tsv"
            tsv_folder.mkdir()
            (tsv_folder / "S1.tsv").write_text(header + "SYNREF.1\t15\tT\tC\t80\t20\t0.2\t100\t0.001\tTRUE\n",
                                               encoding="utf-8")
            (tsv_folder / "S2.tsv").write_text(header + "SYNREF.1\t15\tT\tC\t60\t40\t0.4\t100\t0.001\tTRUE\n",
                                               encoding="utf-8")
            combined = combine_ivar_tables(tsv_folder)
            output = Path(folder) / "pooled_snvs.tsv"
            combined.to_csv(output, sep="\t", index=False)
            data, _ = load_processed_snvs(output)
        self.assertEqual(list(combined["SAMPLE"]), ["S1", "S2"])
        self.assertEqual(data.attrs["input_format"], IVAR_FORMAT)
        self.assertEqual(list(data["Alt_Frequency"]), [0.2, 0.4])


if __name__ == "__main__":
    unittest.main()
