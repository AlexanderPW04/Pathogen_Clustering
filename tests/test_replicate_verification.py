"""Phase 3 audit of technical-replicate identification (false-positive replicate claims).

A shared Specimen_ID, date, and location are not enough: the two samples must
document distinct libraries (Library_ID), and exact duplicate records are not
replicates. Undocumented pairs are downgraded to "unverified" and get no rho.
"""

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from genomics import overdispersion as od
from genomics.overdispersion_analysis import run_overdispersion_analysis
from tests.test_overdispersion import evidence_rows, phase1_like_evidence
from wastewater.processed_snv import load_processed_snvs

EXAMPLES = Path(__file__).resolve().parents[1] / "examples" / "integrated"


def pair_problems(table):
    return {specimen: problem for specimen, _, problem in od.replicate_pairs(table)}


class ReplicateClaimTest(unittest.TestCase):
    def base(self):
        rows = []
        for sample, alt, depth in [("A1", 100, 500), ("A2", 112, 480)]:
            rows.append({"Sample_ID": sample, "Specimen_ID": "A", "Position": 5, "Alt": "G",
                         "Alt_Read_Count": alt, "Total_Depth": depth})
        return evidence_rows(rows)

    def test_documented_distinct_libraries_are_a_usable_pair(self):
        self.assertEqual(pair_problems(self.base())["A"], "")

    def test_matching_metadata_without_library_ids_is_unverified(self):
        table = self.base().drop(columns=["Library_ID"])
        self.assertEqual(pair_problems(table)["A"], od.MISSING_LIBRARY)
        blank = self.base()
        blank.loc[1, "Library_ID"] = ""
        self.assertEqual(pair_problems(blank)["A"], od.MISSING_LIBRARY)

    def test_same_library_listed_twice_is_not_a_replicate(self):
        table = self.base()
        table["Library_ID"] = "LIB_1"
        self.assertEqual(pair_problems(table)["A"], od.SAME_LIBRARY)

    def test_identical_records_are_duplicates_not_replicates(self):
        table = self.base()
        table.loc[1, ["Alt_Read_Count", "Total_Depth"]] = [100, 500]
        self.assertEqual(pair_problems(table)["A"], od.DUPLICATE_RECORDS)

    def test_contradicting_metadata_still_wins_over_library_ids(self):
        table = self.base()
        table.loc[1, "Collection_Date"] = pd.Timestamp("2026-03-09")
        self.assertIn("different collection dates", pair_problems(table)["A"])


class CalibrationDowngradeTest(unittest.TestCase):
    def test_unverified_pair_gets_no_rho_and_a_reason(self):
        evidence = phase1_like_evidence(["A1", "A2"], specimen="SPEC").drop(columns=["Library_ID"])
        out = run_overdispersion_analysis(evidence, None, [0.02])
        row = out["calibration"].iloc[0]
        self.assertEqual(row["Status"], od.UNVERIFIED_PAIR)
        self.assertTrue(out["evidence"]["Rho"].isna().all())
        self.assertEqual(out["summary"]["calibrated_specimens"], 0)
        self.assertIn("could not be estimated", out["summary"]["headline"])
        self.assertTrue(any("unverified" in reason for reason in out["summary"]["reasons"]))

    def test_duplicated_library_file_is_not_calibrated(self):
        once = phase1_like_evidence(["A1"], specimen="SPEC")
        twice = pd.concat([once, once.assign(Sample_ID="A2", Library_ID="LIB_A2")], ignore_index=True)
        out = run_overdispersion_analysis(twice, None, [])
        self.assertEqual(out["calibration"].iloc[0]["Status_Detail"], od.DUPLICATE_RECORDS)
        self.assertTrue(out["evidence"]["Rho"].isna().all())

    def test_documented_pair_is_still_calibrated(self):
        evidence = phase1_like_evidence(["A1", "A2"], specimen="SPEC")
        out = run_overdispersion_analysis(evidence, None, [])
        self.assertIn(out["calibration"].iloc[0]["Status"], od.USABLE_STATUSES)


class LibraryLoaderTest(unittest.TestCase):
    def test_library_id_is_read_from_the_table_or_metadata(self):
        snvs = pd.read_csv(EXAMPLES / "synthetic_replicate_snvs.csv")
        self.assertTrue(snvs["Library_ID"].ne("").all())
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "snvs.csv"
            snvs.to_csv(path, index=False)
            data, _ = load_processed_snvs(path)
            self.assertEqual(data.groupby("Sample_ID")["Library_ID"].first()["REP_A2"], "LIB_REP_A2")
            meta = Path(folder) / "meta.csv"
            snvs[["Sample_ID", "Library_ID"]].drop_duplicates().to_csv(meta, index=False)
            snvs.drop(columns=["Library_ID"]).to_csv(path, index=False)
            data, _ = load_processed_snvs(path, meta)
            self.assertEqual(data.groupby("Sample_ID")["Library_ID"].first()["REP_B1"], "LIB_REP_B1")
            bad = snvs.copy()
            bad.loc[bad.index[0], "Library_ID"] = "OTHER"
            bad.to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, "one Library_ID"):
                load_processed_snvs(path)

    def test_phase2_style_file_without_library_ids_is_unverified(self):
        snvs = pd.read_csv(EXAMPLES / "synthetic_replicate_snvs.csv").drop(columns=["Library_ID"])
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "snvs.csv"
            snvs.to_csv(path, index=False)
            data, _ = load_processed_snvs(path)
        self.assertNotIn("Library_ID", data)
        evidence = data.rename(columns={"Alt_Depth": "Alt_Read_Count", "Read_Depth": "Total_Depth",
                                         "Reference_Genome": "Reference_ID"})
        evidence["Variant_Key"] = (evidence["Reference_ID"] + ":" + evidence["Position"].astype(str) + ":"
                                   + evidence["Ref"] + ">" + evidence["Alt"])
        problems = pair_problems(evidence)
        self.assertEqual(problems["SPEC_A"], od.MISSING_LIBRARY)
        self.assertIn("different collection dates", problems["SPEC_D"])


if __name__ == "__main__":
    unittest.main()
