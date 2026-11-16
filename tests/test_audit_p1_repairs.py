"""Regression tests for the four P1 repairs of the comprehensive Sprint 3 + 5 audit.

docs/SPRINT3_SPRINT5_COMPREHENSIVE_AUDIT.md section 14:
  P1-1  reading an upload is bounded (cell and column limits; one review at a time)
  P1-2  downloaded tables never contain live spreadsheet formulas
  P1-3  what the analysis would refuse is shown at review, and a failed run keeps the uploads
  P1-4  segment/chromosome columns keep records on their own sequence

All data are small SYNTHETIC tables written by the tests.
"""

import io
import os
import re
import threading
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch

import openpyxl
import pandas as pd

import app as app_module
from app import STAGED, app
from genomics import input_inspection
from genomics.excel_tables import WorkbookTooLarge
from genomics.file_formats import read_table
from genomics.input_inspection import ServerBusy, inspect_file
from genomics.spreadsheet_safety import escape_csv_bytes, escape_frame, escape_text
from tests.test_phase4a_import import REFERENCE_SEQUENCE, Workspace, long_table

FORMULA = '=HYPERLINK("https://example.invalid/","click")'


def page_text(page):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", page))


class SpreadsheetSafetyTest(unittest.TestCase):
    """P1-2."""

    def test_only_formula_like_text_is_escaped(self):
        for value in [FORMULA, "@SUM(A1)", "+cmd|x", "-cmd|x", "\tx", "\rx"]:
            self.assertEqual(escape_text(value), "'" + value)
        for value in ["S1", "-", "+", "-0.25", "+3", "-1e-5", "", "C10029T", 3, -2.5, None]:
            self.assertEqual(escape_text(value), value)

    def test_frames_keep_numbers_strands_and_untouched_tables(self):
        table = pd.DataFrame({"Sample_ID": [FORMULA, "S2"], "Strand": ["-", "+"], "Shift": ["-0.25", "1"],
                              "Value": [-1.5, 2.0]})
        escaped = escape_frame(table)
        self.assertEqual(list(escaped["Sample_ID"]), ["'" + FORMULA, "S2"])
        for column in ["Strand", "Shift", "Value"]:
            self.assertEqual(list(escaped[column]), list(table[column]))
        self.assertEqual(table.loc[0, "Sample_ID"], FORMULA)  # the analysis table itself is not changed
        clean = pd.DataFrame({"a": ["x", "-"], "b": [1, -2]})
        self.assertIs(escape_frame(clean), clean)
        self.assertEqual(escape_csv_bytes(b"a,b\nx,-1\n"), b"a,b\nx,-1\n")
        self.assertIn(b"'=1+1", escape_csv_bytes(b"a,b\n=1+1,2\n"))


class SafeDownloadsWebTest(Workspace):
    """P1-2 end to end: a formula in a sample ID reaches no download as a formula."""

    def setUp(self):
        super().setUp()
        os.environ["AI_ENRICHMENT_ENABLED"] = "false"
        self.client = app.test_client()
        self.addCleanup(STAGED.clear)

    def test_csv_and_workbook_downloads_hold_text_not_formulas(self):
        table = long_table()
        table.loc[table["Sample_ID"].eq("S1"), "Sample_ID"] = FORMULA
        files = [(io.BytesIO(table.to_csv(index=False).encode()), "snv.csv"),
                 (io.BytesIO((">SYN.1 SYNTHETIC\n" + REFERENCE_SEQUENCE + "\n").encode()), "ref.fasta")]
        page = self.client.post("/review", data={"files": files}, content_type="multipart/form-data")
        token = re.search(r'name="token" value="([^"]+)"', page.get_data(as_text=True)).group(1)
        with redirect_stdout(io.StringIO()):
            result = self.client.post("/run", data={"token": token, "role_0": "snv", "role_1": "reference"})
        html = result.get_data(as_text=True)
        links = dict((name, link) for link, name in re.findall(r'href="(/results/[^"]+/([^"/]+))"', html))
        observations = self.client.get(links["canonical_variant_observations.csv"]).get_data(as_text=True)
        self.assertIn("'" + FORMULA.replace('"', '""'), observations)
        self.assertNotIn("\n" + FORMULA.replace('"', '""'), observations)
        workbook = openpyxl.load_workbook(io.BytesIO(self.client.get(links["analysis_results.xlsx"]).data))
        formulas = [cell.coordinate for sheet in workbook for row in sheet.iter_rows() for cell in row
                    if cell.data_type == "f"]
        self.assertEqual(formulas, [])


class UploadBoundsTest(Workspace):
    """P1-1."""

    def test_wide_and_large_text_tables_are_refused_before_parsing(self):
        wide = self.write_text("wide.csv", ",".join(f"c{i}" for i in range(50)) + "\n" + ",".join("1" * 50) + "\n")
        with patch.dict(os.environ, {"MAX_TABLE_COLUMNS": "40"}), self.assertRaises(WorkbookTooLarge) as raised:
            read_table(wide)
        self.assertIn("50 columns", str(raised.exception))
        tall = self.write_text("tall.csv", "a,b\n" + "1,2\n" * 100)
        with patch.dict(os.environ, {"MAX_WORKSHEET_CELLS": "150"}), self.assertRaises(WorkbookTooLarge):
            read_table(tall)
        with patch.dict(os.environ, {"MAX_WORKSHEET_CELLS": "200"}):
            self.assertEqual(read_table(tall).shape, (100, 2))  # within the limit: unchanged
        with patch.dict(os.environ, {"MAX_WORKSHEET_CELLS": "150"}):
            info = inspect_file(tall, "tall.csv")
        self.assertTrue(info["refused"])
        self.assertIn("larger than this server can process safely", info["problems"][0])

    def test_rows_wider_than_the_header_give_the_same_worksheet(self):
        path = self.folder / "wide_rows.xlsx"
        book = openpyxl.Workbook()
        sheet = book.active
        sheet.append(["Sample_ID", "Position", "Ref", "Alt"])
        sheet.append(["S1", 2, "C", "T", "extra", "cells", "beyond", "the", "header"])
        sheet.append(["S2", 7, "G", "A"])
        book.save(path)
        table = read_table(path)
        # The worksheet's own width (9 columns here) is kept; cells are never shifted or dropped.
        self.assertEqual(list(table.columns[:4]), ["Sample_ID", "Position", "Ref", "Alt"])
        self.assertEqual(table.iloc[:, :4].values.tolist(), [["S1", "2", "C", "T"], ["S2", "7", "G", "A"]])
        self.assertEqual(table.iloc[0, 4:].tolist(), ["extra", "cells", "beyond", "the", "header"])

    def test_a_second_review_waits_and_is_then_told_the_server_is_busy(self):
        held, release = threading.Event(), threading.Event()

        def hold():
            with input_inspection.INSPECTION_LOCK:
                held.set()
                release.wait(5)

        worker = threading.Thread(target=hold)
        worker.start()
        held.wait(5)
        try:
            with patch.object(input_inspection, "INSPECTION_WAIT_SECONDS", 0.05), self.assertRaises(ServerBusy):
                inspect_file(self.reference, "ref.fasta")
        finally:
            release.set()
            worker.join()
        self.assertEqual(inspect_file(self.reference, "ref.fasta")["kind"], "FASTA sequences")


class ReviewShowsRunProblemsTest(Workspace):
    """P1-3."""

    def test_alternative_reads_above_depth_block_the_review(self):
        table = long_table()
        table.loc[0, "Alt_Depth"] = 150
        review, result = self.run_files([{"path": self.write("snv.csv", table), "role": "snv"},
                                         {"path": self.reference, "role": "reference"}])
        self.assertFalse(review["ready"])
        self.assertIsNone(result)
        self.assertTrue(any("would refuse" in text and "larger than the total read depth" in text
                            for text in review["problems"]))

    def test_unreadable_metadata_dates_are_named_at_review(self):
        variants = long_table().drop(columns=["Collection_Date", "Location_Name"])
        metadata = pd.DataFrame({"Sample_ID": ["S1", "S2", "S3"], "date": ["5-19", "5-28", "6-09"]})
        review, result = self.run_files([{"path": self.write("snv.csv", variants), "role": "snv"},
                                         {"path": self.write("meta.csv", metadata), "role": "pooled_metadata"},
                                         {"path": self.reference, "role": "reference"}])
        self.assertFalse(review["ready"])
        problem = next(text for text in review["problems"] if "would refuse" in text)
        self.assertIn("'5-19'", problem)
        self.assertIn("never guessed", problem)


class FailedRunKeepsUploadsTest(unittest.TestCase):
    """P1-3: when a run produces nothing, the review is shown again and the files stay."""

    def setUp(self):
        os.environ["AI_ENRICHMENT_ENABLED"] = "false"
        self.client = app.test_client()
        self.addCleanup(STAGED.clear)

    def test_failed_run_returns_to_step_2_with_the_files(self):
        table = long_table()
        files = [(io.BytesIO(table.to_csv(index=False).encode()), "snv.csv"),
                 (io.BytesIO((">SYN.1 SYNTHETIC\n" + REFERENCE_SEQUENCE + "\n").encode()), "ref.fasta")]
        page = self.client.post("/review", data={"files": files}, content_type="multipart/form-data")
        token = re.search(r'name="token" value="([^"]+)"', page.get_data(as_text=True)).group(1)
        form = {"token": token, "role_0": "snv", "role_1": "reference"}
        with patch.object(app_module, "run_integrated_analysis", side_effect=ValueError("Synthetic failure.")), \
                redirect_stdout(io.StringIO()):
            page = self.client.post("/run", data=form).get_data(as_text=True)
        words = page_text(page)
        self.assertIn("The analysis could not be completed: Synthetic failure.", words)
        self.assertIn("Step 3 · Analyses these data support", words)
        self.assertIn(token, STAGED)
        with redirect_stdout(io.StringIO()):
            page = self.client.post("/run", data=form).get_data(as_text=True)
        self.assertIn("Results workbook (XLSX)", page)  # the same files run once the cause is gone
        self.assertNotIn(token, STAGED)


class SequenceColumnsTest(Workspace):
    """P1-4."""

    def segmented(self, column):
        table = long_table().drop(columns=["Reference_Genome"])
        table.insert(1, column, ["SYN.1", "SYN.2", "SYN.1", "SYN.2", "SYN.1", "SYN.2"])
        return table

    def test_a_segment_column_is_the_reference_column(self):
        review, result = self.run_files([{"path": self.write("snv.csv", self.segmented("Segment")), "role": "snv",
                                          "options": {"declared_reference": "SYN.1"}},
                                         {"path": self.reference, "role": "reference"}])
        self.assertTrue(review["ready"])
        verification = result["tables"]["reference_verification.csv"].set_index("Reference_ID")
        self.assertEqual(verification.loc["SYN.1", "Records"], 3)
        self.assertEqual(verification.loc["SYN.2", "Status"], "Supplied but unverified")
        self.assertEqual(set(result["tables"]["canonical_variant_observations.csv"]["Reference_ID"]),
                         {"SYN.1", "SYN.2"})

    def test_one_stated_reference_is_refused_for_several_sequences(self):
        review, result = self.run_files([{"path": self.write("snv.csv", self.segmented("Molecule")), "role": "snv",
                                          "options": {"declared_reference": "SYN.1"}},
                                         {"path": self.reference, "role": "reference"}])
        self.assertFalse(review["ready"])
        self.assertTrue(any("Column 'Molecule' names 2 different sequences" in text for text in review["problems"]))


if __name__ == "__main__":
    unittest.main()
