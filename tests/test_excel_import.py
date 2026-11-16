"""Sprint 3 Excel uploads: every workbook either runs or gets a controlled, specific message.

Workbooks are generated here with openpyxl (invented values). The layouts copy
the published Fontenele et al. (2021) supplements: Table 2 has a title row
above its headers, Table 3 has two rows above and lacks frequencies and depths.
The real files are checked separately in tests/test_optional_fontenele_table.py.
"""

import io
import os
import re
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import app as application  # noqa: E402
from app import STAGED, app  # noqa: E402
from genomics.column_detection import explain_snv_columns  # noqa: E402
from genomics.excel_tables import WorkbookTooLarge, describe_sheets, find_header_row, read_sheet  # noqa: E402
from wastewater.processed_snv import PAPER_FORMAT, load_processed_snvs  # noqa: E402

HEADERS = ["Sample_ID", "Reference_Genome", "Position", "Ref", "Alt", "Alt_Frequency", "Read_Depth",
           "Collection_Date", "Location_Name"]
ROWS = [["W1", "NC_045512.2", 241, "C", "T", 0.25, 100, "2026-03-01", "Site A"],
        ["W1", "NC_045512.2", 3037, "C", "T", 0.5, 200, "2026-03-01", "Site A"],
        ["W2", "NC_045512.2", 241, "C", "T", 0.75, 300.0, "2026-03-08", "Site B"]]
PAPER_HEADERS = ["SAMPLE", "COLLECTION_DATE", "STATE", "LOCATION", "POSITION", "REF", "ALT", "ALT_FREQ", "TOTAL_DP"]
PAPER_TITLE = ("Supplementary Table 2: SNVs in 52 wastewater samples of SARS-CoV-2 "
               "(invented rows for this test).")


def workbook_bytes(sheets):
    """{sheet name: rows} as an .xlsx file in memory."""

    book = openpyxl.Workbook()
    book.remove(book.active)
    for name, rows in sheets.items():
        sheet = book.create_sheet(name)
        for row in rows:
            sheet.append(row)
    buffer = io.BytesIO()
    book.save(buffer)
    return buffer.getvalue()


def text(page):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", page)).replace("&#39;", "'")


class WorkbookReadingTest(unittest.TestCase):
    def path_for(self, sheets):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        path = Path(folder.name) / "book.xlsx"
        path.write_bytes(workbook_bytes(sheets))
        return path

    def test_cells_stay_text_as_written(self):
        table = read_sheet(self.path_for({"Data": [HEADERS] + ROWS}))
        self.assertEqual(list(table.columns), HEADERS)
        self.assertEqual(table.loc[2, "Read_Depth"], "300")  # a whole-number float loses its .0
        self.assertEqual(table.loc[0, "Alt_Frequency"], "0.25")
        self.assertEqual(table.attrs["header_row"], 0)

    def test_header_row_is_found_below_a_title(self):
        rows = [["A title", None], [None, None], HEADERS] + ROWS
        table = read_sheet(self.path_for({"Data": rows}))
        self.assertEqual(table.attrs["header_row"], 2)
        self.assertEqual(len(table), 3)
        self.assertIn("A title", table.attrs["preamble"])

    def test_blank_rows_are_skipped_and_headers_stay_unique(self):
        rows = [["Sample_ID", "Date", "Date", None], ["W1", "a", "b", "c"], [None, None, None, None], ["W2", "d", "e", ""]]
        table = read_sheet(self.path_for({"Data": rows}))
        self.assertEqual(list(table.columns), ["Sample_ID", "Date", "Date.1", "Unnamed: 3"])
        self.assertEqual(list(table["Sample_ID"]), ["W1", "W2"])

    def test_first_row_is_the_header_when_nothing_is_recognized(self):
        self.assertEqual(find_header_row([["Patient", "Height"], ["1", "2"]]), 0)

    def test_oversized_worksheet_is_refused_before_it_fills_memory(self):
        path = self.path_for({"Data": [HEADERS] + ROWS * 40})
        with patch.dict(os.environ, {"MAX_WORKSHEET_CELLS": "100"}), self.assertRaises(WorkbookTooLarge) as raised:
            read_sheet(path)
        self.assertIn("larger than this server can process safely", str(raised.exception))

    def test_describe_sheets_reports_recognized_headers(self):
        sheets = describe_sheets(self.path_for({"Notes": [["free text"]], "Calls": [HEADERS] + ROWS}))
        self.assertEqual([sheet["name"] for sheet in sheets], ["Notes", "Calls"])
        self.assertGreater(sheets[1]["score"], sheets[0]["score"])

    def test_published_layout_is_recognized_and_keeps_its_rules(self):
        rows = [[None, None, None, None, PAPER_TITLE], PAPER_HEADERS,
                ["348", "9_Jun_2020", "AZ", "TP01", 201, "t", "C", 0.0485437, 618],
                ["271", "3_Jun_2020", "AZ", "TP06", 202, "t", "A", 0.14, 50]]
        data, _ = load_processed_snvs(self.path_for({"All detected SNVs per sample": rows}))
        self.assertEqual(data.attrs["input_format"], PAPER_FORMAT)
        self.assertEqual(list(data["Position"]), [201, 202])
        self.assertEqual(list(data["Ref"]), ["T", "T"])  # the workbook's lower-case bases
        self.assertEqual(data.loc[0, "Location_Name"], "AZ: TP01")

    def test_invalid_values_name_the_first_bad_cell(self):
        rows = [HEADERS, ROWS[0], ["W1", "NC_045512.2", "abc", "C", "T", 0.5, 200, "2026-03-01", "Site A"]]
        with self.assertRaises(ValueError) as raised:
            load_processed_snvs(self.path_for({"Data": rows}))
        self.assertIn("1 position(s) are not positive integers (first at row 3: 'abc')", str(raised.exception))

    def test_explanations_match_what_was_detected(self):
        self.assertIn("does not contain the required variant-position column",
                      explain_snv_columns(["Sample_ID", "Reads", "Coverage"]))
        self.assertIn("but not the processed SNV values", explain_snv_columns(["POS", "REF", "ALT", "Gene"]))
        self.assertIsNone(explain_snv_columns(HEADERS))


class UploadTest(unittest.TestCase):
    def setUp(self):
        patches = [patch.dict(os.environ, {"GEONAMES_USERNAME": ""})]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        self.addCleanup(STAGED.clear)
        self.client = app.test_client()

    def review(self, content, name="book.xlsx", **form):
        response = self.client.post("/review", data={"files": [(io.BytesIO(content), name)], **form},
                                    content_type="multipart/form-data")
        page = response.get_data(as_text=True)
        token = re.search(r'name="token" value="([^"]+)"', page)
        return page, token.group(1) if token else None

    def run_review(self, token, **form):
        with redirect_stdout(io.StringIO()):
            return self.client.post("/run", data={"token": token, **form}).get_data(as_text=True)

    def test_valid_workbook_and_csv_give_the_same_result(self):
        page, token = self.review(workbook_bytes({"Calls": [HEADERS] + ROWS}))
        self.assertNotIn("disabled>Run analysis", page)
        self.assertIn("Preview of book.xlsx", page)
        from_excel = self.run_review(token)
        csv = "\n".join(",".join(str(value) for value in row) for row in [HEADERS] + ROWS).encode()
        _, token = self.review(csv, "calls.csv")
        from_csv = self.run_review(token)
        for page in (from_excel, from_csv):
            self.assertIn("2 pooled samples analyzed", page)
            self.assertNotIn("Analysis could not be completed", page)

    def test_unsupported_workbook_says_which_values_are_missing(self):
        rows = [["Supplementary Table 3. Novel SNVs"], ["A) SNV information"],
                ["POS", "Gene", "REF", "ALT", "Wastewater Sample"], [279, "ORF1A", "t", "C", 253]]
        page, _ = self.review(workbook_bytes({"novel": rows}))
        words = text(page)
        self.assertIn("contains variant positions but not the processed SNV values", words)
        self.assertIn("allele frequency (Alt_Frequency), read depth (Read_Depth)", words)
        self.assertIn("another table may hold the per-sample values", words)
        self.assertNotIn("Traceback", page)

    def test_several_data_worksheets_are_a_question_not_a_guess(self):
        other = [HEADERS, ["X9", "NC_045512.2", 500, "A", "G", 0.5, 40, "2026-04-01", "Site C"]]
        content = workbook_bytes({"Run 1": [HEADERS] + ROWS, "Run 2": other})
        page, token = self.review(content)
        self.assertIn("Multiple worksheets were detected", text(page))
        self.assertIn('name="sheet_0"', page)
        response = self.client.post("/review", data={"token": token, "sheet_0": "Run 2"})
        page = response.get_data(as_text=True)
        self.assertNotIn("Multiple worksheets were detected", text(page))
        self.assertIn("Worksheet 'Run 2' was read", text(page))
        self.assertIn("X9", page)

    def test_a_notes_sheet_does_not_get_in_the_way(self):
        page, _ = self.review(workbook_bytes({"Read me": [["About this file"]], "Calls": [HEADERS] + ROWS}))
        self.assertIn("Worksheet 'Calls' was read", text(page))
        self.assertNotIn("Multiple worksheets were detected", text(page))

    def test_missing_columns_on_the_snv_page(self):
        rows = [[name for name in HEADERS if name != "Alt_Frequency"], ["W1", "NC_045512.2", 241, "C", "T", 100, "", ""]]
        response = self.client.post("/wastewater", data={"snv_file": (io.BytesIO(workbook_bytes({"Calls": rows})), "a.xlsx")},
                                    content_type="multipart/form-data")
        self.assertIn("allele frequency (Alt_Frequency)", text(response.get_data(as_text=True)))

    def test_invalid_numbers_are_caught_at_review(self):
        rows = [HEADERS, ["W1", "NC_045512.2", "abc", "C", "T", "high", 100, "2026-03-01", "Site A"]]
        page, _ = self.review(workbook_bytes({"Calls": rows}))
        words = text(page)
        self.assertIn("Column 'Position' should hold variant positions", words)
        self.assertIn("Column 'Alt_Frequency' should hold allele frequencies", words)
        self.assertIn("disabled>Run analysis", page)

    def test_missing_dates_and_sites_do_not_stop_the_analysis(self):
        columns = HEADERS[:7]
        page, token = self.review(workbook_bytes({"Calls": [columns] + [row[:7] for row in ROWS]}))
        self.assertNotIn("disabled>Run analysis", page)
        self.assertIn("2 pooled samples analyzed", self.run_review(token))

    def test_a_dataset_too_large_for_this_server_is_refused_clearly_and_not_trimmed(self):
        content = workbook_bytes({"Calls": [HEADERS] + ROWS})
        with patch.dict(os.environ, {"MAX_OBSERVATION_CELLS": "3"}):
            page, _ = self.review(content)
        words = text(page)
        self.assertIn("exceeds the available memory for this server: 2 variants x 2 samples is 4 combinations", words)
        self.assertIn("No samples or variants are removed automatically", words)
        self.assertIn("disabled>Run analysis", page)
        with patch.dict(os.environ, {"MAX_OBSERVATION_CELLS": "4"}):
            page, _ = self.review(content)
        self.assertNotIn("exceeds the available memory", text(page))

    def test_running_out_of_memory_is_a_message_and_the_server_stays_up(self):
        _, token = self.review(workbook_bytes({"Calls": [HEADERS] + ROWS}))
        with patch("app.run_integrated_analysis", side_effect=MemoryError()), \
                self.assertLogs(app.logger, "ERROR"):
            page = self.run_review(token)
        self.assertIn("exceeds the available memory for this server", text(page))
        self.assertEqual(self.client.get("/health").status_code, 200)
        self.assertTrue(application.ANALYSIS_SLOT.acquire(blocking=False))  # the slot was released
        application.ANALYSIS_SLOT.release()

    def test_unexpected_analysis_errors_are_logged_and_stay_private(self):
        _, token = self.review(workbook_bytes({"Calls": [HEADERS] + ROWS}))
        with patch("app.attach_enrichment", side_effect=KeyError("private detail")), \
                self.assertLogs(app.logger, "ERROR") as logs:
            page = self.run_review(token)
        self.assertIn("unexpected server error (KeyError)", page)
        self.assertNotIn("private detail", page)
        self.assertIn("private detail", "\n".join(logs.output))  # kept for the developer

    def test_one_analysis_at_a_time_and_the_files_are_kept_for_a_retry(self):
        _, token = self.review(workbook_bytes({"Calls": [HEADERS] + ROWS}))
        self.assertTrue(application.ANALYSIS_SLOT.acquire(blocking=False))
        try:
            page = self.run_review(token)
        finally:
            application.ANALYSIS_SLOT.release()
        self.assertIn("Another analysis is running on this server", text(page))
        self.assertIn("2 pooled samples analyzed", self.run_review(token))  # same upload, no re-upload

    def test_failure_pages_are_readable_and_offer_a_way_back(self):
        with patch.dict(app.config, {"MAX_CONTENT_LENGTH": 10}):
            response = self.client.post("/review", data={"files": [(io.BytesIO(b"x" * 100), "a.csv")]},
                                        content_type="multipart/form-data")
        self.assertEqual(response.status_code, 413)
        self.assertIn("Upload too large", response.get_data(as_text=True))
        self.assertIn('href="/"', response.get_data(as_text=True))
        self.assertEqual(self.client.get("/health").status_code, 200)

    def test_a_corrupt_workbook_is_refused_on_its_own(self):
        page, _ = self.review(b"PK\x03\x04 this is not a real workbook", "broken.xlsx")
        self.assertIn("Could not read this file", text(page))
        self.assertEqual(self.client.get("/health").status_code, 200)


if __name__ == "__main__":
    unittest.main()
