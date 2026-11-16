"""Read Excel workbooks as plain text tables, without guessing silently.

Published workbooks rarely start with a clean table: the Fontenele et al.
(2021) Supplementary Table 2 has a title row above its headers, and its
Table 3 has two. This module streams one worksheet in a single pass (no second
parse, little memory), finds the header row by counting recognized column
names, and returns every cell as text, exactly as pandas would have shown it,
so nothing is converted or guessed.

When several worksheets look like data, nothing is chosen for the user: the
review page asks which one to use.
"""

import os
from datetime import date, datetime, time

import openpyxl
import pandas as pd

from genomics.column_detection import ROLE_FIELDS, squash

HEADER_SCAN_ROWS = 15
PREVIEW_SHEET_ROWS = 40
MIN_DATA_HEADERS = 3  # recognized headers that make a worksheet "look like data"
KNOWN_HEADERS = {alias for fields in ROLE_FIELDS.values() for spec in fields for alias in spec["aliases"]}


def max_cells():
    """Largest worksheet accepted (cells); larger tables would not fit on a small server."""

    return int(os.getenv("MAX_WORKSHEET_CELLS", "4000000"))


class WorkbookTooLarge(ValueError):
    """The worksheet is larger than this server can process safely."""


def cell_text(value):
    """One cell as text: blanks stay blank, and whole-number floats lose their '.0'."""

    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    if isinstance(value, (datetime, date, time)):
        return str(value)
    return str(value)


def header_score(cells):
    """How many cells of a row are recognized column names."""

    return sum(squash(cell) in KNOWN_HEADERS for cell in cells if cell)


def find_header_row(rows):
    """Index of the row with the most recognized headers (the first row when none match)."""

    best_index, best_score = 0, 0
    for index, cells in enumerate(rows[:HEADER_SCAN_ROWS]):
        score = header_score(cells)
        if score > best_score:
            best_index, best_score = index, score
    return best_index


def locate_header(head, chosen=None):
    """The user's header row (0-based) when given and present, else the detected one."""

    if chosen is not None:
        if chosen >= len(head):
            raise ValueError(f"The worksheet has only {len(head)} rows, so row {chosen + 1} cannot be the header.")
        return chosen
    return find_header_row([[cell.strip() for cell in row] for row in head])


def open_workbook(path):
    return openpyxl.load_workbook(path, read_only=True, data_only=True)


def describe_sheets(path):
    """[{name, header_row, score, headers}] from the first rows of every visible worksheet."""

    workbook = open_workbook(path)
    try:
        sheets = []
        for sheet in workbook.worksheets:
            if sheet.sheet_state != "visible":
                continue
            rows = []
            for values in sheet.iter_rows(max_row=PREVIEW_SHEET_ROWS, values_only=True):
                rows.append([cell_text(value).strip() for value in values])
            header_row = find_header_row(rows)
            headers = [cell for cell in rows[header_row] if cell] if rows else []
            sheets.append({"name": sheet.title, "header_row": header_row, "headers": headers,
                           "score": header_score(headers), "empty": not any(any(row) for row in rows)})
        return sheets
    finally:
        workbook.close()


def best_sheet(sheets):
    """The worksheet with the most recognized headers (ties: the first)."""

    usable = [sheet for sheet in sheets if not sheet["empty"]] or sheets
    return max(usable, key=lambda sheet: sheet["score"])["name"] if usable else None


def data_sheets(sheets):
    """Worksheets that look like data tables; more than one means the user must choose."""

    return [sheet for sheet in sheets if sheet["score"] >= MIN_DATA_HEADERS]


def unique_headers(cells):
    """Blank headers get 'Unnamed: n'; repeated ones get '.1', '.2' (as pandas names them)."""

    names, seen = [], {}
    for index, cell in enumerate(cells):
        name = cell or f"Unnamed: {index}"
        count = seen.get(name, 0)
        seen[name] = count + 1
        names.append(name if count == 0 else f"{name}.{count}")
    return names


def read_sheet(path, sheet=None, nrows=None, header_row=None):
    """One worksheet as a text DataFrame; attrs say which sheet and header row were used.

    Fully blank rows are skipped: they are spacing, not records. `header_row`
    (0-based) is the user's choice of header row; without it the row with the
    most recognized column names in the first rows is used.
    """

    workbook = open_workbook(path)
    try:
        if sheet is None:
            sheet = best_sheet(describe_sheets(path))
        if sheet not in workbook.sheetnames:
            raise ValueError(f"The workbook has no worksheet named '{sheet}'.")
        limit = max_cells()
        chosen, header_row = header_row, None
        scan = HEADER_SCAN_ROWS if chosen is None else max(HEADER_SCAN_ROWS, chosen + 1)
        head, records, width = [], [], 0
        preamble = []
        for values in workbook[sheet].iter_rows(values_only=True):
            cells = [cell_text(value) for value in values]
            if header_row is None and len(head) < scan:
                head.append(cells)
                continue
            if header_row is None:
                header_row = locate_header(head, chosen)
                preamble = head[:header_row]
                width = len(head[header_row])
                records = [row for row in head[header_row + 1:] if any(row)]
            if any(cells):
                # Cells beyond the header are never used (rows are cut to the header width
                # below), so they are not stored either; the cell limit then counts what is kept.
                records.append(cells[:width])
            if (len(records) + header_row + 1) * max(width, 1) > limit:
                raise WorkbookTooLarge(
                    f"Worksheet '{sheet}' is larger than this server can process safely (more than "
                    f"{limit:,} cells). Split it into smaller files, or run the analysis locally.")
            if nrows is not None and len(records) >= nrows:
                break
        if header_row is None:  # the sheet is shorter than the scanned rows
            if not head:
                raise ValueError(f"Worksheet '{sheet}' is empty.")
            header_row = locate_header(head, chosen)
            preamble = head[:header_row]
            width = len(head[header_row])
            records = [row for row in head[header_row + 1:] if any(row)]
        columns = unique_headers([cell.strip() for cell in head[header_row]])
        if nrows is not None:
            records = records[:nrows]
        rows = [(row + [""] * width)[:width] for row in records]
        table = pd.DataFrame(rows, columns=columns, dtype=object)
        table.attrs.update(sheet=sheet, header_row=header_row,
                           preamble=" ".join(cell for row in preamble for cell in row if cell).strip())
        return table
    finally:
        workbook.close()
