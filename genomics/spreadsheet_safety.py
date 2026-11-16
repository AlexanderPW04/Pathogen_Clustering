"""Keep downloaded tables from running as spreadsheet formulas (audit P1-2).

Excel, LibreOffice, and Google Sheets treat a cell that starts with ``=``,
``+``, ``-``, or ``@`` as a formula, both in an XLSX workbook and when a CSV
file is opened. An uploaded sample ID such as ``=HYPERLINK(...)`` would then
become a live formula in the downloads. The usual defence (OWASP "CSV
injection") is to put an apostrophe in front of such text, which spreadsheets
show as plain text.

Only text is changed, and only text that could start a formula. Numbers,
including negative numbers written as text ("-0.25"), and a lone "+" or "-"
(for example a gene's strand) are left exactly as they are, so no value
used by an analysis changes. The escaping is applied to the downloads only.
"""

from io import BytesIO

import numpy as np
import pandas as pd

ALWAYS_RISKY = ("=", "@", "\t", "\r")
SIGNS = ("+", "-")
PREFIX = "'"


def is_number(text):
    try:
        float(text)
    except ValueError:
        return False
    return True


def escape_text(value):
    """One value, escaped if it is text that a spreadsheet would read as a formula."""

    if not isinstance(value, str) or not value:
        return value
    first = value[0]
    if first in ALWAYS_RISKY or (first in SIGNS and len(value) > 1 and not is_number(value)):
        return PREFIX + value
    return value


def escape_frame(table):
    """The table with formula-like text escaped; the same object when nothing needs it."""

    escaped = None
    for position, column in enumerate(table.columns):
        values = table.iloc[:, position]
        if not (pd.api.types.is_object_dtype(values) or pd.api.types.is_string_dtype(values)
                or isinstance(values.dtype, pd.CategoricalDtype)):
            continue
        texts = values.astype(object)
        try:
            first = texts.str[:1]  # NaN for values that are not text
        except AttributeError:  # no text in this column
            continue
        candidates = np.flatnonzero(first.isin(ALWAYS_RISKY + SIGNS).to_numpy())
        if not len(candidates):
            continue
        cells = texts.to_numpy(copy=True)
        changed = False
        for index in candidates:
            value = escape_text(cells[index])
            if value is not cells[index]:
                cells[index] = value
                changed = True
        if not changed:
            continue
        if escaped is None:
            escaped = table.copy()
        escaped.isetitem(position, cells)
    renamed = {name: escape_text(name) for name in table.columns if isinstance(name, str)
               and escape_text(name) != name}
    if renamed:
        escaped = (table.copy() if escaped is None else escaped).rename(columns=renamed)
    return table if escaped is None else escaped


def safe_csv(table):
    """CSV bytes of a table, as `table.to_csv(index=False)`, with formula-like text escaped."""

    return escape_frame(table).to_csv(index=False).encode("utf-8")


def escape_csv_bytes(content):
    """Escape formula-like text in CSV bytes written elsewhere; unchanged bytes when nothing needs it."""

    try:
        table = pd.read_csv(BytesIO(content), dtype=str, keep_default_na=False)
    except (ValueError, UnicodeDecodeError, pd.errors.ParserError):
        return content
    escaped = escape_frame(table)
    if escaped is table:
        return content
    ending = "\r\n" if b"\r\n" in content else "\n"
    return escaped.to_csv(index=False, lineterminator=ending).encode("utf-8")
