"""Read pooled-sample collection dates without guessing the day/month order.

A date such as 03/04/2021 can mean 3 April or 4 March. This module never
picks one silently. Each value is read on its own with these rules:

- Year first (2021-04-03, 2021/4/3, optionally followed by a time such as
  "00:00:00" from a spreadsheet date cell) is always year-month-day.
- Month names (3 Apr 2021, April 3, 2021, 03_Apr_2021) are unambiguous.
- Numeric day/month/year values are unambiguous when one of the first two
  numbers is above 12 (13/04/2021 must be day first) or when both are equal.
- Otherwise the value is ambiguous. It is kept as text, flagged, and left out
  of every calculation that needs an exact date, unless the user chose the
  order ("day_first" or "month_first"). A value that contradicts the chosen
  order (04/13/2021 with day first) is invalid.
- A two-digit year (03/04/21) is ambiguous in its century as well, so it is
  never resolved automatically.
- A year or a year and month alone (2021, 2021-04) is a partial date: no day
  is assumed.

`parse_collection_dates` returns the parsed dates and a status for each value,
so the original text and the reason it was or was not used are both kept.
"""

import re
from datetime import datetime

import pandas as pd


AUTOMATIC = "auto"
DAY_FIRST = "day_first"
MONTH_FIRST = "month_first"
DATE_ORDERS = {
    AUTOMATIC: "Unresolved: ambiguous day/month dates are flagged and not used",
    DAY_FIRST: "Day first (DD/MM/YYYY)",
    MONTH_FIRST: "Month first (MM/DD/YYYY)",
}

PARSED = "Parsed"
PARSED_WITH_ORDER = "Parsed with the selected day/month order"
NOT_SUPPLIED = "Not supplied"
AMBIGUOUS = "Ambiguous day/month order (not used)"
TWO_DIGIT_YEAR = "Two-digit year (not used)"
PARTIAL = "Partial date without a day (not used)"
INVALID = "Invalid"

YEAR_FIRST = re.compile(r"^(\d{4})[-/.](\d{1,2})[-/.](\d{1,2})(?:[T ]\d{1,2}:\d{2}(?::\d{2}(?:\.\d+)?)?Z?)?$")
NUMERIC = re.compile(r"^(\d{1,2})[-/.](\d{1,2})[-/.](\d{2}|\d{4})$")
PARTIAL_DATE = re.compile(r"^\d{4}(?:[-/.]\d{1,2})?$")
NAMED_FORMATS = ["%d %b %Y", "%d %B %Y", "%b %d %Y", "%B %d %Y", "%Y %b %d", "%Y %B %d"]


def build_date(year, month, day):
    try:
        return pd.Timestamp(datetime(int(year), int(month), int(day)))
    except ValueError:
        return None


def parse_one(text, date_order=AUTOMATIC, fixed_format=None):
    """Return (Timestamp or None, status) for one collection-date value."""

    text = "" if text is None or (not isinstance(text, str) and pd.isna(text)) else str(text).strip()
    if not text or text.casefold() in {"nan", "nat", "none"}:
        return None, NOT_SUPPLIED
    if fixed_format:
        try:
            return pd.Timestamp(datetime.strptime(text, fixed_format)), PARSED
        except ValueError:
            return None, INVALID
    match = YEAR_FIRST.match(text)
    if match:
        date = build_date(*match.groups())
        return (date, PARSED) if date is not None else (None, INVALID)
    if PARTIAL_DATE.match(text):
        return None, PARTIAL
    match = NUMERIC.match(text)
    if match:
        first, second, year = match.groups()
        first, second = int(first), int(second)
        if len(year) == 2:
            return None, TWO_DIGIT_YEAR
        if date_order == DAY_FIRST:
            date = build_date(year, second, first)
            return (date, PARSED_WITH_ORDER) if date is not None else (None, INVALID)
        if date_order == MONTH_FIRST:
            date = build_date(year, first, second)
            return (date, PARSED_WITH_ORDER) if date is not None else (None, INVALID)
        if first == second or (first > 12) != (second > 12):
            day, month = (first, second) if first > 12 else (second, first)
            date = build_date(year, month, day)
            return (date, PARSED) if date is not None else (None, INVALID)
        if first > 12 and second > 12:
            return None, INVALID
        return None, AMBIGUOUS
    words = " ".join(re.sub(r"[_,/\-.]", " ", text).split())
    for date_format in NAMED_FORMATS:
        try:
            return pd.Timestamp(datetime.strptime(words, date_format)), PARSED
        except ValueError:
            continue
    return None, INVALID


def parse_collection_dates(values, date_order=AUTOMATIC, fixed_format=None):
    """Parse a column of dates; return a DataFrame with Date, Status, and Original.

    The index matches `values`. Unused values have Date = NaT, never a guess.
    """

    if date_order not in DATE_ORDERS:
        raise ValueError("The day/month order must be 'auto', 'day_first', or 'month_first'.")
    values = pd.Series(values)
    originals = values.fillna("").astype(str).str.strip()
    cache = {}
    dates, statuses = [], []
    for text in originals:
        if text not in cache:
            cache[text] = parse_one(text, date_order, fixed_format)
        date, status = cache[text]
        dates.append(date if date is not None else pd.NaT)
        statuses.append(status)
    return pd.DataFrame({"Date": pd.to_datetime(pd.Series(dates, index=values.index, dtype=object)),
                         "Status": statuses, "Original": originals}, index=values.index)


def day_first_evidence(originals):
    """Example values that can only be read day first or only month first."""

    day_first, month_first = [], []
    for text in pd.Series(originals).dropna().astype(str).unique():
        match = NUMERIC.match(text.strip())
        if not match or len(match.group(3)) != 4:
            continue
        first, second = int(match.group(1)), int(match.group(2))
        if first > 12 >= second:
            day_first.append(text)
        elif second > 12 >= first:
            month_first.append(text)
    return day_first[:3], month_first[:3]


def date_warnings(parsed, samples=None):
    """Plain-language warnings for unused dates; `samples` gives the sample of each row."""

    messages = []
    status = parsed["Status"]

    def where(mask):
        examples = list(dict.fromkeys(parsed.loc[mask, "Original"]))[:3]
        sample_text = ""
        if samples is not None:
            sample_text = f" in {pd.Series(samples)[mask.to_numpy()].nunique()} sample(s)"
        return f"{int(mask.sum())} record(s){sample_text} (for example {', '.join(repr(x) for x in examples)})"

    ambiguous = status.eq(AMBIGUOUS)
    if ambiguous.any():
        day_first, month_first = day_first_evidence(parsed["Original"])
        hint = ""
        if day_first and not month_first:
            hint = f" Other dates in the file can only be read day first (for example {day_first[0]!r})."
        elif month_first and not day_first:
            hint = f" Other dates in the file can only be read month first (for example {month_first[0]!r})."
        messages.append(
            f"{where(ambiguous)} have a collection date whose day and month order is ambiguous. "
            "The original text is kept, but these dates are not used for time-based analyses. "
            "Choose the day/month order under Advanced settings to use them." + hint)
    if status.eq(TWO_DIGIT_YEAR).any():
        messages.append(f"{where(status.eq(TWO_DIGIT_YEAR))} use a two-digit year; the century and the "
                        "day/month order are not assumed, so these dates are not used. Use YYYY-MM-DD.")
    if status.eq(PARTIAL).any():
        messages.append(f"{where(status.eq(PARTIAL))} give only a year or a year and month; no day is "
                        "assumed, so these dates are not used for time-based analyses.")
    if status.eq(PARSED_WITH_ORDER).any():
        messages.append(f"{int(status.eq(PARSED_WITH_ORDER).sum())} numeric collection date(s) were read "
                        "with the day/month order you selected.")
    return messages
