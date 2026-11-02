"""Can a publication document a sample's collection? A simple date check.

A paper cannot report the collection of a specimen that was collected after
it was published. PubMed dates are often partial ("2022 Apr", "2022 Mar-Apr",
"2022 Spring"), so each date is turned into the latest day it can mean, and
a sample's collection date into the earliest day it can mean. Only when the
publication certainly ended before the collection certainly started is the
paper ruled out; anything uncertain is left to the other checks.
"""

import re

import pandas as pd

MONTHS = {name: number for number, name in enumerate(
    ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"], start=1)}


def latest_day(text):
    """The last day a PubMed-style date can mean ('2022 Mar-Apr' -> 2022-04-30), or None."""

    text = str(text or "").strip()
    year = re.search(r"\b(19|20)\d{2}\b", text)
    if not year:
        return None
    rest = text[year.end():].casefold()
    months = [MONTHS[word[:3]] for word in re.findall(r"[a-z]+", rest) if word[:3] in MONTHS]
    if not months:
        return pd.Timestamp(int(year.group()), 12, 31)  # year only, or a season
    month = months[-1]
    day = re.search(r"\b(\d{1,2})\b", rest)
    if day and len(months) == 1 and 1 <= int(day.group(1)) <= 31:
        try:
            return pd.Timestamp(int(year.group()), month, int(day.group(1)))
        except ValueError:
            pass
    return pd.Timestamp(int(year.group()), month, 1) + pd.offsets.MonthEnd(0)


def earliest_day(value):
    """The first day a collection date can mean ('2026' -> 2026-01-01), or None."""

    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return None
    if isinstance(value, pd.Timestamp):
        return value.normalize()
    text = str(value).strip()
    if re.fullmatch(r"\d{4}", text):
        return pd.Timestamp(int(text), 1, 1)
    if re.fullmatch(r"\d{4}-\d{1,2}", text):
        year, month = text.split("-")
        return pd.Timestamp(int(year), int(month), 1)
    parsed = pd.to_datetime(text, errors="coerce")
    return None if pd.isna(parsed) else parsed.normalize()


def chronology(publication_dates, collection_date):
    """(predates, explanation) for one publication and one sample.

    publication_dates holds the print and electronic dates (either may be
    blank); the earliest one counts, because the paper existed from then on.
    predates is True only when the paper was certainly published before the
    sample was collected; None when it cannot be checked.
    """

    ends = [day for day in (latest_day(text) for text in publication_dates) if day is not None]
    published = next((text for text in publication_dates if latest_day(text) == min(ends)), "") if ends else ""
    collected = earliest_day(collection_date)
    if not ends:
        return None, "Not checked: the publication date is unknown."
    if collected is None:
        return None, "Not checked: the sample has no collection date."
    if min(ends) < collected:
        return True, (f"Published {published}, before this sample was collected "
                      f"({collected.date().isoformat()}); it cannot document this collection event.")
    return False, f"Published {published}, not before the sample's collection ({collected.date().isoformat()})."
