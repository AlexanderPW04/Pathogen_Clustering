"""Deterministic checks of model output against the retrieved article text.

A model's quotation is accepted only if it is found in the cited passage;
the proposed value must appear in that quotation; the sample must be named in
it for sample-specific evidence; and the result is compared with the original
metadata. The model's own evidence label can only lower the final status.
"""

import re
import unicodedata

import pandas as pd

from literature.geography import compare_coordinates, compare_date, compare_location
from literature.identifiers import mentions

DIRECT = "Directly stated"
TABLE = "Supported by an explicit sample-to-site table"
INFERRED = "Inferred from broader context"
AMBIGUOUS = "Ambiguous"
CONTRADICTORY = "Contradictory"
UNSUPPORTED = "Unsupported"
NOT_COLLECTION = "Not a collection location"
STUDY_LEVEL = "Study-level context"
STRENGTH = {DIRECT: 4, TABLE: 4, INFERRED: 2, AMBIGUOUS: 1}
MODEL_LABELS = {"Directly stated": DIRECT, "Sample-to-site table": TABLE,
                "Inferred from broader context": INFERRED, "Ambiguous": AMBIGUOUS}
GEOGRAPHIC_FIELDS = {"Location", "Coordinates", "Site_Identifier"}
CHANGES_ANALYSIS = {"Fills missing value", "More specific than original"}
MONTHS = {name: number for number, name in enumerate(
    ["january", "february", "march", "april", "may", "june", "july", "august", "september",
     "october", "november", "december"], start=1)}
MONTH = r"(january|february|march|april|may|june|july|august|september|october|november|december)"


def normalized(text):
    """Casefolded text with unified quotes, dashes, and spaces, plus an index map."""

    characters, positions = [], []
    previous_space = True
    for index, character in enumerate(unicodedata.normalize("NFKC", text)):
        character = {"‘": "'", "’": "'", "“": '"', "”": '"', "–": "-",
                     "—": "-", " ": " "}.get(character, character).casefold()
        if character.isspace():
            if previous_space:
                continue
            character = " "
        previous_space = character == " "
        characters.append(character)
        positions.append(index)
    return "".join(characters).strip(), positions


def locate_quote(quote, text):
    """(start, end) character offsets of the quote in the passage, or None."""

    needle, _ = normalized(quote)
    haystack, positions = normalized(text)
    if len(needle) < 10:
        return None
    start = haystack.find(needle)
    if start < 0:
        return None
    return positions[start], positions[start + len(needle) - 1] + 1


def dates_in(text):
    """Every full calendar date written in the text (ISO, '5 March 2026', 'March 5, 2026', '1 and 14 March 2026')."""

    found = []
    lowered = text.casefold()
    for year, month, day in re.findall(r"\b(\d{4})-(\d{2})-(\d{2})\b", lowered):
        found.append((int(year), int(month), int(day)))
    pattern = r"\b(\d{1,2})(?:\s*(?:-|and|to|–)\s*(\d{1,2}))?\s+" + MONTH + r",?\s+(\d{4})"
    for first, second, month, year in re.findall(pattern, lowered):
        for day in [first, second]:
            if day:
                found.append((int(year), MONTHS[month], int(day)))
    pattern = MONTH + r"\s+(\d{1,2})(?:\s*(?:-|and|to|–)\s*(\d{1,2}))?,?\s+(\d{4})"
    for month, first, second, year in re.findall(pattern, lowered):
        for day in [first, second]:
            if day:
                found.append((int(year), MONTHS[month], int(day)))
    # '1 March and 20 April 2026': the first date shares the year written after the second.
    pattern = r"\b(\d{1,2})\s+" + MONTH + r"\s*(?:-|and|to|–|until|through)\s*(\d{1,2})\s+" + MONTH + r",?\s+(\d{4})"
    for first, month, second, second_month, year in re.findall(pattern, lowered):
        found.append((int(year), MONTHS[month], int(first)))
        found.append((int(year), MONTHS[second_month], int(second)))
    dates = set()
    for year, month, day in found:
        try:
            dates.add(pd.Timestamp(year=year, month=month, day=day))
        except ValueError:
            continue
    return dates


def parse_iso_date(text):
    try:
        return pd.Timestamp(str(text).strip()).normalize()
    except (ValueError, TypeError):
        return None


def coordinates_from(value):
    numbers = re.findall(r"-?\d+(?:\.\d+)?", str(value))
    return (float(numbers[0]), float(numbers[1])) if len(numbers) == 2 else None


def value_in_quote(fact, quote):
    """Check that the proposed value is written in the quotation; return (ok, detail)."""

    field, value = fact["Field"], str(fact["Value"]).strip()
    text = quote.casefold()
    if field in {"Location", "Site_Identifier"}:
        first = re.split(r"[,;:]", value)[0].strip().casefold()
        return (bool(first) and first in text), f"'{first}' must appear in the quotation"
    if field == "Coordinates":
        pair = coordinates_from(value)
        if pair is None:
            return False, "coordinates must be two numbers"
        written = [float(number) for number in re.findall(r"-?\d+\.\d+", quote)]
        return all(any(abs(number - item) < 1e-9 for number in written) for item in pair), \
            "both coordinates must be written in the quotation"
    if field == "Collection_Date":
        day = parse_iso_date(value)
        return (day is not None and day in dates_in(quote)), "the date must be written in the quotation"
    if field == "Collection_Date_Range":
        ends = [parse_iso_date(part) for part in value.split("/")]
        ok = len(ends) == 2 and all(ends) and ends[0] <= ends[1] and set(ends) <= dates_in(quote)
        return ok, "both ends of the range must be written in the quotation"
    return True, ""


def verify_fact(fact, passages, request_samples, all_aliases):
    """Return the fact with deterministic Evidence_Status, Approvable, and notes.

    passages: {Passage_ID: passage} sent in the request.
    request_samples: {Sample_ID: {"Aliases": [...], "Row": metadata row}}.
    all_aliases: {Sample_ID: aliases} for every sample in the dataset.
    """

    checked = {**fact, "Evidence_Start": None, "Evidence_End": None, "Metadata_Comparison": "Not compared",
               "Approvable": False, "Check_Notes": []}
    notes = checked["Check_Notes"]

    def finish(status, note=None):
        if note:
            notes.append(note)
        checked["Evidence_Status"] = status
        checked["Check_Notes"] = " ".join(notes)
        return checked

    sample_id = str(fact["Sample_ID"])
    if sample_id != "STUDY" and sample_id not in request_samples:
        return finish(UNSUPPORTED, "The sample was not part of this request.")
    passage = passages.get(fact["Passage_ID"])
    if passage is None:
        return finish(UNSUPPORTED, "The cited passage was not among those sent to the model.")
    span = locate_quote(fact["Evidence_Quote"], passage["Text"])
    if span is None:
        return finish(UNSUPPORTED, "The quotation was not found in the cited passage.")
    checked["Evidence_Start"], checked["Evidence_End"] = span
    checked["Evidence_Section"] = passage["Section"]
    ok, detail = value_in_quote(fact, fact["Evidence_Quote"])
    if not ok:
        return finish(UNSUPPORTED, f"The proposed value is not in the quotation ({detail}).")
    if sample_id == "STUDY":
        return finish(STUDY_LEVEL, "Study-level fact; it is not assigned to any sample.")
    if fact["Field"] in GEOGRAPHIC_FIELDS and fact["Location_Role"] != "Sample collection site":
        return finish(NOT_COLLECTION, f"The text describes a {fact['Location_Role'].lower()}, "
                                      "not where the sample was collected.")

    quote = fact["Evidence_Quote"]
    names_sample = any(mentions(quote, alias) for alias in request_samples[sample_id]["Aliases"])
    others = [other for other, aliases in all_aliases.items()
              if other != sample_id and any(mentions(quote, alias) for alias in aliases)]
    if names_sample and passage["Section"] == "TABLE":
        status = TABLE
    elif names_sample:
        status = DIRECT
    else:
        status = INFERRED
        notes.append("The quotation does not name this sample.")
    model_status = MODEL_LABELS.get(fact["Evidence_Type"], AMBIGUOUS)
    if STRENGTH[model_status] < STRENGTH[status]:
        status = model_status
        notes.append(f"The model labeled the evidence '{fact['Evidence_Type']}'.")
    if fact.get("Ambiguity", "").strip() and STRENGTH[status] > STRENGTH[AMBIGUOUS]:
        status = AMBIGUOUS
        notes.append(f"Model note: {fact['Ambiguity'].strip()}")
    if others and status != TABLE and STRENGTH[status] > STRENGTH[AMBIGUOUS]:
        status = AMBIGUOUS
        notes.append(f"The same quotation also names {', '.join(sorted(others))}; the sample-to-site "
                     "mapping needs a person to check it.")

    comparison, note = compare_with_original(fact, request_samples[sample_id]["Row"])
    checked["Metadata_Comparison"] = comparison
    if note:
        notes.append(note)
    if comparison == "Conflicts with original" and status in {DIRECT, TABLE}:
        status = CONTRADICTORY
    checked["Approvable"] = status in {DIRECT, TABLE} and comparison in CHANGES_ANALYSIS
    if status in {DIRECT, TABLE} and comparison not in CHANGES_ANALYSIS and comparison != "Conflicts with original":
        notes.append("Nothing to change: the proposal does not add information to the original metadata.")
    return finish(status)


def compare_with_original(fact, row):
    field, value = fact["Field"], str(fact["Value"])
    if field == "Location":
        return compare_location(row, value, fact["Geographic_Level"])
    if field == "Coordinates":
        pair = coordinates_from(value)
        return compare_coordinates(row, *pair)
    if field == "Collection_Date":
        return compare_date(row, parse_iso_date(value))
    if field == "Collection_Date_Range":
        start, end = (parse_iso_date(part) for part in value.split("/"))
        return compare_date(row, start, end)
    return "Fills missing value", "Descriptive field; it does not change any calculation."

