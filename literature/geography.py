"""Compare proposed locations and dates with the original sample metadata.

The checks decide whether a proposal is more informative than the original
value and whether it is consistent with the original country and state. A
proposal outside the original state or country is a contradiction for a
person to review, never a silent replacement.
"""

import re

import pandas as pd

from geocode_locations import US_NAMES, US_STATES, is_state_only_location
from validate_inputs import is_missing_value

STATE_ABBREVIATIONS = {
    "al": "alabama", "ak": "alaska", "az": "arizona", "ar": "arkansas", "ca": "california",
    "co": "colorado", "ct": "connecticut", "de": "delaware", "fl": "florida", "ga": "georgia",
    "hi": "hawaii", "id": "idaho", "il": "illinois", "in": "indiana", "ia": "iowa", "ks": "kansas",
    "ky": "kentucky", "la": "louisiana", "me": "maine", "md": "maryland", "ma": "massachusetts",
    "mi": "michigan", "mn": "minnesota", "ms": "mississippi", "mo": "missouri", "mt": "montana",
    "ne": "nebraska", "nv": "nevada", "nh": "new hampshire", "nj": "new jersey", "nm": "new mexico",
    "ny": "new york", "nc": "north carolina", "nd": "north dakota", "oh": "ohio", "ok": "oklahoma",
    "or": "oregon", "pa": "pennsylvania", "ri": "rhode island", "sc": "south carolina",
    "sd": "south dakota", "tn": "tennessee", "tx": "texas", "ut": "utah", "vt": "vermont",
    "va": "virginia", "wa": "washington", "wv": "west virginia", "wi": "wisconsin", "wy": "wyoming",
    "dc": "district of columbia", "pr": "puerto rico",
}
LEVEL_RANK = {"Not geographic": 0, "Country": 1, "State or province": 2, "County": 3,
              "City or town": 4, "Facility or site": 5, "Coordinates": 6}
PRECISE_COORDINATES = 6


def clean(text):
    return " ".join(str(text).split()).casefold()


def parse_location(text):
    """{'country', 'state', 'level', 'name'} for labels like 'USA: AZ, Tempe' or 'Tempe, Arizona'."""

    if text is None or (not isinstance(text, str) and pd.isna(text)) or is_missing_value(text):
        return {"country": None, "state": None, "level": 0, "name": ""}
    parts = [part.strip() for part in re.split(r"[:,;]", clean(text)) if part.strip()]
    country = state = None
    other = []
    for part in parts:
        if part in US_NAMES:
            country = "usa"
        elif part in US_STATES:
            state = part
        elif part in STATE_ABBREVIATIONS and state is None and len(parts) > 1:
            state = STATE_ABBREVIATIONS[part]
        else:
            other.append(part)
    if state and not country:
        country = "usa"
    if other:
        level = 3 if all("county" in part for part in other) else 4
    elif state:
        level = 2
    else:
        level = 1 if country else 4  # an unrecognized name is treated as a named place
    return {"country": country, "state": state, "level": level, "name": other[0] if other else ""}


def original_geography(row):
    """Level and context of a sample's original geography (6 = reliable coordinates)."""

    parsed = parse_location(row.get("Location_Name"))
    precision = str(row.get("Geographic_Precision", ""))
    has_coordinates = pd.notna(row.get("Latitude")) and pd.notna(row.get("Longitude"))
    coarse = "too coarse" in precision or is_state_only_location(row.get("Location_Name"))
    if has_coordinates and not coarse:
        parsed["level"] = PRECISE_COORDINATES
    parsed["coarse"] = coarse
    return parsed


def compare_location(row, value, level):
    """Return (comparison, note) for a proposed location."""

    original = original_geography(row)
    proposed = parse_location(value)
    if original["level"] == PRECISE_COORDINATES:
        return "Original coordinates kept", "The sample already has coordinates; they are never overwritten."
    if original["state"] and proposed["state"] and original["state"] != proposed["state"]:
        return "Conflicts with original", (f"Proposed state '{proposed['state'].title()}' differs from the "
                                           f"original '{original['state'].title()}'.")
    if original["country"] and proposed["country"] and original["country"] != proposed["country"]:
        return "Conflicts with original", "The proposed country differs from the original."
    note = ""
    if original["state"] and not proposed["state"]:
        note = "The proposal does not name a state; GeoNames checks the state when it is applied."
    rank = LEVEL_RANK.get(level, 0)
    if original["level"] == 0:
        return "Fills missing value", note
    if rank > original["level"]:
        return "More specific than original", note
    if rank == original["level"] and proposed["name"] == original["name"]:
        return "Same as original", note
    return "Not more specific than original", note


def compare_coordinates(row, latitude, longitude):
    original = original_geography(row)
    if original["level"] == PRECISE_COORDINATES:
        return "Original coordinates kept", "The sample already has coordinates; they are never overwritten."
    if not (-90 <= latitude <= 90 and -180 <= longitude <= 180):
        return "Conflicts with original", "Coordinates are outside the valid range."
    return "Fills missing value", "Coordinates are shown as reported in the publication."


def compare_date(row, start, end=None):
    """Compare a proposed date (or range start..end) with the original collection date."""

    original = pd.to_datetime(row.get("Collection_Date"), errors="coerce")
    if pd.isna(original):
        return "Fills missing value", "" if end is None else "A range is kept as a range, not an exact date."
    original = original.normalize()
    if end is None:
        if original == start:
            return "Same as original", ""
        return "Conflicts with original", f"The original date is {original.date()}."
    if start <= original <= end:
        return "Same as original", "The original date lies inside the documented range."
    return "Conflicts with original", f"The original date {original.date()} is outside the range."
