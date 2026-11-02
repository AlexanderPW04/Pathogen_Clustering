"""Combine a metadata table's geographic columns into one location per row.

NCBI Virus tables describe one place at several levels: Geo_Location
("USA: SD"), Country ("USA"), and a column named "USA" holding the US state
("SD"). These are levels of the same hierarchy, not competing candidates:

- The specific location column (Location_Name, e.g. from Geo_Location) is
  used as given and never overwritten.
- A row without one gets "Country: State" or "Country" from the broader
  columns, labeled with that level and its source.
- Country and state columns are kept unchanged beside it.
- A row whose columns disagree (another country or state) is flagged for
  review and treated as too coarse for kilometer screening.

The level ("country", "state", or "conflict") lets the Sprint 2 pipeline label
precision correctly even for non-US places that the text alone cannot reveal.
"""

import pandas as pd

from geocode_locations import US_NAMES
from literature.geography import STATE_ABBREVIATIONS
from validate_inputs import is_missing_value

COUNTRY = "country"
STATE = "state"
CONFLICT = "conflict"
COARSE_LEVELS = {COUNTRY, STATE, CONFLICT}
# Headers (compared without case, spaces, underscores, or hyphens). NCBI Virus
# exports put the US state code in a column literally named "USA".
COUNTRY_HEADERS = {"country", "countryname", "nation"}
STATE_HEADERS = {"state", "province", "stateprovince", "admin1", "division", "usa"}


def squash(name):
    return "".join(character for character in str(name).casefold() if character.isalnum())


def find_column(columns, headers):
    matches = [column for column in columns if squash(column) in headers]
    return matches[0] if len(matches) == 1 else None


def text(value):
    if value is None or (not isinstance(value, str) and pd.isna(value)) or is_missing_value(value):
        return ""
    return " ".join(str(value).split())


def country_key(value):
    key = text(value).casefold()
    return "usa" if key in US_NAMES else key


def state_key(value):
    key = text(value).casefold()
    return STATE_ABBREVIATIONS.get(key, key)


def row_location(location, country, state, location_header):
    """(location, level, source, conflict note) for one row."""

    location, country, state = text(location), text(country), text(state)
    if not location:
        if country and state:
            return f"{country}: {state}", STATE, "Built from country and state columns", ""
        if country:
            return country, COUNTRY, "Built from the country column", ""
        if state:
            return state, STATE, "Built from the state column", ""
        return "", "", "", ""
    head, separator, rest = location.partition(":")
    area = rest.split(",")[0].strip() if separator else ""
    problems = []
    if country and separator and country_key(head) != country_key(country):
        problems.append(f"{location_header} '{location}' names a different country than '{country}'")
    if state and area and state_key(area) != state_key(state):
        problems.append(f"{location_header} '{location}' names a different state than '{state}'")
    source = f"{location_header} column"
    if problems:
        return location, CONFLICT, source, "; ".join(problems)
    if country and country_key(location) == country_key(country):
        return location, COUNTRY, source, ""
    if separator and area and "," not in rest and state and state_key(area) == state_key(state):
        return location, STATE, source, ""
    return location, "", source, ""


def location_hierarchy(table, location_column="Location_Name", location_header=None):
    """Per-row locations from a table; returns (DataFrame, summary) without changing `table`.

    The DataFrame has Location_Name, Location_Level, Location_Name_Source,
    and Location_Conflict, aligned with the table's index.
    """

    country_column = find_column(table.columns, COUNTRY_HEADERS)
    state_column = find_column(table.columns, STATE_HEADERS)
    header = location_header or location_column
    rows = []
    for index in table.index:
        location = table.at[index, location_column] if location_column in table else None
        country = table.at[index, country_column] if country_column else None
        state = table.at[index, state_column] if state_column else None
        rows.append(row_location(location, country, state, header))
    result = pd.DataFrame(rows, index=table.index,
                          columns=["Location_Name", "Location_Level", "Location_Name_Source", "Location_Conflict"])
    summary = {"country_column": country_column, "state_column": state_column,
               "filled": int(result["Location_Name_Source"].str.startswith("Built").sum()),
               "country_level": int((result["Location_Level"] == COUNTRY).sum()),
               "state_level": int((result["Location_Level"] == STATE).sum()),
               "conflicts": int((result["Location_Level"] == CONFLICT).sum())}
    return result, summary


def apply_location_hierarchy(table):
    """Fill missing locations from broader columns and record each row's level and source.

    Explicit location values are never changed; country and state columns stay as they are.
    """

    if not ({"Location_Name"} & set(table.columns)) and not (
            find_column(table.columns, COUNTRY_HEADERS) or find_column(table.columns, STATE_HEADERS)):
        return table
    result, _ = location_hierarchy(table)
    table = table.copy()
    original = table["Location_Name"] if "Location_Name" in table else pd.Series("", index=table.index)
    blank = original.map(lambda value: text(value) == "")
    table["Location_Name"] = original.astype(object).where(~blank, result["Location_Name"].replace("", pd.NA))
    table["Location_Level"] = result["Location_Level"].replace("", pd.NA)
    table["Location_Name_Source"] = result["Location_Name_Source"].replace("", pd.NA)
    return table
