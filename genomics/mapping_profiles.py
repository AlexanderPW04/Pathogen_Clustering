"""Save and reuse a confirmed column mapping as a small JSON profile.

Sprint 5 Phase 4A. A profile says *how to read* files of one documented
format: which role each file has, which column supplies each standard field,
and the interpretation choices (scale, numbering, stated reference, depth
definition, base-count minimum, wide layout, header row, worksheet, file-name
sample IDs). It never contains sequencing records, sample values, or IDs.

A profile is applied to a new upload only where it fits: a file must have the
same layout, worksheet, and every column the profile maps. A file that does
not fit is listed with the missing columns and keeps its own detected mapping;
nothing is applied to an incompatible schema silently.

An uploaded profile is untrusted input (final polish 4A). Every key and value
is checked before anything is used: only known keys, roles, standard field
names, and option values that the review form itself would accept (the same
`read_file_options` rules: header row 1-50, position numbering 0 or 1,
fraction or percent scale, ...). A profile that breaks a rule is refused as a
whole, with the reason; it is never partly applied.
"""

import json
import re
from datetime import date

from genomics.collection_dates import DATE_ORDERS
from genomics.column_detection import ROLE_FIELDS
from genomics.import_options import OPTION_NAMES, read_file_options

PROFILE_FORMAT = "genomic-surveillance-mapping-profile"
PROFILE_VERSION = 1
MAX_PROFILE_BYTES = 64 * 1024
PER_FILE_OPTIONS = [name for name in OPTION_NAMES]
GLOBAL_OPTIONS = ["date_order"]
TOP_LEVEL_KEYS = {"profile_format", "profile_version", "created", "description", "options", "files"}
FILE_KEYS = {"role", "kind", "format", "layout", "sheet", "fields", "options"}
PROFILE_ROLES = {"fasta", "metadata", "reference", "snv", "pooled_metadata", "positions", "masked_sites",
                 "signatures", "annotation"}  # the review page's roles (except "ignore", never saved)
MAX_PROFILE_FILES = 20  # the most files one upload may hold
MAX_TEXT = 200


def build_profile(review, options):
    """A JSON-serializable profile from the current review (mappings and choices only)."""

    files = []
    for entry in review["files"]:
        info, role = entry["info"], entry["role"]
        if info["refused"] or role in {None, "ignore"}:
            continue
        fields = {match["field"]: match["column"] for match in entry.get("matches", [])
                  if match["status"] == "matched"}
        file_options = {key: value for key, value in (entry.get("options") or {}).items() if key in PER_FILE_OPTIONS}
        files.append({
            "role": role,
            "kind": info["kind"],
            "format": info.get("format"),
            "layout": info.get("layout"),
            "sheet": info.get("sheet"),
            "fields": fields,
            "options": file_options,
        })
    return {
        "profile_format": PROFILE_FORMAT,
        "profile_version": PROFILE_VERSION,
        "created": date.today().isoformat(),
        "description": "How to read these files (column mapping and confirmed choices). Contains no records.",
        "options": {key: options.get(key) for key in GLOBAL_OPTIONS if options.get(key)},
        "files": files,
    }


def profile_json(review, options):
    return json.dumps(build_profile(review, options), indent=2).encode("utf-8")


def read_profile(content):
    """Parse and validate a profile; raises ValueError with a plain message."""

    if len(content) > MAX_PROFILE_BYTES:
        raise ValueError("The mapping profile is larger than 64 KB, so it is not a profile made by this application.")
    try:
        profile = json.loads(content.decode("utf-8-sig"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        raise ValueError("The mapping profile is not valid JSON.") from None
    if not isinstance(profile, dict) or profile.get("profile_format") != PROFILE_FORMAT:
        raise ValueError("This JSON file is not a mapping profile saved by this application.")
    if profile.get("profile_version") != PROFILE_VERSION:
        raise ValueError(f"Mapping profile version {profile.get('profile_version')} is not supported.")
    if not isinstance(profile.get("files"), list):
        raise ValueError("The mapping profile lists no files.")
    validate_profile(profile)
    return profile


def short_text(value, where):
    if value is None:
        return None
    if not isinstance(value, str) or len(value) > MAX_TEXT or any(ord(ch) < 32 for ch in value):
        raise ValueError(f"The mapping profile's {where} must be plain text of at most {MAX_TEXT} characters.")
    return value


def validate_options(options, where):
    """Per-file options must be exactly what the review form accepts (no value is changed or dropped)."""

    if not isinstance(options, dict):
        raise ValueError(f"The mapping profile's {where} options must be a JSON object.")
    unknown = sorted(set(options) - set(PER_FILE_OPTIONS))
    if unknown:
        raise ValueError(f"The mapping profile's {where} has unknown option(s): {', '.join(map(str, unknown))}.")
    form = {}
    for key, value in options.items():
        if isinstance(value, bool):
            form[f"opt_0_{key}"] = "on" if value else ""
        elif isinstance(value, (int, str)) and not isinstance(value, float):
            form[f"opt_0_{key}"] = str(value)
        else:
            raise ValueError(f"The mapping profile's {where} option '{key}' has an invalid value.")
    accepted = read_file_options(form, 0)
    for key, value in options.items():
        if value is False or value == "":
            continue  # an option that is off
        same = (accepted.get(key) is True) if isinstance(value, bool) else (
            key in accepted and str(accepted[key]) == str(value).strip())
        if not same:
            rule = {"header_row": "a whole number from 1 to 50",
                    "base_count_minimum": "a whole number from 1 to 1000",
                    "position_base": "0 or 1",
                    "frequency_scale": "'fraction' or 'percent'"}.get(key, "a value the review page offers")
            raise ValueError(f"The mapping profile's {where} option '{key}' is {value!r}; it must be {rule}.")


def validate_profile(profile):
    """Refuse a profile with unknown keys, roles, fields, or option values (see the module notes)."""

    unknown = sorted(set(profile) - TOP_LEVEL_KEYS)
    if unknown:
        raise ValueError("The mapping profile has unknown field(s): " + ", ".join(map(str, unknown)) + ".")
    short_text(profile.get("description"), "description")
    created = profile.get("created")
    if created is not None and not (isinstance(created, str) and re.fullmatch(r"\d{4}-\d{2}-\d{2}", created)):
        raise ValueError("The mapping profile's creation date must look like 2026-11-16.")
    options = profile.get("options") or {}
    if not isinstance(options, dict) or set(options) - set(GLOBAL_OPTIONS):
        raise ValueError("The mapping profile's general options may contain only 'date_order'.")
    if "date_order" in options and options["date_order"] not in DATE_ORDERS:
        raise ValueError("The mapping profile's date order must be one of: " + ", ".join(DATE_ORDERS) + ".")
    files = profile["files"]
    if not files or len(files) > MAX_PROFILE_FILES:
        raise ValueError(f"A mapping profile describes 1 to {MAX_PROFILE_FILES} files.")
    for number, item in enumerate(files, start=1):
        where = f"entry {number}"
        if not isinstance(item, dict):
            raise ValueError(f"The mapping profile's {where} must be a JSON object.")
        unknown = sorted(set(item) - FILE_KEYS)
        if unknown:
            raise ValueError(f"The mapping profile's {where} has unknown field(s): {', '.join(map(str, unknown))}.")
        role = item.get("role")
        if role not in PROFILE_ROLES:
            raise ValueError(f"The mapping profile's {where} has an unknown role {str(role)[:40]!r}.")
        for key in ["kind", "format", "layout", "sheet"]:
            short_text(item.get(key), f"{where} {key}")
        fields = item.get("fields") or {}
        if not isinstance(fields, dict) or len(fields) > 60:
            raise ValueError(f"The mapping profile's {where} columns must be a JSON object of at most 60 fields.")
        allowed = {spec["name"] for spec in ROLE_FIELDS.get(role, [])}
        for field_name, column in fields.items():
            if field_name not in allowed:
                raise ValueError(f"The mapping profile's {where} maps an unknown field {str(field_name)[:40]!r} for "
                                 f"the role '{role}'.")
            short_text(column, f"{where} column for {field_name}")
            if not column:
                raise ValueError(f"The mapping profile's {where} maps {field_name} to an empty column name.")
        validate_options(item.get("options") or {}, where)


def fits(item, info):
    """(True, "") when the profile entry fits this inspected file, else (False, reason)."""

    if info["refused"]:
        return False, "the file was not accepted"
    if item.get("kind") and info["kind"] != item["kind"] and not (
            {info["kind"], item["kind"]} <= {"Table", "Excel workbook"}):
        return False, f"file kind differs ({info['kind']} instead of {item['kind']})"
    if item.get("role") == "reference" and info.get("records") != 1:
        return False, "a reference FASTA must hold exactly one sequence"
    if item.get("layout") and info.get("layout") != item["layout"]:
        return False, f"layout differs ({info.get('layout')} instead of {item['layout']})"
    if item.get("format") and info.get("format") != item["format"]:
        return False, f"format differs ({info.get('format')} instead of {item['format']})"
    if item.get("sheet") and info.get("sheets") and item["sheet"] not in info["sheets"]:
        return False, f"the workbook has no worksheet '{item['sheet']}'"
    columns = set(info.get("columns", []))
    missing = [column for column in item.get("fields", {}).values() if column not in columns]
    if missing:
        return False, "missing column(s) " + ", ".join(f"'{name}'" for name in missing)
    return True, ""


def apply_profile(profile, inspections):
    """Return (roles, choices, file_options, global_options, report) for files the profile fits.

    Each profile entry is used for at most one file per role order, except
    multi-file roles (several SNV files may share one entry).
    """

    roles, choices, file_options, report = {}, {}, {}, []
    used = set()
    for number, item in enumerate(profile["files"]):
        applied = []
        reasons = []
        for index, info in enumerate(inspections):
            if index in used:
                continue
            ok, reason = fits(item, info)
            if not ok:
                reasons.append(f"{info['name']}: {reason}")
                continue
            if item.get("role") == "snv" or not applied:
                applied.append(index)
                if item.get("role") != "snv":
                    break
        for index in applied:
            used.add(index)
            roles[index] = item.get("role")
            choices[index] = {field: column for field, column in item.get("fields", {}).items()}
            file_options[index] = {key: value for key, value in item.get("options", {}).items()
                                   if key in PER_FILE_OPTIONS}
            report.append(f"Profile entry {number + 1} ({item.get('role')}) applied to "
                          f"{inspections[index]['name']}.")
        if not applied:
            report.append(f"Profile entry {number + 1} ({item.get('role')}) was not applied: "
                          + ("; ".join(reasons[:4]) if reasons else "no remaining file") + ".")
    global_options = {key: value for key, value in (profile.get("options") or {}).items() if key in GLOBAL_OPTIONS}
    return roles, choices, file_options, global_options, report
