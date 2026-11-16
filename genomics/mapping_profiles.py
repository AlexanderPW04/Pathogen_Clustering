"""Save and reuse a confirmed column mapping as a small JSON profile.

Sprint 5 Phase 4A. A profile says *how to read* files of one documented
format: which role each file has, which column supplies each standard field,
and the interpretation choices (scale, numbering, stated reference, depth
definition, base-count minimum, wide layout, header row, worksheet, file-name
sample IDs). It never contains sequencing records, sample values, or IDs.

A profile is applied to a new upload only where it fits: a file must have the
same layout and every column the profile maps. A file that does not fit is
listed with the missing columns and keeps its own detected mapping; nothing
is applied to an incompatible schema silently.
"""

import json
from datetime import date

from genomics.import_options import OPTION_NAMES

PROFILE_FORMAT = "genomic-surveillance-mapping-profile"
PROFILE_VERSION = 1
MAX_PROFILE_BYTES = 64 * 1024
PER_FILE_OPTIONS = [name for name in OPTION_NAMES]
GLOBAL_OPTIONS = ["date_order"]


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
    return profile


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
