"""Recognize two processed-table layouts and rewrite them as long SNV records.

Sprint 5 Phase 4A. The existing loaders read *long* tables: one row per
sample and variant. Two other layouts are common in published supplements:

1. **Wide sample-by-mutation matrices**: one row per sample, one column per
   mutation, frequencies in the cells (for example the PiGx SARS-CoV-2
   ``data_mutation_plot.csv`` tables of Schumann et al. 2022). A column is a
   mutation only when its header states one single-nucleotide change with its
   position (``C10029T``, ``10029C>T``, ``g.10029C>T``, or the PiGx form
   ``ORF1ab:T3255I::C10029T`` whose nucleotide part follows the last ``::``).
   Blank and "NA" cells mean the mutation was *not reported* for that sample:
   they become no record at all, never a 0% frequency.
2. **Per-position base-count tables with a reference base** (for example the
   inStrain-style Table S4 of Crits-Christoph et al. 2021): Sample, Position,
   Ref, and A/C/G/T read counts, but no Alt column. Each non-reference base
   with at least the user-confirmed minimum number of reads becomes one
   observation with *measured* counts; its frequency is calculated here and
   labelled as calculated, not reported.

Both rewrites need the user's confirmation on the review page. Nothing else
is guessed: other columns are carried along unchanged, and tool-specific
columns (for example inStrain's ``con_base``, ``var_base``, or ``var_freq``,
whose ``var_base`` can be the reference base) are never reinterpreted.
"""

import re

import numpy as np
import pandas as pd

from genomics.column_detection import ROLE_FIELDS, squash

WIDE_LAYOUT = "Wide matrix: one row per sample, one column per mutation"
BASE_COUNT_LAYOUT = "Per-position base counts (A/C/G/T) with a reference base"
LONG_LAYOUT = "Long table: one row per sample and variant"

MISSING_TOKENS = {"", "na", "n/a", "nan", "null", "none"}
MINIMUM_MUTATION_COLUMNS = 3

PATTERNS = [
    re.compile(r"^(?P<ref>[ACGT])(?P<pos>[1-9]\d*)(?P<alt>[ACGT])$"),
    re.compile(r"^(?:g\.)?(?P<pos>[1-9]\d*)(?P<ref>[ACGT])>(?P<alt>[ACGT])$"),
]
# Headers that describe a change but not a single-nucleotide one (reported, never parsed).
CHANGE_LIKE = re.compile(r"::|^[ACGTN-]+\d+[ACGTN-]+$|^(?:g\.)?\d+[ACGTN]*(?:>|del|ins|dup)", re.IGNORECASE)

DERIVED_FREQUENCY = ("Calculated by the importer from A/C/G/T counts (alternative reads / A+C+G+T); "
                     "the source reported no frequency")
WIDE_FREQUENCY = "Reported by the source in a wide sample-by-mutation matrix"


def parse_mutation_header(name):
    """{"Ref", "Position", "Alt"} for a single-nucleotide mutation header, else None."""

    text = str(name).strip()
    nucleotide = text.rsplit("::", 1)[-1].strip()
    for pattern in PATTERNS:
        match = pattern.match(nucleotide.upper() if not nucleotide.startswith("g.") else
                              "g." + nucleotide[2:].upper())
        if match:
            return {"Ref": match["ref"], "Position": int(match["pos"]), "Alt": match["alt"]}
    return None


def detect_wide_matrix(columns):
    """Describe a wide mutation matrix, or None when the headers are not one.

    Returns {"mutation_columns": {header: parsed}, "other_columns": [...],
    "not_interpreted": [...]}; "not_interpreted" are change-like headers that
    are not single-nucleotide changes (insertions, deletions, multi-base).
    """

    parsed, others, skipped = {}, [], []
    for column in columns:
        change = parse_mutation_header(column)
        if change is not None:
            parsed[column] = change
        elif CHANGE_LIKE.search(str(column).strip()):
            skipped.append(column)
        else:
            others.append(column)
    if len(parsed) < MINIMUM_MUTATION_COLUMNS or not others:
        return None
    return {"mutation_columns": parsed, "other_columns": others, "not_interpreted": skipped}


def is_missing(values):
    return np.isin(np.char.lower(np.char.strip(values.astype(str))), list(MISSING_TOKENS))


def wide_to_long(table, layout):
    """One long record per reported (non-missing) cell; returns (table, counts).

    The original header is kept in Source_Mutation_Label, the cell text in
    Alt_Frequency (validated later exactly like any reported frequency).
    """

    mutation_columns = list(layout["mutation_columns"])
    others = layout["other_columns"]
    cells = table[mutation_columns].to_numpy(dtype=object)
    missing = is_missing(cells)
    rows, columns = np.nonzero(~missing)
    labels = np.array(mutation_columns, dtype=object)[columns]
    long = table[others].iloc[rows].reset_index(drop=True)
    long["Position"] = [layout["mutation_columns"][label]["Position"] for label in labels]
    long["Ref"] = [layout["mutation_columns"][label]["Ref"] for label in labels]
    long["Alt"] = [layout["mutation_columns"][label]["Alt"] for label in labels]
    long["Alt_Frequency"] = [str(value).strip() for value in cells[rows, columns]]
    long["Source_Mutation_Label"] = labels
    long["Source_Frequency_Provenance"] = WIDE_FREQUENCY
    keys = [(layout["mutation_columns"][label]["Position"], layout["mutation_columns"][label]["Ref"],
             layout["mutation_columns"][label]["Alt"]) for label in mutation_columns]
    without = np.flatnonzero(~(~missing).any(axis=1))
    counts = {"samples": len(table), "mutation_columns": len(mutation_columns),
              "rows_without_records": list(table[others[0]].iloc[without]) if len(without) else [],
              "cells": int(missing.size), "reported": int((~missing).sum()),
              "not_reported": int(missing.sum()),
              "repeated_mutations": len(keys) - len(set(keys)),
              "not_interpreted": len(layout["not_interpreted"])}
    long.attrs.update(table.attrs)
    return long, counts


def base_count_columns(matches):
    """{field: column} when a table has Sample, Position, Ref, and A/C/G/T but no Alt; else None."""

    found = {match["field"]: match["column"] for match in matches if match["status"] == "matched"}
    needed = {"Sample_ID", "Position", "Ref", "A_Count", "C_Count", "G_Count", "T_Count"}
    if needed <= set(found) and "Alt" not in found and "Alt_Frequency" not in found:
        return found
    return None


def derive_from_base_counts(table, columns, minimum_reads):
    """Rewrite a base-count table as SNV records with measured counts; returns (table, counts).

    For each row, every base other than Ref with at least `minimum_reads`
    reads is one record: Alt_Depth = that base's reads, Ref_Depth = the Ref
    base's reads, Read_Depth = A+C+G+T. Rows whose counts are not whole
    numbers are refused rather than repaired.
    """

    counts = {}
    for base in "ACGT":
        raw = table[columns[f"{base}_Count"]].astype(str).str.strip()
        values = pd.to_numeric(raw, errors="coerce")
        bad = values.isna() | values.lt(0) | values.ne(values.round())
        if bad.any():
            first = raw[bad].iloc[0]
            raise ValueError(f"Column '{columns[f'{base}_Count']}' must hold whole read counts; "
                             f"{int(bad.sum())} value(s) are not (for example '{first[:20]}').")
        counts[base] = values.astype(int).to_numpy()
    ref = table[columns["Ref"]].astype(str).str.strip().str.upper().to_numpy()
    depth = counts["A"] + counts["C"] + counts["G"] + counts["T"]
    used = {columns[name] for name in ("Sample_ID", "Position", "Ref", "A_Count", "C_Count", "G_Count",
                                       "T_Count")}
    kept = [column for column in table.columns if column not in used]
    # Columns that could be mistaken for the derived values keep a Source_ name.
    derived = {"Alt", "Alt_Frequency", "Read_Depth", "Alt_Depth", "Ref_Depth"}
    aliases = {alias for spec in ROLE_FIELDS["snv"] if spec["name"] in derived for alias in spec["aliases"]}
    renamed = {column: f"Source_{column}" for column in kept if squash(column) in aliases}
    frames = []
    for base in "ACGT":
        rows = np.flatnonzero((ref != base) & (counts[base] >= minimum_reads) & (counts[base] > 0))
        if not len(rows):
            continue
        part = pd.DataFrame({
            "Sample_ID": table[columns["Sample_ID"]].iloc[rows].to_numpy(),
            "Position": table[columns["Position"]].iloc[rows].to_numpy(),
            "Ref": ref[rows],
            "Alt": base,
            "Alt_Depth": counts[base][rows],
            "Ref_Depth": [counts[r][i] if r in counts else pd.NA for r, i in zip(ref[rows], rows)],
            "Read_Depth": depth[rows],
            "Alt_Frequency": [repr(float(a) / float(d)) for a, d in zip(counts[base][rows], depth[rows])],
            "Source_Frequency_Provenance": DERIVED_FREQUENCY,
            "Source_Base_Count_Minimum": minimum_reads,
            "Source_Row": rows + 1,
        })
        for column in kept:
            part[renamed.get(column, column)] = table[column].iloc[rows].to_numpy()
        frames.append(part)
    if not frames:
        raise ValueError(f"No non-reference base has at least {minimum_reads} read(s), so the base counts "
                         "contain no variant observations.")
    long = pd.concat(frames, ignore_index=True).sort_values(["Source_Row", "Alt"], kind="stable")
    long = long.reset_index(drop=True)
    summary = {"rows": len(table), "records": len(long),
               "rows_without_variant": int(len(table) - long["Source_Row"].nunique()),
               "rows_with_several": int(long["Source_Row"].value_counts().gt(1).sum()),
               "minimum_reads": minimum_reads, "renamed": renamed, "depth_checks": {}}
    for original, name in renamed.items():
        values = pd.to_numeric(long[name], errors="coerce")
        if values.notna().all() and values.eq(long["Read_Depth"]).all() and squash(original) in {
                alias for spec in ROLE_FIELDS["snv"] if spec["name"] == "Read_Depth" for alias in spec["aliases"]}:
            summary["depth_checks"][original] = "equals A+C+G+T on every derived record"
        elif squash(original) in {alias for spec in ROLE_FIELDS["snv"] if spec["name"] == "Read_Depth"
                                  for alias in spec["aliases"]}:
            differs = int((values.ne(long["Read_Depth"]) | values.isna()).sum())
            summary["depth_checks"][original] = f"differs from A+C+G+T on {differs} derived record(s)"
    long.attrs.update(table.attrs)
    return long, summary
