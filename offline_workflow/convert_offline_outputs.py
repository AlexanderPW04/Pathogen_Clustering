"""Convert offline raw-read workflow outputs into the website's input tables.

This module is part of the OPTIONAL offline workflow (see README.md in this
folder). It does not run BWA, SAMtools, or iVar; it only reads their text
outputs:

1. `combine_ivar_tables`: several per-sample `ivar variants` TSV files ->
   one TSV with a SAMPLE column (accepted as the pooled SNV table).
2. `pileup_to_base_counts`: `samtools mpileup` text -> A/C/G/T counts per
   reference position (accepted as the per-position table, and required for
   the paper's Yue-Clayton PCoA). Bases below the minimum base quality (20,
   as in Fontenele et al. 2021) are not counted.

Run from the project folder, for example:

    python offline_workflow/convert_offline_outputs.py ivar results/ivar_tsv pooled_snvs.tsv
    python offline_workflow/convert_offline_outputs.py pileup sample1.mpileup S1 counts_S1.csv
"""

import sys
from pathlib import Path

import pandas as pd


BASE_COLUMNS = {"A": "A_Count", "C": "C_Count", "G": "G_Count", "T": "T_Count"}


def combine_ivar_tables(folder):
    """Stack <sample>.tsv iVar outputs, adding SAMPLE from each file name."""

    frames = []
    for path in sorted(Path(folder).glob("*.tsv")):
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        table.insert(0, "SAMPLE", path.stem)
        frames.append(table)
    if not frames:
        raise ValueError(f"No .tsv files were found in {folder}.")
    return pd.concat(frames, ignore_index=True)


def count_pileup_bases(reference_base, bases, qualities, minimum_quality):
    """Count A/C/G/T in one mpileup base string, skipping low-quality bases."""

    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    reference_base = reference_base.upper()
    index = 0
    quality_index = 0
    while index < len(bases):
        symbol = bases[index]
        if symbol == "^":
            index += 2  # read start marker plus its mapping-quality character
            continue
        if symbol == "$":
            index += 1  # read end marker; it has no quality value
            continue
        if symbol in "+-":
            # An indel after this position: +2AG or -3TTT. Skip the digits and bases.
            index += 1
            digits = ""
            while index < len(bases) and bases[index].isdigit():
                digits += bases[index]
                index += 1
            index += int(digits) if digits else 0
            continue

        quality = ord(qualities[quality_index]) - 33 if quality_index < len(qualities) else 0
        quality_index += 1
        index += 1
        if symbol in ".,":
            base = reference_base
        elif symbol.upper() in counts:
            base = symbol.upper()
        else:
            continue  # '*' deletion, '<' '>' reference skips, or 'N'
        if base in counts and quality >= minimum_quality:
            counts[base] += 1
    return counts


def pileup_to_base_counts(path, sample_id, minimum_quality=20):
    """Read `samtools mpileup` output into one row per reference position."""

    rows = []
    with Path(path).open(encoding="utf-8") as pileup:
        for line_number, line in enumerate(pileup, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"Line {line_number} does not look like samtools mpileup output.")
            reference_id, position, reference_base, _, bases, qualities = fields[:6]
            counts = count_pileup_bases(reference_base, bases, qualities, minimum_quality)
            row = {"Sample_ID": sample_id, "Reference_ID": reference_id, "Position": int(position)}
            for base, column in BASE_COLUMNS.items():
                row[column] = counts[base]
            rows.append(row)
    return pd.DataFrame(rows, columns=["Sample_ID", "Reference_ID", "Position"] + list(BASE_COLUMNS.values()))


def main(arguments):
    if len(arguments) == 3 and arguments[0] == "ivar":
        combine_ivar_tables(arguments[1]).to_csv(arguments[2], sep="\t", index=False)
        return 0
    if len(arguments) == 4 and arguments[0] == "pileup":
        pileup_to_base_counts(arguments[1], arguments[2]).to_csv(arguments[3], index=False)
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
