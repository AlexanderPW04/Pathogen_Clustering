"""Align comparable sequences with MAFFT, then calculate genomic distances."""

import argparse
import sys
from pathlib import Path

import pandas as pd

from align_sequences import align_samples
from validate_inputs import (
    FASTA_EXTENSIONS,
    METADATA_EXTENSIONS,
    load_and_validate_inputs,
)


RESULT_COLUMNS = [
    "Sample_ID_1",
    "Sample_ID_2",
    "Pathogen",
    "Gene_Segment",
    "Differences",
    "Comparable_Positions",
    "P_Distance",
    "Percent_Identity",
]


def compare_sequences(sequence_1, sequence_2):
    """Compare two already-aligned sequences position by position."""

    if len(sequence_1) != len(sequence_2):
        raise ValueError("sequences have different lengths")

    differences = 0
    comparable_positions = 0

    for base_1, base_2 in zip(sequence_1, sequence_2):
        # Unknown bases (N or IUPAC ambiguity codes) and alignment gaps do not
        # count. For sequences containing only A, C, G, T, N, and -, this is
        # the same rule as the original Sprint 2 comparison.
        if base_1 not in "ACGT" or base_2 not in "ACGT":
            continue

        comparable_positions += 1

        if base_1 != base_2:
            differences += 1

    if comparable_positions == 0:
        p_distance = None
        percent_identity = None
    else:
        p_distance = differences / comparable_positions
        matching_positions = comparable_positions - differences
        percent_identity = matching_positions / comparable_positions * 100

    return {
        "Differences": differences,
        "Comparable_Positions": comparable_positions,
        "P_Distance": p_distance,
        "Percent_Identity": percent_identity,
    }


def compare_samples(data):
    """Compare every pair with the same Pathogen and Gene_Segment."""

    result_rows = []
    errors = []

    groups = data.groupby(
        ["Pathogen", "Gene_Segment"],
        sort=True,
        dropna=False,
    )

    for group_names, group in groups:
        pathogen, gene_segment = group_names
        records = group.to_dict("records")

        # Use two basic loops to create every unique pair once.
        for first_index in range(len(records)):
            for second_index in range(first_index + 1, len(records)):
                sample_1 = records[first_index]
                sample_2 = records[second_index]
                sequence_1 = str(sample_1["Sequence"])
                sequence_2 = str(sample_2["Sequence"])

                try:
                    measurements = compare_sequences(sequence_1, sequence_2)
                except ValueError:
                    errors.append(
                        f"Cannot compare {sample_1['Sample_ID']} and "
                        f"{sample_2['Sample_ID']} ({pathogen}, {gene_segment}): "
                        f"sequence lengths are {len(sequence_1)} and "
                        f"{len(sequence_2)}."
                    )
                    continue

                result_row = {
                    "Sample_ID_1": sample_1["Sample_ID"],
                    "Sample_ID_2": sample_2["Sample_ID"],
                    "Pathogen": pathogen,
                    "Gene_Segment": gene_segment,
                }
                result_row.update(measurements)
                result_rows.append(result_row)

    results = pd.DataFrame(result_rows, columns=RESULT_COLUMNS)
    return results, errors


def run_comparison(data, output_path):
    """Calculate genomic comparisons from validated data and save a CSV."""

    results, comparison_errors = compare_samples(data)

    try:
        results.to_csv(output_path, index=False)
    except Exception as error:
        print(f"ERROR: Could not write the output CSV: {error}", file=sys.stderr)
        return 2

    group_sizes = data.groupby(
        ["Pathogen", "Gene_Segment"],
        sort=True,
        dropna=False,
    ).size()

    group_descriptions = []
    for group_names, count in group_sizes.items():
        pathogen, gene_segment = group_names
        sample_word = "sample" if count == 1 else "samples"
        description = f"{pathogen} / {gene_segment} ({count} {sample_word})"
        group_descriptions.append(description)

    group_text = "; ".join(group_descriptions)
    if not group_text:
        group_text = "none"

    print(f"Number of samples: {len(data)}")
    print(f"Pathogen/gene groups ({len(group_sizes)}): {group_text}")
    print(f"Pairwise comparisons: {len(results)}")

    for error in comparison_errors:
        print(f"ERROR: {error}", file=sys.stderr)

    if comparison_errors:
        return 1
    return 0


def main():
    """Load the two source files and run the sequence comparison."""

    parser = argparse.ArgumentParser(
        description="Align and compare sequences within pathogen/gene groups."
    )
    parser.add_argument("fasta_file", type=Path, help="Input .fasta, .fa, or .fna")
    parser.add_argument("metadata_file", type=Path, help="Input .xlsx or .csv metadata")
    parser.add_argument("--output", "-o", type=Path, help="Output .csv file")
    parser.add_argument("--default-gene-segment", help="Region for unannotated samples")
    parser.add_argument(
        "--confirm-same-region", action="store_true",
        help="Confirm all sequences represent the same region (one pathogen only)",
    )
    arguments = parser.parse_args()

    fasta_path = arguments.fasta_file
    metadata_path = arguments.metadata_file

    if arguments.output:
        output_path = arguments.output
    else:
        output_name = f"{metadata_path.stem}_pairwise_comparisons.csv"
        output_path = metadata_path.with_name(output_name)

    if fasta_path.suffix.lower() not in FASTA_EXTENSIONS:
        print("ERROR: FASTA must use .fasta, .fa, or .fna.", file=sys.stderr)
        return 2

    if not fasta_path.is_file():
        print("ERROR: The FASTA file does not exist.", file=sys.stderr)
        return 2

    if metadata_path.suffix.lower() not in METADATA_EXTENSIONS:
        print("ERROR: Metadata must use .xlsx or .csv.", file=sys.stderr)
        return 2

    if not metadata_path.is_file():
        print("ERROR: The metadata file does not exist.", file=sys.stderr)
        return 2

    if output_path.suffix.lower() != ".csv":
        print("ERROR: The output file must have the .csv extension.", file=sys.stderr)
        return 2

    try:
        data, validation_errors, _ = load_and_validate_inputs(
            fasta_path,
            metadata_path,
            arguments.default_gene_segment,
            arguments.confirm_same_region,
        )
    except Exception as error:
        print(f"ERROR: Could not read the input files: {error}", file=sys.stderr)
        return 2

    if validation_errors:
        for error in validation_errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    try:
        aligned_data = align_samples(data)
    except RuntimeError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    return run_comparison(aligned_data, output_path)


if __name__ == "__main__":
    raise SystemExit(main())
