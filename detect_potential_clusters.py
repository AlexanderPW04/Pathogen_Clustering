"""Detect potential outbreak clusters from enriched sample-pair data."""

import argparse
import math
import sys
from pathlib import Path

import pandas as pd


PAIR_COLUMNS = [
    "Sample_ID_1",
    "Sample_ID_2",
    "Pathogen",
    "Gene_Segment",
    "P_Distance",
    "Differences",
    "Days_Apart",
    "Km_Apart",
    "Same_Location",
]

CLUSTER_COLUMNS = [
    "Cluster_ID",
    "Sample_ID",
    "Pathogen",
    "Gene_Segment",
]


def find_related_pairs(
    pairs,
    max_p_distance,
    max_days_apart,
    max_km_apart,
):
    """Keep pairs that pass the genetic, time, and location criteria."""

    related_rows = []

    for _, pair in pairs.iterrows():
        p_distance = pair["P_Distance"]
        days_apart = pair["Days_Apart"]
        km_apart = pair["Km_Apart"]

        # Genetic and date information are always required.
        if pd.isna(p_distance) or pd.isna(days_apart):
            continue

        genetic_passes = p_distance <= max_p_distance
        time_passes = days_apart <= max_days_apart

        if not genetic_passes or not time_passes:
            continue

        # Coordinates take priority. Location-name matching is the fallback.
        same_location_text = str(pair["Same_Location"]).strip().casefold()
        same_location = same_location_text == "true"

        if pd.notna(km_apart):
            geography_passes = km_apart <= max_km_apart
        else:
            geography_passes = same_location

        if not geography_passes:
            continue

        related_row = {}
        for column in PAIR_COLUMNS:
            related_row[column] = pair[column]
        related_rows.append(related_row)

    return pd.DataFrame(related_rows, columns=PAIR_COLUMNS)


def build_clusters(related_pairs):
    """Use depth-first search to find connected samples within each group."""

    membership_rows = []
    cluster_number = 1

    groups = related_pairs.groupby(
        ["Pathogen", "Gene_Segment"],
        sort=True,
        dropna=False,
    )

    for group_names, group in groups:
        pathogen, gene_segment = group_names

        # Each dictionary key is a sample. Its set contains connected samples.
        connections = {}

        for _, pair in group.iterrows():
            sample_1 = str(pair["Sample_ID_1"]).strip()
            sample_2 = str(pair["Sample_ID_2"]).strip()

            if sample_1 not in connections:
                connections[sample_1] = set()
            if sample_2 not in connections:
                connections[sample_2] = set()

            connections[sample_1].add(sample_2)
            connections[sample_2].add(sample_1)

        visited = set()

        for starting_sample in sorted(connections):
            if starting_sample in visited:
                continue

            # A stack is a simple list of samples still waiting to be checked.
            stack = [starting_sample]
            cluster_members = []

            while stack:
                sample = stack.pop()

                if sample in visited:
                    continue

                visited.add(sample)
                cluster_members.append(sample)

                for neighbor in connections[sample]:
                    if neighbor not in visited:
                        stack.append(neighbor)

            if len(cluster_members) < 2:
                continue

            cluster_id = f"C{cluster_number}"
            cluster_number += 1

            for sample in sorted(cluster_members):
                membership_rows.append(
                    {
                        "Cluster_ID": cluster_id,
                        "Sample_ID": sample,
                        "Pathogen": pathogen,
                        "Gene_Segment": gene_segment,
                    }
                )

    return pd.DataFrame(membership_rows, columns=CLUSTER_COLUMNS)


def run_stage(
    input_path,
    max_p_distance,
    max_days_apart,
    max_km_apart,
):
    """Load enriched pairs, detect potential clusters, and save two CSVs."""

    try:
        pairs = pd.read_csv(input_path)
    except Exception as error:
        print(f"ERROR: Could not read the enriched CSV: {error}", file=sys.stderr)
        return 2

    missing_columns = []
    for column in PAIR_COLUMNS:
        if column not in pairs.columns:
            missing_columns.append(column)

    if missing_columns:
        print(
            "ERROR: Enriched CSV is missing column(s): "
            + ", ".join(missing_columns),
            file=sys.stderr,
        )
        return 1

    numeric_columns = ["P_Distance", "Days_Apart", "Km_Apart"]
    for column in numeric_columns:
        pairs[column] = pd.to_numeric(pairs[column], errors="coerce")

    related_pairs = find_related_pairs(
        pairs,
        max_p_distance,
        max_days_apart,
        max_km_apart,
    )
    cluster_membership = build_clusters(related_pairs)

    related_name = f"{input_path.stem}_related_pairs.csv"
    cluster_name = f"{input_path.stem}_cluster_membership.csv"
    related_path = input_path.with_name(related_name)
    clusters_path = input_path.with_name(cluster_name)

    try:
        related_pairs.to_csv(related_path, index=False)
        cluster_membership.to_csv(clusters_path, index=False)
    except Exception as error:
        print(f"ERROR: Could not write an output CSV: {error}", file=sys.stderr)
        return 2

    number_of_clusters = cluster_membership["Cluster_ID"].nunique()

    print(f"Total pairs evaluated: {len(pairs)}")
    print(f"Related pairs passing thresholds: {len(related_pairs)}")
    print(f"Potential outbreak clusters: {number_of_clusters}")
    print(f"Samples included in clusters: {len(cluster_membership)}")
    return 0


def main():
    """Read the CSV path and prototype thresholds from the command line."""

    parser = argparse.ArgumentParser(
        description="Detect potential outbreak clusters using prototype thresholds."
    )
    parser.add_argument("enriched_csv", type=Path, help="Enriched pairwise .csv file")
    parser.add_argument("--max-p-distance", type=float, required=True)
    parser.add_argument("--max-days-apart", type=float, required=True)
    parser.add_argument("--max-km-apart", type=float, required=True)
    arguments = parser.parse_args()

    input_path = arguments.enriched_csv
    thresholds = [
        arguments.max_p_distance,
        arguments.max_days_apart,
        arguments.max_km_apart,
    ]

    if input_path.suffix.lower() != ".csv" or not input_path.is_file():
        print("ERROR: The input must be an existing .csv file.", file=sys.stderr)
        return 2

    for threshold in thresholds:
        if threshold < 0 or not math.isfinite(threshold):
            print(
                "ERROR: Thresholds must be finite and cannot be negative.",
                file=sys.stderr,
            )
            return 2

    return run_stage(input_path, *thresholds)


if __name__ == "__main__":
    raise SystemExit(main())
