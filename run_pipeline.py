"""Run the genomic epidemiology pipeline from FASTA and metadata inputs."""

import argparse
import math
import sys
from pathlib import Path

from add_temporal_geographic import run_stage as add_time_and_geography
from align_sequences import align_samples, find_mafft
from compare_genetic_distances import run_comparison
from detect_potential_clusters import run_stage as detect_clusters
from geocode_locations import fill_missing_coordinates
from validate_inputs import (
    FASTA_EXTENSIONS,
    METADATA_EXTENSIONS,
    load_and_validate_inputs,
)


def output_paths(metadata_path):
    """Create the four output paths from the metadata filename."""

    genomic_name = f"{metadata_path.stem}_pairwise_comparisons.csv"
    genomic_path = metadata_path.with_name(genomic_name)

    enriched_name = f"{genomic_path.stem}_temporal_geographic.csv"
    enriched_path = genomic_path.with_name(enriched_name)

    related_name = f"{enriched_path.stem}_related_pairs.csv"
    related_path = enriched_path.with_name(related_name)

    cluster_name = f"{enriched_path.stem}_cluster_membership.csv"
    clusters_path = enriched_path.with_name(cluster_name)

    return genomic_path, enriched_path, related_path, clusters_path


def region_assignment_path(metadata_path):
    """Name the small CSV that explains each sample's region label."""

    return metadata_path.with_name(f"{metadata_path.stem}_region_assignments.csv")


def run_pipeline(
    fasta_path,
    metadata_path,
    max_p_distance,
    max_days_apart,
    max_km_apart,
    default_gene_segment=None,
    confirm_same_region=False,
    complete_genomes=False,
    run_summary=None,
):
    """Run the backend stages and stop when a stage reports a problem."""

    output_files = output_paths(metadata_path)
    genomic_path = output_files[0]
    enriched_path = output_files[1]
    related_path = output_files[2]
    clusters_path = output_files[3]

    if complete_genomes and (default_gene_segment or confirm_same_region):
        print(
            "ERROR: Complete-genome mode cannot be combined with another "
            "genomic-region override.",
            file=sys.stderr,
        )
        return 2

    print("1/6 Validating and matching FASTA with metadata...")
    try:
        data, validation_errors, _ = load_and_validate_inputs(
            fasta_path,
            metadata_path,
            default_gene_segment,
            confirm_same_region,
            complete_genomes,
        )
    except Exception as error:
        print(f"ERROR: Could not validate the input files: {error}", file=sys.stderr)
        return 2

    if validation_errors:
        for error in validation_errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    # Check the required external program before optional network requests.
    try:
        mafft_path = find_mafft()
    except RuntimeError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    print("2/6 Filling missing coordinates with GeoNames when available...")
    data = fill_missing_coordinates(data)

    print("3/6 Aligning each pathogen/gene group with MAFFT...")
    try:
        aligned_data = align_samples(data, mafft_path)
    except RuntimeError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    print("4/6 Calculating genetic distances...")
    genetic_status = run_comparison(aligned_data, genomic_path)
    if genetic_status != 0:
        return 1

    print("5/6 Adding temporal and geographic distances...")
    enrichment_status = add_time_and_geography(
        data,
        genomic_path,
        enriched_path,
    )
    if enrichment_status != 0:
        return 1

    print("6/6 Detecting potential clusters...")
    cluster_status = detect_clusters(
        enriched_path,
        max_p_distance,
        max_days_apart,
        max_km_apart,
    )
    if cluster_status != 0:
        return 1

    assignments_path = region_assignment_path(metadata_path)
    try:
        data[["Sample_ID", "Pathogen", "Gene_Segment", "Gene_Segment_Source"]].to_csv(
            assignments_path, index=False
        )
    except Exception as error:
        print(f"ERROR: Could not save region assignments: {error}", file=sys.stderr)
        return 2

    print("Final output files:")
    print(f"Related pairs: {related_path.resolve()}")
    print(f"Cluster membership: {clusters_path.resolve()}")
    print(f"Region assignments: {assignments_path.resolve()}")
    if run_summary is not None:
        group_sizes = data.groupby(["Pathogen", "Gene_Segment"]).size()
        coarse = data["Geographic_Precision"].map(
            lambda value: "too coarse" in str(value)
        )
        missing_coordinates = data["Latitude"].isna() | data["Longitude"].isna()
        geonames_locations = {
            " ".join(str(name).split()).casefold()
            for name in data.loc[
                data["Coordinate_Source"] == "GeoNames", "Location_Name"
            ]
        }
        coarse_locations = {
            " ".join(str(name).split()).casefold()
            for name in data.loc[coarse, "Location_Name"]
        }
        run_summary.update(
            aligned_sequences=int(group_sizes[group_sizes > 1].sum()),
            geonames_resolved=len(geonames_locations),
            coarse_locations=len(coarse_locations),
            insufficient_geography_samples=int((coarse | missing_coordinates).sum()),
        )
    return 0


def main():
    """Read source paths and thresholds from the command line."""

    parser = argparse.ArgumentParser(
        description="Run the complete backend pipeline."
    )
    parser.add_argument("fasta_file", type=Path, help="Input .fasta, .fa, or .fna")
    parser.add_argument("metadata_file", type=Path, help="Input .xlsx or .csv metadata")
    parser.add_argument("--max-p-distance", type=float, required=True)
    parser.add_argument("--max-days-apart", type=float, required=True)
    parser.add_argument("--max-km-apart", type=float, required=True)
    parser.add_argument("--default-gene-segment", help="Region for unannotated samples")
    parser.add_argument(
        "--confirm-same-region", action="store_true",
        help="Confirm all sequences represent the same region (one pathogen only)",
    )
    parser.add_argument(
        "--complete-genomes", action="store_true",
        help="Confirm all sequences are complete genomes of one pathogen",
    )
    arguments = parser.parse_args()

    fasta_path = arguments.fasta_file
    metadata_path = arguments.metadata_file
    thresholds = [
        arguments.max_p_distance,
        arguments.max_days_apart,
        arguments.max_km_apart,
    ]

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

    for threshold in thresholds:
        invalid_threshold = threshold < 0 or not math.isfinite(threshold)
        if invalid_threshold:
            print(
                "ERROR: Thresholds must be finite and cannot be negative.",
                file=sys.stderr,
            )
            return 2

    options = {
        "default_gene_segment": arguments.default_gene_segment,
        "confirm_same_region": arguments.confirm_same_region,
    }
    if arguments.complete_genomes:
        options["complete_genomes"] = True

    return run_pipeline(
        fasta_path,
        metadata_path,
        *thresholds,
        **options,
    )


if __name__ == "__main__":
    raise SystemExit(main())
