"""Add temporal and geographic measurements to genomic pairwise results."""

import argparse
import math
import sys
from pathlib import Path

import pandas as pd

from geocode_locations import is_state_only_location
from validate_inputs import METADATA_EXTENSIONS, load_and_validate_metadata


def haversine_km(latitude_1, longitude_1, latitude_2, longitude_2):
    """Calculate the great-circle distance between two coordinates."""

    earth_radius_km = 6371.0

    # Python's trigonometry functions use radians instead of degrees.
    latitude_1 = math.radians(latitude_1)
    longitude_1 = math.radians(longitude_1)
    latitude_2 = math.radians(latitude_2)
    longitude_2 = math.radians(longitude_2)

    latitude_change = latitude_2 - latitude_1
    longitude_change = longitude_2 - longitude_1

    first_part = math.sin(latitude_change / 2) ** 2
    second_part = (
        math.cos(latitude_1)
        * math.cos(latitude_2)
        * math.sin(longitude_change / 2) ** 2
    )
    haversine_value = first_part + second_part

    # Keep tiny floating-point errors inside the valid range of 0 to 1.
    haversine_value = max(0.0, haversine_value)
    haversine_value = min(1.0, haversine_value)

    central_angle = 2 * math.atan2(
        math.sqrt(haversine_value),
        math.sqrt(1 - haversine_value),
    )
    return earth_radius_km * central_angle


def enrich_pairs(metadata, pairs):
    """Add days, kilometers, and same-location values to every sample pair."""

    enriched = pairs.copy()
    samples = metadata.copy()
    samples["Sample_ID"] = samples["Sample_ID"].astype("string").str.strip()
    samples = samples.set_index("Sample_ID")

    days_apart = []
    kilometers_apart = []
    same_location = []
    latitudes_1 = []
    longitudes_1 = []
    latitudes_2 = []
    longitudes_2 = []
    coordinate_sources_1 = []
    coordinate_sources_2 = []
    geographic_precisions_1 = []
    geographic_precisions_2 = []

    for _, pair in enriched.iterrows():
        sample_id_1 = str(pair["Sample_ID_1"]).strip()
        sample_id_2 = str(pair["Sample_ID_2"]).strip()

        sample_1_exists = sample_id_1 in samples.index
        sample_2_exists = sample_id_2 in samples.index

        if not sample_1_exists or not sample_2_exists:
            raise ValueError(
                f"Could not find metadata for pair {sample_id_1} and {sample_id_2}."
            )

        sample_1 = samples.loc[sample_id_1]
        sample_2 = samples.loc[sample_id_2]

        # Calculate the absolute number of days between collection dates.
        date_1 = sample_1["Collection_Date"]
        date_2 = sample_2["Collection_Date"]

        if pd.notna(date_1) and pd.notna(date_2):
            difference = date_1 - date_2
            days_apart.append(abs(difference.days))
        else:
            days_apart.append(None)

        # Calculate kilometers only when all four coordinates are available.
        latitude_1 = sample_1["Latitude"]
        longitude_1 = sample_1["Longitude"]
        latitude_2 = sample_2["Latitude"]
        longitude_2 = sample_2["Longitude"]

        # Preserve the actual coordinates used, including GeoNames estimates.
        latitudes_1.append(latitude_1)
        longitudes_1.append(longitude_1)
        latitudes_2.append(latitude_2)
        longitudes_2.append(longitude_2)

        location_1 = sample_1["Location_Name"]
        location_2 = sample_2["Location_Name"]
        precision_1 = sample_1.get("Geographic_Precision")
        precision_2 = sample_2.get("Geographic_Precision")
        if pd.isna(precision_1):
            if pd.notna(latitude_1) and pd.notna(longitude_1):
                precision_1 = "Provided coordinates (precision unverified)"
            else:
                precision_1 = "Named location only (precision unverified)"
        if pd.isna(precision_2):
            if pd.notna(latitude_2) and pd.notna(longitude_2):
                precision_2 = "Provided coordinates (precision unverified)"
            else:
                precision_2 = "Named location only (precision unverified)"

        coarse_1 = is_state_only_location(location_1) or "too coarse" in precision_1
        coarse_2 = is_state_only_location(location_2) or "too coarse" in precision_2
        if coarse_1:
            precision_1 = "State/administrative area (too coarse for distance cutoff)"
        if coarse_2:
            precision_2 = "State/administrative area (too coarse for distance cutoff)"
        geographic_precisions_1.append(precision_1)
        geographic_precisions_2.append(precision_2)

        coordinates_are_available = (
            pd.notna(latitude_1)
            and pd.notna(longitude_1)
            and pd.notna(latitude_2)
            and pd.notna(longitude_2)
            and not coarse_1
            and not coarse_2
        )

        if coordinates_are_available:
            distance = haversine_km(
                float(latitude_1),
                float(longitude_1),
                float(latitude_2),
                float(longitude_2),
            )
            kilometers_apart.append(round(distance, 3))
        else:
            kilometers_apart.append(None)

        # Keep the origin of each coordinate visible in the enriched CSV.
        source_1 = sample_1.get("Coordinate_Source")
        source_2 = sample_2.get("Coordinate_Source")
        if pd.isna(source_1) and pd.notna(latitude_1) and pd.notna(longitude_1):
            source_1 = "Provided"
        if pd.isna(source_2) and pd.notna(latitude_2) and pd.notna(longitude_2):
            source_2 = "Provided"
        coordinate_sources_1.append(source_1)
        coordinate_sources_2.append(source_2)

        # Compare location text without capitalization or extra-space differences.
        if pd.isna(location_1) or pd.isna(location_2) or coarse_1 or coarse_2:
            same_location.append(None)
        else:
            location_1 = " ".join(str(location_1).split()).casefold()
            location_2 = " ".join(str(location_2).split()).casefold()
            same_location.append(location_1 == location_2)

    enriched["Days_Apart"] = pd.array(days_apart, dtype="Int64")
    enriched["Km_Apart"] = kilometers_apart
    enriched["Same_Location"] = pd.array(same_location, dtype="boolean")
    enriched["Latitude_1"] = latitudes_1
    enriched["Longitude_1"] = longitudes_1
    enriched["Latitude_2"] = latitudes_2
    enriched["Longitude_2"] = longitudes_2
    enriched["Coordinate_Source_1"] = coordinate_sources_1
    enriched["Coordinate_Source_2"] = coordinate_sources_2
    enriched["Geographic_Precision_1"] = geographic_precisions_1
    enriched["Geographic_Precision_2"] = geographic_precisions_2
    return enriched


def run_stage(metadata, pairs_path, output_path):
    """Load pairwise results, enrich them with metadata, and save a new CSV."""

    try:
        pairs = pd.read_csv(pairs_path)
    except Exception as error:
        print(f"ERROR: Could not read the pairwise CSV: {error}", file=sys.stderr)
        return 2

    required_columns = ["Sample_ID_1", "Sample_ID_2"]
    missing_columns = []

    for column in required_columns:
        if column not in pairs.columns:
            missing_columns.append(column)

    if missing_columns:
        print(
            "ERROR: Pairwise CSV is missing column(s): "
            + ", ".join(missing_columns),
            file=sys.stderr,
        )
        return 1

    try:
        enriched = enrich_pairs(metadata, pairs)
        enriched.to_csv(output_path, index=False)
    except Exception as error:
        print(f"ERROR: Could not create the enriched CSV: {error}", file=sys.stderr)
        return 2

    valid_dates = enriched["Days_Apart"].notna().sum()
    valid_distances = enriched["Km_Apart"].notna().sum()

    print(f"Pairs processed: {len(enriched)}")
    print(f"Valid date comparisons: {valid_dates}")
    print(f"Valid geographic-distance comparisons: {valid_distances}")
    return 0


def main():
    """Read command-line paths and run the enrichment stage."""

    parser = argparse.ArgumentParser(
        description="Add temporal and geographic values to pairwise results."
    )
    parser.add_argument(
        "metadata_file",
        type=Path,
        help="Validated .xlsx or .csv metadata",
    )
    parser.add_argument(
        "pairwise_csv",
        type=Path,
        help="Genomic pairwise .csv file",
    )
    parser.add_argument("--output", "-o", type=Path, help="New enriched .csv file")
    arguments = parser.parse_args()

    metadata_path = arguments.metadata_file
    pairs_path = arguments.pairwise_csv

    if arguments.output:
        output_path = arguments.output
    else:
        output_name = f"{pairs_path.stem}_temporal_geographic.csv"
        output_path = pairs_path.with_name(output_name)

    if metadata_path.suffix.lower() not in METADATA_EXTENSIONS:
        print("ERROR: Metadata must use .xlsx or .csv.", file=sys.stderr)
        return 2

    if not metadata_path.is_file():
        print("ERROR: The metadata file does not exist.", file=sys.stderr)
        return 2

    if pairs_path.suffix.lower() != ".csv" or not pairs_path.is_file():
        print("ERROR: Pairwise input must be an existing .csv file.", file=sys.stderr)
        return 2

    if output_path.suffix.lower() != ".csv":
        print("ERROR: The output file must have the .csv extension.", file=sys.stderr)
        return 2

    if output_path.resolve() == pairs_path.resolve():
        print("ERROR: The output cannot overwrite the pairwise CSV.", file=sys.stderr)
        return 2

    try:
        metadata, validation_errors, _ = load_and_validate_metadata(metadata_path)
    except Exception as error:
        print(f"ERROR: Could not read the metadata file: {error}", file=sys.stderr)
        return 2

    # This standalone stage can leave dates or coordinates blank.
    allowed_error_starts = (
        "Missing Collection_Date",
        "Missing Latitude",
        "Missing Longitude",
    )
    remaining_errors = []

    for error in validation_errors:
        if not error.startswith(allowed_error_starts):
            remaining_errors.append(error)

    if remaining_errors:
        for error in remaining_errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    # The complete runner supplies FASTA Sample_ID values directly.
    # Standalone use assumes pairwise IDs match metadata accessions.
    metadata["Sample_ID"] = metadata["Accession"]
    return run_stage(metadata, pairs_path, output_path)


if __name__ == "__main__":
    raise SystemExit(main())
