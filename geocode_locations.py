"""Fill missing sample coordinates with approximate GeoNames place coordinates."""

import math
import os
import sys

import pandas as pd
import requests


GEONAMES_URL = "https://secure.geonames.org/searchJSON"


def is_state_only_location(location):
    """Recognize names like 'USA: ND' without treating them as collection sites."""

    if pd.isna(location):
        return False
    country, separator, area = " ".join(str(location).split()).partition(":")
    return (
        bool(separator)
        and country.casefold() in {"usa", "us", "united states"}
        and len(area.strip()) == 2
        and area.strip().isalpha()
    )


def fill_missing_coordinates(data):
    """Keep provided coordinates and try GeoNames for incomplete locations."""

    samples = data.copy()
    samples["Coordinate_Source"] = pd.NA
    samples["Geographic_Precision"] = "Named location only (precision unverified)"

    no_location = samples["Location_Name"].isna() | samples["Location_Name"].astype(
        "string").str.strip().eq("").fillna(True)
    samples.loc[no_location, "Geographic_Precision"] = "No location given"

    provided = samples["Latitude"].notna() & samples["Longitude"].notna()
    samples.loc[provided, "Coordinate_Source"] = "Provided"
    samples.loc[provided, "Geographic_Precision"] = (
        "Provided coordinates (precision unverified)"
    )

    state_only = samples["Location_Name"].map(is_state_only_location)
    samples.loc[state_only, "Geographic_Precision"] = (
        "State-level location (too coarse for distance cutoff)"
    )
    if state_only.any():
        print(
            f"State-level locations: {int(state_only.sum())}; their coordinates "
            "and same-location names will not satisfy the geographic cutoff.",
            file=sys.stderr,
        )

    missing_count = int((~provided & ~state_only).sum())
    if missing_count == 0:
        return samples

    username = os.getenv("GEONAMES_USERNAME", "").strip()
    if not username:
        print(
            "GeoNames lookup skipped: GEONAMES_USERNAME is not set. "
            "Specific place names may use the Same_Location fallback.",
            file=sys.stderr,
        )
        return samples

    cache = {}
    geocoded_count = 0
    coarse_count = 0

    for index, sample in samples.loc[~provided & ~state_only].iterrows():
        location = sample["Location_Name"]
        if pd.isna(location):
            continue

        query = " ".join(str(location).split())
        if not query:
            continue

        cache_key = query.casefold()
        if cache_key not in cache:
            try:
                response = requests.get(
                    GEONAMES_URL,
                    params={"q": query, "maxRows": 1, "username": username},
                    timeout=8,
                    verify=True,
                )
                response.raise_for_status()
                result = response.json()
            except (requests.RequestException, ValueError):
                print(
                    "GeoNames lookup unavailable. Remaining missing coordinates "
                    "may use the Same_Location fallback for specific places.",
                    file=sys.stderr,
                )
                break

            if not isinstance(result, dict) or "status" in result:
                print(
                    "GeoNames lookup unavailable. Remaining missing coordinates "
                    "may use the Same_Location fallback for specific places.",
                    file=sys.stderr,
                )
                break

            cache[cache_key] = None
            places = result.get("geonames", [])
            if places:
                try:
                    latitude = float(places[0]["lat"])
                    longitude = float(places[0]["lng"])
                except (KeyError, TypeError, ValueError, IndexError):
                    continue

                valid_latitude = math.isfinite(latitude) and -90 <= latitude <= 90
                valid_longitude = math.isfinite(longitude) and -180 <= longitude <= 180
                if valid_latitude and valid_longitude:
                    feature_code = str(places[0].get("fcode", "")).upper()
                    feature_class = str(places[0].get("fcl", "")).upper()
                    if feature_class == "A" or feature_code.startswith(
                        ("ADM", "PCL")
                    ):
                        precision = "Administrative area (too coarse for distance cutoff)"
                    else:
                        precision = "GeoNames place (approximate)"
                    cache[cache_key] = (latitude, longitude, precision)

        coordinates = cache[cache_key]
        if coordinates is None:
            continue

        # A partial original coordinate cannot be used for a distance, so both
        # values are replaced together with the approximate place coordinates.
        samples.at[index, "Latitude"] = coordinates[0]
        samples.at[index, "Longitude"] = coordinates[1]
        samples.at[index, "Coordinate_Source"] = "GeoNames"
        samples.at[index, "Geographic_Precision"] = coordinates[2]
        geocoded_count += 1
        if "too coarse" in coordinates[2]:
            coarse_count += 1

    unresolved_count = int(
        (samples["Latitude"].isna() | samples["Longitude"].isna()).sum()
    )
    print(
        f"GeoNames coordinates added: {geocoded_count}; "
        f"samples still missing coordinates: {unresolved_count}."
    )
    if coarse_count:
        print(
            f"GeoNames administrative-area estimates: {coarse_count}; these "
            "are too coarse for the distance cutoff and same-location fallback.",
            file=sys.stderr,
        )
    return samples


def check_geonames():
    """Make one short live lookup without revealing the configured username."""

    username = os.getenv("GEONAMES_USERNAME", "").strip()
    if not username:
        return {
            "status": "Not configured",
            "message": "GEONAMES_USERNAME is not set on this server.",
        }

    try:
        response = requests.get(
            GEONAMES_URL,
            params={"q": "Tempe, Arizona, USA", "maxRows": 1, "username": username},
            timeout=4,
            verify=True,
        )
        response.raise_for_status()
    except requests.HTTPError:
        return {"status": "Check failed", "message": "GeoNames returned an API error."}
    except requests.RequestException:
        return {"status": "Unavailable", "message": "Could not connect to GeoNames."}

    try:
        result = response.json()
        if not isinstance(result, dict) or "status" in result:
            return {"status": "Check failed", "message": "GeoNames rejected the lookup."}
        place = result["geonames"][0]
        latitude = float(place["lat"])
        longitude = float(place["lng"])
    except (ValueError, TypeError, KeyError, IndexError):
        return {
            "status": "Check failed",
            "message": "GeoNames did not return valid coordinates.",
        }

    valid_latitude = math.isfinite(latitude) and -90 <= latitude <= 90
    valid_longitude = math.isfinite(longitude) and -180 <= longitude <= 180
    if not (valid_latitude and valid_longitude):
        return {
            "status": "Check failed",
            "message": "GeoNames did not return valid coordinates.",
        }

    return {
        "status": "Operational",
        "message": "A Tempe, Arizona lookup returned valid approximate place coordinates.",
    }
