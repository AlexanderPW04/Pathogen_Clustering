"""Turn approved proposals into derived metadata and recompute affected results.

The original metadata is never modified: the derived table keeps every
original value in an Original_* column beside the enriched one. Only date
and distance screening is recomputed, from the genetic distances already
calculated, so MAFFT and variant calling are not run again.
"""

import pandas as pd

from add_temporal_geographic import enrich_pairs
from detect_potential_clusters import build_clusters, enough_overlap, find_related_pairs
from geocode_locations import GeoNamesUnavailable
from literature.geography import STATE_ABBREVIATIONS, original_geography
from literature.verification import coordinates_from, parse_iso_date

TRACKED = ["Location_Name", "Latitude", "Longitude", "Coordinate_Source", "Geographic_Precision",
           "Collection_Date"]
EXTRA_FIELDS = ["Site_Identifier", "Sample_Type", "Collection_Method", "Study_Identifier", "Study_Context"]
FIELD_ORDER = {"Location": 0, "Coordinates": 1, "Collection_Date": 2, "Collection_Date_Range": 3}
PAIR_GEOGRAPHY = ["Days_Apart", "Km_Apart", "Same_Location", "Latitude_1", "Longitude_1", "Latitude_2",
                  "Longitude_2", "Coordinate_Source_1", "Coordinate_Source_2", "Geographic_Precision_1",
                  "Geographic_Precision_2"]
LITERATURE_PLACE = "GeoNames place for an approved literature location (approximate)"
NAME_ONLY = "Approved literature place name (no coordinates; precision unverified)"
PUBLICATION_COORDINATES = "Publication-reported site coordinates (precision unverified)"


def apply_approved(samples, facts, decisions, geocoder=None):
    """Return (derived metadata, application log) for approved, approvable facts.

    geocoder(query) returns a dict from geocode_locations.lookup_place, None
    when no place matches, or raises GeoNamesUnavailable.
    """

    derived = samples.copy()
    for column in TRACKED:
        if column not in derived:
            derived[column] = pd.NA
        derived[f"Original_{column}"] = derived[column]
    for column in ["Enriched_Date_Range_Start", "Enriched_Date_Range_End", "Enrichment_Sources"] + \
            [f"Enriched_{field}" for field in EXTRA_FIELDS]:
        derived[column] = pd.NA
    derived["Coordinate_Source"] = derived["Coordinate_Source"].astype(object)
    derived["Geographic_Precision"] = derived["Geographic_Precision"].astype(object)
    log = []
    approved = [fact for fact in facts if fact.get("Approvable")
                and decisions.get(fact["Fact_ID"]) == "approve"]
    approved.sort(key=lambda fact: FIELD_ORDER.get(fact["Field"], 9))
    index_of = {str(sample_id): index for index, sample_id in derived["Sample_ID"].items()}
    geonames_note = None
    for fact in approved:
        index = index_of.get(str(fact["Sample_ID"]))
        if index is None:
            continue
        row = samples.loc[index]
        field, value = fact["Field"], str(fact["Value"]).strip()
        outcome = "Applied"
        if field == "Location":
            outcome, geonames_note = apply_location(derived, index, row, value, geocoder, geonames_note)
        elif field == "Coordinates":
            latitude, longitude = coordinates_from(value)
            derived.at[index, "Latitude"] = latitude
            derived.at[index, "Longitude"] = longitude
            derived.at[index, "Coordinate_Source"] = "Publication (approved literature evidence)"
            derived.at[index, "Geographic_Precision"] = PUBLICATION_COORDINATES
        elif field == "Collection_Date":
            derived.at[index, "Collection_Date"] = parse_iso_date(value)
        elif field == "Collection_Date_Range":
            start, end = value.split("/")
            derived.at[index, "Enriched_Date_Range_Start"] = start.strip()
            derived.at[index, "Enriched_Date_Range_End"] = end.strip()
            outcome = "Applied (kept as a range; exact-day screening is unchanged)"
        else:
            derived.at[index, f"Enriched_{field}"] = value
            outcome = "Applied (descriptive; no calculation uses it)"
        if outcome.startswith("Applied"):
            source = f"{field} from {fact['Article_ID']} ({fact['Passage_ID']})"
            previous = derived.at[index, "Enrichment_Sources"]
            derived.at[index, "Enrichment_Sources"] = source if pd.isna(previous) else f"{previous}; {source}"
        log.append({"Fact_ID": fact["Fact_ID"], "Sample_ID": fact["Sample_ID"], "Field": field,
                    "Value": value, "Outcome": outcome})
    return derived, log


def apply_location(derived, index, row, value, geocoder, geonames_note):
    """Set an approved place name; resolve it with GeoNames when it is configured."""

    original = original_geography(row)
    query = value
    if original["state"] and original["state"] not in value.casefold():
        query = f"{value}, {original['state'].title()}"
    place = None
    if geocoder is not None and geonames_note is None:
        try:
            place = geocoder(query)
        except GeoNamesUnavailable as problem:
            geonames_note = str(problem)
    if place is not None and original["country"] == "usa" and place["country_code"] not in {"", "US"}:
        return (f"Not applied: GeoNames placed '{query}' in {place['country']}, which conflicts with the "
                f"original country."), geonames_note
    if place is not None and original["state"]:
        found_state = place["state"].casefold()
        if found_state and STATE_ABBREVIATIONS.get(found_state, found_state) != original["state"]:
            return (f"Not applied: GeoNames placed '{query}' in {place['state']}, which conflicts with the "
                    f"original state."), geonames_note
    derived.at[index, "Location_Name"] = value
    # Coarse original coordinates (for example a state centroid) must not be
    # reused with the more specific name.
    derived.at[index, "Latitude"] = place["latitude"] if place else pd.NA
    derived.at[index, "Longitude"] = place["longitude"] if place else pd.NA
    if place:
        precision = LITERATURE_PLACE if "too coarse" not in place["precision"] else place["precision"]
        derived.at[index, "Coordinate_Source"] = "GeoNames (approved literature location)"
        derived.at[index, "Geographic_Precision"] = precision
        return f"Applied (GeoNames: {place['name']}, {place['state']}; approximate)", geonames_note
    derived.at[index, "Coordinate_Source"] = pd.NA
    derived.at[index, "Geographic_Precision"] = NAME_ONLY
    reason = geonames_note or ("GeoNames found no match" if geocoder else "GeoNames not used")
    return f"Applied (place name only: {reason})", geonames_note


def unresolved_geography(pairs, thresholds):
    """Pairs meeting the genetic and date cutoffs whose geography cannot be judged."""

    count = 0
    for _, pair in pairs.iterrows():
        if pd.isna(pair["P_Distance"]) or pd.isna(pair["Days_Apart"]) or not enough_overlap(pair):
            continue
        if pair["P_Distance"] > thresholds[0] or pair["Days_Apart"] > thresholds[1]:
            continue
        if pd.isna(pair["Km_Apart"]) and str(pair["Same_Location"]).casefold() != "true":
            count += 1
    return count


def rescreen(original, derived, pairs, thresholds):
    """Before/after date and distance screening from the existing genetic distances."""

    genetic = pairs[[column for column in pairs.columns if column not in PAIR_GEOGRAPHY]]
    before = enrich_pairs(original, genetic)
    after = enrich_pairs(derived, genetic)
    related_before = find_related_pairs(before, *thresholds)
    related_after = find_related_pairs(after, *thresholds)
    key = ["Sample_ID_1", "Sample_ID_2", "Pathogen", "Gene_Segment"]
    before_keys = {tuple(map(str, row)) for row in related_before[key].itertuples(index=False)}
    after_keys = {tuple(map(str, row)) for row in related_after[key].itertuples(index=False)}
    changed_samples = set(derived.loc[derived["Enrichment_Sources"].notna(), "Sample_ID"].astype(str))
    rows = []
    for (_, old), (_, new) in zip(before.iterrows(), after.iterrows()):
        pair_key = tuple(str(old[column]) for column in key)
        was, now = pair_key in before_keys, pair_key in after_keys
        differs = any(str(old[column]) != str(new[column]) for column in ["Days_Apart", "Km_Apart", "Same_Location"])
        if was == now:
            change = "Values changed; screening result unchanged" if differs else "Unchanged"
        else:
            change = "Newly meets all screening criteria" if now else "No longer meets all criteria"
        enriched = sorted({pair_key[0], pair_key[1]} & changed_samples)
        rows.append({
            "Sample_ID_1": pair_key[0], "Sample_ID_2": pair_key[1], "Pathogen": pair_key[2],
            "Gene_Segment": pair_key[3], "P_Distance": old["P_Distance"],
            "Days_Apart_Before": old["Days_Apart"], "Days_Apart_After": new["Days_Apart"],
            "Km_Apart_Before": old["Km_Apart"], "Km_Apart_After": new["Km_Apart"],
            "Same_Location_Before": old["Same_Location"], "Same_Location_After": new["Same_Location"],
            "Meets_Criteria_Before": was, "Meets_Criteria_After": now, "Change": change,
            "Reason": (f"Approved enrichment of {', '.join(enriched)}" if enriched and change != "Unchanged"
                       else ""),
        })
    table = pd.DataFrame(rows)
    summary = {
        "Pairs meeting all screening criteria": (len(related_before), len(related_after)),
        "Potential clusters": (build_clusters(related_before)["Cluster_ID"].nunique(),
                               build_clusters(related_after)["Cluster_ID"].nunique()),
        "Pairs with a kilometre distance": (int(before["Km_Apart"].notna().sum()),
                                            int(after["Km_Apart"].notna().sum())),
        "Genetic + date matches with unresolved geography": (unresolved_geography(before, thresholds),
                                                             unresolved_geography(after, thresholds)),
    }
    return table, summary


def precision_counts(samples):
    """How many samples have each geographic precision label."""

    labels = samples["Geographic_Precision"].fillna("No location given").astype(str)
    return labels.value_counts().to_dict()
