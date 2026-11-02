"""Run the integrated individual + pooled analysis and collect its results.

This module contains no Flask code. The website and the command line both
call `run_integrated_analysis()`, which:

1. loads the reference genome (when supplied);
2. builds reference-based variant profiles for individual sequences;
3. imports pooled SNV records, optional coverage, and applies paper QC;
4. runs the pooled analyses (SNVs per sample, depth regression, genome
   positions, temporal trends, population structure);
5. compares the two evidence sources when both are present and compatible;
6. draws the figures and records the provenance of every step.

Each analysis that lacks the information it needs is skipped with a message.
"""

import platform
from importlib.metadata import PackageNotFoundError, version

import numpy as np
import pandas as pd

from geocode_locations import is_state_only_location
from genomics import figures
from genomics.diversity_analysis import EXPLORATORY_ORDINATION, population_structure
from genomics.pooled_observations import (
    build_depth_lookup,
    load_pooled_metadata,
    load_position_table,
    pooled_sample_records,
    pooled_variant_records,
)
from genomics.quality_control import (
    PAPER_MODE,
    PAPER_SETTINGS,
    apply_pooled_quality_control,
    load_masked_sites,
)
from genomics.references import build_reference_registry, load_reference_fasta
from genomics.sequence_variants import (
    MINIMUM_CALLABLE_POSITIONS,
    MINIMUM_MAPPING_IDENTITY,
    build_individual_profiles,
)
from genomics.temporal_analysis import (
    depth_regression,
    frequency_trends,
    location_date_summary,
    observation_matrix,
    variants_per_sample,
)
from genomics.unified_schema import INDIVIDUAL, POOLED
from genomics.variant_matching import (
    CATEGORY_ORDER,
    INCOMPARABLE,
    category_counts,
    compare_evidence,
    compare_signatures,
    load_signatures,
)
from wastewater.processed_snv import load_processed_snvs


DEFAULT_SETTINGS = {
    "qc_mode": PAPER_MODE,
    "date_window_days": 14,
    **PAPER_SETTINGS,
}

TRAJECTORY_VARIANTS = 6
MATRIX_VARIANTS = 40
HEATMAP_VARIANTS = 40


def package_version(name):
    try:
        return version(name)
    except PackageNotFoundError:
        return "not installed"


def first_text(row, columns):
    for column in columns:
        value = row.get(column)
        if pd.notna(value) and str(value).strip():
            return str(value).strip()
    return ""


def position_summary(pooled_variants, individual_variants, reference_id):
    """One row per variant on one reference: samples reporting, annotations."""

    rows = {}
    if pooled_variants is not None:
        included = pooled_variants[pooled_variants["Analysis_Included"]
                                   & pooled_variants["Reference_ID"].eq(reference_id)]
        for key, group in included.groupby("Variant_Key"):
            first = group.iloc[0].to_dict()
            rows[key] = {
                "Variant_Key": key, "Position": int(first["Position"]), "Ref": first["Ref"],
                "Alt": first["Alt"], "Pooled_Samples": group["Sample_ID"].nunique(),
                "Individual_Carriers": 0,
                "Median_Reported_Frequency": float(group["Alt_Frequency"].median()),
                "Gene_Or_Annotation": first_text(first, ["Published_Annotation", "Source_Gene_Feature"]),
                "Published_Novelty": first_text(first, ["Published_Novel_Global"]),
                "Published_Amino_Acid_Change": (
                    f"{first_text(first, ['Published_Ref_AA'])}>{first_text(first, ['Published_Alt_AA'])}"
                    if first_text(first, ["Published_Ref_AA"]) else ""
                ),
            }
    if individual_variants is not None:
        on_reference = individual_variants[individual_variants["Reference_ID"].eq(reference_id)]
        for key, group in on_reference.groupby("Variant_Key"):
            first = group.iloc[0]
            row = rows.setdefault(key, {
                "Variant_Key": key, "Position": int(first["Position"]), "Ref": first["Ref"],
                "Alt": first["Alt"], "Pooled_Samples": 0, "Median_Reported_Frequency": np.nan,
                "Gene_Or_Annotation": "", "Published_Novelty": "", "Published_Amino_Acid_Change": "",
            })
            row["Individual_Carriers"] = group["Sample_ID"].nunique()

    table = pd.DataFrame(list(rows.values()), columns=[
        "Variant_Key", "Position", "Ref", "Alt", "Pooled_Samples", "Individual_Carriers",
        "Median_Reported_Frequency", "Gene_Or_Annotation", "Published_Novelty",
        "Published_Amino_Acid_Change",
    ]).sort_values(["Position", "Alt"]).reset_index(drop=True)
    pooled_alts = table[table["Pooled_Samples"] > 0].groupby("Position")["Alt"].nunique()
    table["Multiple_Alt_At_Position"] = table["Position"].map(pooled_alts).fillna(0).gt(1)
    return table


def annotation_ranges(summary):
    """Approximate annotated regions from the positions of annotated records."""

    ranges = []
    annotated = summary[summary["Gene_Or_Annotation"].ne("")]
    for name, group in annotated.groupby("Gene_Or_Annotation"):
        ranges.append((name, int(group["Position"].min()), int(group["Position"].max()) + 1))
    return sorted(ranges, key=lambda item: item[1])


def add_figure(result, name, title, caption, make):
    """Draw one figure; a drawing problem becomes a visible message, not a crash."""

    try:
        result["figures"][name] = {"png": make(), "title": title, "caption": caption}
    except Exception as error:  # noqa: BLE001 - shown to the user as a limitation
        result["notes"].setdefault("figures", []).append(f"{title} could not be drawn: {error}")


def run_integrated_analysis(individual_data=None, mafft_path=None, reference_path=None,
                            snv_path=None, pooled_metadata_path=None, position_path=None,
                            masked_sites_path=None, signature_path=None, settings=None,
                            input_labels=None, relationship=None):
    """Return a dictionary of tables, figures, summaries, notes, and provenance.

    `relationship` optionally carries the Sprint 2 pairwise results
    ({"pairs", "thresholds", "candidates", "clusters"}) so their figures are
    drawn with the rest of the dashboard; they are not recalculated here.
    """

    settings = {**DEFAULT_SETTINGS, **(settings or {})}
    labels = input_labels or {}
    result = {"tables": {}, "figures": {}, "notes": {}, "summary": {}, "reference": None,
              "individual": None, "pooled": None, "population": None, "comparison": None,
              "settings": settings, "relationship": relationship}

    def note(section, text):
        result["notes"].setdefault(section, []).append(text)

    # 1. Reference genome.
    reference = None
    registry = {}
    if reference_path:
        reference = load_reference_fasta(reference_path)
        registry = build_reference_registry([reference])
        result["reference"] = reference

    # 2. Individual sequences -> reference-coordinate profiles.
    individual = None
    if individual_data is not None:
        if reference is None:
            note("variation", "No reference FASTA was supplied, so individual sequences have "
                              "genetic distances but no reference-coordinate variant profile.")
        else:
            individual = build_individual_profiles(
                individual_data, reference, mafft_path, labels.get("fasta", "individual FASTA")
            )
            result["individual"] = individual
            unmapped = individual["samples"][
                ~individual["samples"]["Variant_Profile_Status"].eq("Mapped to reference")
            ]
            for row in unmapped.to_dict("records"):
                note("variation", f"{row['Sample_ID']}: {row['Variant_Profile_Status']}. Its "
                                  "Sprint 2 genetic distances are kept, but it is not used in "
                                  "cross-source comparisons.")
            result["tables"]["individual_variants.csv"] = individual["variants"]
            result["tables"]["individual_samples.csv"] = individual["samples"]
            if not individual["indels"].empty:
                result["tables"]["individual_indels.csv"] = individual["indels"]

    # 3. Pooled observations and quality control.
    pooled = None
    if snv_path:
        snv_data, warnings = load_processed_snvs(snv_path, pooled_metadata_path)
        input_format = snv_data.attrs.get("input_format", "Processed SNV table")
        metadata = load_pooled_metadata(pooled_metadata_path) if pooled_metadata_path else None
        position_table = load_position_table(position_path) if position_path else None
        lengths = {reference["Reference_ID"]: reference["Length"]} if reference else {}
        source = labels.get("snv", "pooled SNV table")
        variants = pooled_variant_records(snv_data, source)
        samples = pooled_sample_records(snv_data, source, metadata, position_table, lengths,
                                        settings["coverage_depth"])
        masked = load_masked_sites(masked_sites_path) if masked_sites_path else None
        variants, samples, report = apply_pooled_quality_control(
            variants, samples, settings["qc_mode"], masked, settings, reference, registry
        )
        pooled = {"variants": variants, "samples": samples, "report": report,
                  "metadata_supplied": metadata is not None,
                  "warnings": warnings, "input_format": input_format,
                  "position_table": position_table,
                  "depth_lookup": build_depth_lookup(position_table)}
        result["pooled"] = pooled
        for warning in warnings:
            note("quality", warning)
        if position_table is None:
            note("quality", "No per-position depth table was supplied: unreported pooled variants "
                            "are 'coverage unknown', never 0%.")
        result["tables"]["pooled_variants_with_qc.csv"] = variants
        result["tables"]["pooled_samples_with_qc.csv"] = samples
        result["tables"]["paper_method_reproduction.csv"] = report

    # 4. Pooled analyses.
    if pooled is not None:
        variants = pooled["variants"]
        samples = pooled["samples"]
        per_sample = variants_per_sample(variants, samples)
        result["tables"]["variants_per_sample.csv"] = per_sample
        add_figure(result, "variants_per_sample", "SNVs per pooled sample",
                   "Analogue of the paper's Fig 2A. Hatched bars were excluded by quality control.",
                   lambda: figures.variants_per_sample_figure(per_sample))
        add_figure(result, "qc_outcomes", "Quality checks applied to pooled records",
                   "Each paper filter is applied only to values the files supply; hatched segments "
                   "are checks that could not be verified, never counted as passes.",
                   lambda: figures.qc_outcome_figure(variants))
        add_figure(result, "depth_distribution", "Read-depth distribution",
                   "Left: depth at each reported variant (dashed line: the paper's 20x minimum). "
                   "Right, when a per-position table is supplied: depth at every listed position.",
                   lambda: figures.depth_distribution_figure(
                       variants, settings["min_variant_depth"], pooled["position_table"]))

        fit = depth_regression(per_sample)
        result["regression"] = fit
        if fit is None:
            note("quality", "Depth-versus-SNV regression needs a mean genome depth for at least "
                            "three included pooled samples (from metadata or a depth table).")
        else:
            add_figure(result, "depth_vs_snvs", "Sequencing depth versus detected SNVs",
                       "Analogue of the paper's Fig 2B. The band is a 95% bootstrap interval of "
                       "the fitted line; the p-value is a permutation test.",
                       lambda: figures.depth_regression_figure(per_sample, fit))

        main_reference = (variants.loc[variants["Analysis_Included"], "Reference_ID"].mode()
                          if variants["Analysis_Included"].any() else pd.Series(dtype=str))
        reference_id = (reference["Reference_ID"] if reference is not None
                        else (main_reference.iloc[0] if not main_reference.empty else None))
        if reference_id is not None:
            summary = position_summary(variants, individual["variants"] if individual else None,
                                       reference_id)
            result["tables"]["variant_positions.csv"] = summary
            if not summary.empty:
                add_figure(result, "genome_positions", "Variants along the reference genome",
                           "Analogue of the paper's Fig 3. Annotation bands span the annotated "
                           "records only; they are not complete gene boundaries.",
                           lambda: figures.genome_position_figure(summary, annotation_ranges(summary)))

        matrix = observation_matrix(variants, samples, pooled["depth_lookup"],
                                    settings["min_variant_depth"])
        result["tables"]["pooled_observation_matrix.csv"] = matrix
        trends = frequency_trends(matrix)
        result["tables"]["temporal_trends.csv"] = trends
        result["tables"]["site_date_summary.csv"] = location_date_summary(variants, samples)
        # Figures show one reference at a time so positions are never mixed.
        figure_matrix = matrix[matrix["Reference_ID"].eq(reference_id)] if reference_id else matrix
        other_references = sorted(set(matrix["Reference_ID"]) - {reference_id}) if reference_id else []
        if other_references:
            note("temporal", "Figures show reference " + str(reference_id) + "; records on "
                             + ", ".join(other_references) + " are in the downloadable tables only.")
        if not figure_matrix.empty:
            matrix_for_figures = figure_matrix
            counts = (figure_matrix[figure_matrix["Status"].eq("Variant observed")]
                      .groupby("Variant_Key")["Sample_ID"].nunique().sort_values(ascending=False))
            repeated = [key for key, value in counts.items() if value >= 2][:TRAJECTORY_VARIANTS]
            heatmap_keys = sorted(list(counts.index[:HEATMAP_VARIANTS]),
                                  key=lambda key: (int(key.split(":")[1]), key))
            if repeated:
                add_figure(result, "trajectories", "Variant frequencies over time",
                           "Lines join reported frequencies only. Unreported observations sit in "
                           "the gray strip by reason; they are not plotted as 0%.",
                           lambda: figures.trajectory_figure(matrix_for_figures, repeated))
            else:
                note("temporal", "No variant was reported in two or more included pooled samples.")
            add_figure(result, "heatmap", "Sample-by-variant heatmap",
                       f"Up to {HEATMAP_VARIANTS} variants reported in the most samples.",
                       lambda: figures.heatmap_figure(matrix_for_figures, heatmap_keys))
            add_figure(result, "observation_status", "Missing and unknown observations",
                       "Share of the analyzed variant positions with a reported frequency, adequate "
                       "depth without a report, low or no coverage, unknown coverage, or a QC-excluded "
                       "call. Unknown is never treated as 0%.",
                       lambda: figures.observation_status_figure(matrix_for_figures))
            if trends["Change_In_Reported_Frequency"].notna().any():
                add_figure(result, "frequency_changes", "Changes between collection dates",
                           "First versus last reported frequency for each variant reported at least "
                           "twice at one site, with the number of days between them. Intervals are "
                           "irregular, so the changes are not rates.",
                           lambda: figures.frequency_change_figure(trends))
            else:
                note("temporal", "No variant has two reported frequencies at the same site, so "
                                 "changes between dates cannot be shown.")

        population = population_structure(samples, variants, pooled["position_table"])
        result["population"] = population
        note("population", population["message"])
        if population["available"]:
            result["tables"]["pcoa_coordinates.csv"] = population["coordinates"]
            result["tables"]["dissimilarity_matrix.csv"] = population["distances"].reset_index(
                names="Sample_ID")
            if population.get("diversity") is not None:
                result["tables"]["nucleotide_diversity.csv"] = population["diversity"]
            exploratory = population["method"] == EXPLORATORY_ORDINATION
            pcoa = population["pcoa"]
            if pcoa["negative_eigenvalues"]:
                note("population", f"{pcoa['negative_eigenvalues']} negative eigenvalue(s) "
                                   f"(magnitude {pcoa['negative_magnitude_ratio']:.1%} of the positive "
                                   "sum) show the distances are not exactly Euclidean; axes use "
                                   "positive eigenvalues only and no correction was applied.")
            add_figure(result, "pcoa", "Population structure (PCoA)",
                       population["message"],
                       lambda: figures.pcoa_figure(population["coordinates"], pcoa,
                                                   population["method"], exploratory))

    # 5. Cross-source comparison.
    if individual is not None and pooled is not None:
        cells = []
        comparison = compare_evidence(individual, pooled, reference, registry,
                                      settings["min_variant_depth"], settings["date_window_days"],
                                      cells=cells)
        counts = category_counts(comparison)
        result["comparison"] = {"table": comparison, "counts": counts}
        result["tables"]["cross_source_comparison.csv"] = comparison
        result["tables"]["cross_source_evidence_by_sample.csv"] = pd.DataFrame(
            cells, columns=["Variant_Key", "Sample_ID", "Sample_Type", "Status", "Alt_Frequency"])
        add_figure(result, "cross_source", "Cross-source evidence categories",
                   "Each reference-based variant is counted once.",
                   lambda: figures.cross_source_figure(counts))
        matrix_keys = evidence_matrix_keys(comparison)
        if matrix_keys and cells:
            individual_ids = sorted({cell["Sample_ID"] for cell in cells
                                     if cell["Sample_Type"] == INDIVIDUAL})
            pooled_rows = pooled["samples"][pooled["samples"]["Sample_ID"].isin(
                {cell["Sample_ID"] for cell in cells if cell["Sample_Type"] == POOLED})]
            pooled_rows = pooled_rows.sort_values(["Collection_Date", "Sample_ID"])
            pooled_labels = [f"{row['Sample_ID']} ({row['Collection_Date'].date()})"
                             if pd.notna(row["Collection_Date"]) else str(row["Sample_ID"])
                             for row in pooled_rows.to_dict("records")]
            add_figure(result, "evidence_matrix", "Variant-overlap matrix across both sources",
                       f"Up to {MATRIX_VARIANTS} comparable variants (shared first). Each column is "
                       "one sample; pooled samples are ordered by collection date. Detection, "
                       "non-detection with adequate coverage, and unknown coverage are drawn "
                       "differently.",
                       lambda: figures.evidence_matrix_figure(
                           cells, matrix_keys, individual_ids, list(pooled_rows["Sample_ID"]),
                           pooled_labels))
            add_figure(result, "cross_source_positions", "Reference positions by evidence category",
                       "Every comparable variant at its reference position, one row per category.",
                       lambda: figures.cross_source_positions_figure(comparison, reference["Length"]))
        if not comparison.empty:
            add_figure(result, "cross_source_timeline", "When variants were seen in each source",
                       "Triangles: individual sequences carrying the variant. Circles: pooled "
                       "samples reporting it (larger = higher frequency).",
                       lambda: figures.shared_timeline_figure(
                           comparison, individual["variants"],
                           pooled["variants"][pooled["variants"]["Analysis_Included"]]))
        if signature_path:
            signatures = load_signatures(signature_path)
            result["tables"]["variant_signature_overlap.csv"] = compare_signatures(
                signatures, reference, registry, individual, pooled, settings["min_variant_depth"])
    elif individual is not None or pooled is not None:
        note("cross_source", "Cross-source comparison needs individual sequences with a "
                             "reference FASTA and a pooled SNV table on a compatible reference.")

    if individual is not None and individual["calls"] is not None:
        add_figure(result, "callable_positions", "Callable and uncallable reference positions",
                   "How each reference position was called in each individual sequence after "
                   "alignment to the reference (N, ambiguity codes, deletions, and unsequenced ends "
                   "are never treated as reference bases).",
                   lambda: figures.callable_positions_figure(
                       individual["samples"], individual["calls"], reference["Length"]))

    if pooled is not None or individual is not None:
        add_figure(result, "quality", "Sequencing quality by source",
                   "Pooled genome coverage at 10x (paper threshold >90%) and individual callable "
                   "reference positions.",
                   lambda: figures.quality_figure(
                       pooled["samples"] if pooled else None,
                       individual["samples"] if individual else None,
                       settings["min_genome_coverage_percent"]))

    geography = geographic_summary(individual_data, pooled["samples"] if pooled else None)
    result["tables"]["geographic_context.csv"] = geography
    note("geography", "Locations are context, not proof of transmission. State-only names and "
                      "sewershed site labels are never treated as precise collection points, and "
                      "pooled site labels are not geocoded.")
    located = geography_samples(individual_data, pooled["samples"] if pooled else None)
    if not located.empty:
        result["tables"]["sample_locations.csv"] = located
        add_figure(result, "geography", "Geographic precision",
                   "Exact source-provided coordinates, approximate GeoNames place coordinates, broad "
                   "administrative areas, named sites without coordinates, and missing locations. "
                   "The coordinate panel is a plain latitude/longitude plot, not a map; broad areas "
                   "are never drawn as points.",
                   lambda: figures.geography_figure(located))
        if located["Collection_Date"].notna().any():
            add_figure(result, "sampling_timeline", "Sampling dates and sites",
                       "Every sample at its collection date, grouped by site, for both sources. "
                       "Gaps are real gaps in sampling; nothing is interpolated.",
                       lambda: figures.sampling_timeline_figure(located))
        else:
            note("temporal", "No collection dates were supplied, so no sampling timeline is drawn.")
    result["completeness"] = completeness_rows(individual_data, individual, pooled)
    if result["completeness"]:
        add_figure(result, "completeness", "Metadata completeness",
                   "Share of samples in each source that carry each kind of information. Fields "
                   "the uploaded files cannot provide are labeled instead of drawn as 0%.",
                   lambda: figures.completeness_figure(result["completeness"]))
    if relationship is not None and relationship.get("pairs") is not None:
        relationship_figures(result, relationship, note)

    result["summary"] = build_summary(result, individual_data)
    result["analyses"] = analysis_status(result, individual_data, relationship)
    result["integration"] = integration_summary(result, individual_data)
    result["tables"]["analysis_status.csv"] = pd.DataFrame(
        result["analyses"], columns=["Analysis", "Status", "Detail"])
    result["tables"]["provenance.csv"] = provenance(result, settings, labels, mafft_path)
    return result


def evidence_matrix_keys(comparison, limit=None):
    """Comparable variants for the evidence matrix: shared first, then by position."""

    limit = limit or MATRIX_VARIANTS
    priority = {category: index for index, category in enumerate(CATEGORY_ORDER)}
    data = comparison[~comparison["Category"].eq(INCOMPARABLE)].copy()
    if data.empty:
        return []
    data["Priority"] = data["Category"].map(priority)
    chosen = data.sort_values(["Priority", "Position"]).head(limit)
    return list(chosen.sort_values(["Position", "Variant_Key"])["Variant_Key"])


def precision_category(row):
    """Sort one sample's location into a precision class without upgrading it."""

    location = row.get("Location_Name")
    location = "" if location is None or pd.isna(location) else str(location).strip()
    precision = str(row.get("Geographic_Precision") or "")
    source = str(row.get("Coordinate_Source") or "")
    has_point = pd.notna(row.get("Latitude")) and pd.notna(row.get("Longitude"))
    if "too coarse" in precision or is_state_only_location(location):
        return "Broad administrative area"
    if has_point and source == "GeoNames":
        return "Approximate (GeoNames place)"
    if has_point:
        return "Source-provided coordinates"
    if location:
        return "Named place or site, no coordinates"
    return "No location"


def geography_samples(individual_data, pooled_samples):
    """One row per sample with its location precision class and dates."""

    columns = ["Sample_ID", "Sample_Type", "Collection_Date", "Location_Name", "Latitude",
               "Longitude", "Coordinate_Source", "Geographic_Precision", "Precision_Category",
               "Analysis_Included"]
    frames = []
    if individual_data is not None:
        frame = individual_data.copy()
        frame["Sample_Type"] = INDIVIDUAL
        frame["Analysis_Included"] = True
        frames.append(frame)
    if pooled_samples is not None:
        frames.append(pooled_samples.copy())
    if not frames:
        return pd.DataFrame(columns=columns)
    rows = []
    for frame in frames:
        for row in frame.to_dict("records"):
            rows.append({column: row.get(column, pd.NA) for column in columns})
            rows[-1]["Precision_Category"] = precision_category(row)
    table = pd.DataFrame(rows, columns=columns)
    table["Collection_Date"] = pd.to_datetime(table["Collection_Date"], errors="coerce")
    table["Latitude"] = pd.to_numeric(table["Latitude"], errors="coerce")
    table["Longitude"] = pd.to_numeric(table["Longitude"], errors="coerce")
    return table


def completeness_rows(individual_data, individual, pooled):
    """Counts of samples per source with dates, locations, coordinates, and coverage."""

    rows = []
    if individual_data is not None:
        categories = [precision_category(row) for row in individual_data.to_dict("records")]
        location = individual_data["Location_Name"].astype("string").str.strip().fillna("")
        rows.append({
            "Source": "Individual sequences", "Samples": len(individual_data),
            "With_Collection_Date": int(individual_data["Collection_Date"].notna().sum()),
            "With_Location": int(location.ne("").sum()),
            "With_Coordinates": sum(category in {"Source-provided coordinates",
                                                 "Approximate (GeoNames place)"}
                                    for category in categories),
            "With_Genome_Coverage": (int(individual["samples"]["Genome_Coverage_Percent"].notna().sum())
                                     if individual is not None else "needs a reference FASTA"),
        })
    if pooled is not None:
        samples = pooled["samples"]
        location = samples["Location_Name"].astype("string").str.strip().fillna("")
        coverage = int(samples["Genome_Coverage_Percent"].notna().sum())
        rows.append({
            "Source": "Pooled samples", "Samples": len(samples),
            "With_Collection_Date": int(samples["Collection_Date"].notna().sum()),
            "With_Location": int(location.ne("").sum()),
            "With_Coordinates": (int((samples["Latitude"].notna() & samples["Longitude"].notna()).sum())
                                 if pooled.get("metadata_supplied") else "needs pooled metadata"),
            "With_Genome_Coverage": coverage if coverage else "not supplied",
        })
    return rows


def relationship_figures(result, relationship, note):
    """Figures for the original Sprint 2 pairwise screening (values not recalculated)."""

    pairs = relationship["pairs"]
    thresholds = relationship["thresholds"]
    result["tables"]["pairwise_distances.csv"] = pairs
    if pairs.empty:
        note("relationships", "Fewer than two comparable sequences: no pairwise distances.")
        return
    add_figure(result, "genetic_distances", "Genetic-distance distribution",
               "Left: p-distance of every compared pair (dashed line: your cutoff). Right: genetic "
               "versus temporal distance; the shaded box meets both the genetic and date cutoffs, "
               "and the marker shows whether geography could be checked.",
               lambda: figures.distance_figure(pairs, thresholds))
    if relationship.get("clusters") is None or relationship.get("candidates") is None:
        return  # no screening lists supplied, so no network is drawn
    strict = []
    for cluster in relationship["clusters"]:
        members = set(cluster["samples"])
        strict += [(row["Sample_ID_1"], row["Sample_ID_2"]) for row in pairs.to_dict("records")
                   if row["Sample_ID_1"] in members and row["Sample_ID_2"] in members
                   and passes_all(row, thresholds)]
    candidates = [(pair["sample_1"], pair["sample_2"]) for pair in relationship["candidates"]]
    samples = sorted(set(pairs["Sample_ID_1"].astype(str)) | set(pairs["Sample_ID_2"].astype(str)))
    if len(samples) > 40:
        note("relationships", "The relationship network is drawn only for 40 or fewer sequences.")
        return
    add_figure(result, "relationship_network", "Sample-relationship network",
               "Each point is an individual sequence. Solid lines join pairs that meet the genetic, "
               "date, and distance cutoffs; dashed lines join pairs that meet the genetic and date "
               "cutoffs but whose geography could not be checked. Screening signals only.",
               lambda: figures.relationship_network_figure(samples, strict, candidates))


def passes_all(row, thresholds):
    """Mirror of the Sprint 2 strict rule, used only to draw lines for existing clusters."""

    max_p, max_days, max_km = thresholds
    if pd.isna(row["P_Distance"]) or pd.isna(row["Days_Apart"]):
        return False
    if row["P_Distance"] > max_p or row["Days_Apart"] > max_days:
        return False
    if pd.notna(row["Km_Apart"]):
        return row["Km_Apart"] <= max_km
    return str(row.get("Same_Location")).strip().casefold() == "true"


def analysis_status(result, individual_data, relationship):
    """What actually ran in this analysis, and why anything was skipped."""

    rows = []

    def add(name, status, detail):
        rows.append({"Analysis": name, "Status": status, "Detail": detail})

    pairs = relationship.get("pairs") if relationship else None
    if pairs is not None:
        add("Sequence comparison (MAFFT, p-distance)", "Completed", f"{len(pairs)} sequence pairs compared")
        dated = int(pd.to_numeric(pairs["Days_Apart"], errors="coerce").notna().sum())
        add("Temporal screening of sequence pairs", "Completed" if dated else "Skipped",
            f"{dated} of {len(pairs)} pairs have two collection dates")
        km = int(pd.to_numeric(pairs["Km_Apart"], errors="coerce").notna().sum())
        add("Geographic-distance screening", "Completed" if km else "Limited",
            f"{km} of {len(pairs)} pairs have two sufficiently precise coordinates")
    elif individual_data is None:
        add("Sequence comparison (MAFFT, p-distance)", "Skipped", "No individual sequences supplied")
    individual = result["individual"]
    if individual is not None:
        mapped = int(individual["samples"]["Variant_Profile_Status"].eq("Mapped to reference").sum())
        add("Reference-coordinate variant profiles", "Completed",
            f"{mapped} of {len(individual['samples'])} sequences mapped to the reference")
    elif individual_data is not None:
        add("Reference-coordinate variant profiles", "Skipped", "No reference FASTA supplied")
    pooled = result["pooled"]
    if pooled is not None:
        add("Pooled variant QC and frequencies", "Completed",
            f"{int(pooled['variants']['Analysis_Included'].sum())} of {len(pooled['variants'])} "
            "records retained")
        coverage = int(pooled["samples"]["Genome_Coverage_Percent"].notna().sum())
        add("Paper coverage filter (>90% at 10x)", "Completed" if coverage else "Not evaluable",
            f"Genome coverage known for {coverage} of {len(pooled['samples'])} pooled samples")
        add("Depth-versus-SNV regression (Fig 2B analogue)",
            "Completed" if result.get("regression") else "Skipped",
            "Needs mean depth for at least three included pooled samples"
            if not result.get("regression") else f"n = {result['regression']['n']}")
        add("Temporal frequency trends", "Completed" if "trajectories" in result["figures"] else "Limited",
            "Trajectories need a variant reported in two or more included samples"
            if "trajectories" not in result["figures"] else "Reported frequencies over collection dates")
        population = result["population"]
        if population is not None:
            add("Population structure", "Completed" if population["available"] else "Skipped",
                population["method"] if population["available"] else population["message"])
        add("Coverage-aware non-detection", "Completed" if pooled["position_table"] is not None else "Skipped",
            "Per-position depth supplied" if pooled["position_table"] is not None else
            "No per-position depth: unreported variants stay 'coverage unknown' (never 0%)")
    else:
        add("Pooled variant QC and frequencies", "Skipped", "No pooled variant file supplied")
    if result["comparison"] is not None:
        add("Cross-source comparison", "Completed",
            f"{len(result['comparison']['table'])} variants placed in evidence categories")
    else:
        add("Cross-source comparison", "Skipped",
            "Needs individual sequences, a reference FASTA, and pooled variants on the same reference")
    return rows


def integration_summary(result, individual_data):
    """Plain facts about how the two sources were connected in this analysis."""

    reference = result["reference"]
    individual = result["individual"]
    pooled = result["pooled"]
    facts = []
    if reference is not None:
        facts.append(("Shared reference", f"{reference['Reference_ID']}, {reference['Length']} bases, "
                                          f"SHA-256 {reference['Checksum'][:16]}..."))
    if individual is not None:
        mapped = individual["samples"]["Variant_Profile_Status"].eq("Mapped to reference")
        facts.append(("Individual sequences on the reference",
                      f"{int(mapped.sum())} of {len(mapped)} aligned to the reference with MAFFT "
                      f"(at least {MINIMUM_MAPPING_IDENTITY:.0%} identity over "
                      f"{MINIMUM_CALLABLE_POSITIONS} or more callable positions)"))
    if pooled is not None and reference is not None:
        variants = pooled["variants"]
        same_label = variants["Reference_ID"].eq(reference["Reference_ID"])
        checks = variants["QC_Reference_Base"].astype(str)
        facts.append(("Pooled records on the reference",
                      f"{int(same_label.sum())} of {len(variants)} name {reference['Reference_ID']}; "
                      f"Ref base agrees with the reference for {int(checks.eq('Pass').sum())} and "
                      f"disagrees for {int(checks.str.startswith('Fail').sum())}"))
        other = sorted(set(variants.loc[~same_label, "Reference_ID"]))
        if other:
            facts.append(("Kept apart", "Records on " + ", ".join(other) + " are incomparable: a "
                                        "different label is never assumed to be the same genome"))
    if result["comparison"] is not None:
        counts = result["comparison"]["counts"].set_index("Category")["Variants"]
        facts.append(("Variants by evidence", "; ".join(f"{name.split(':')[0]} {int(value)}"
                                                       for name, value in counts.items() if value)))
    date_ranges = []
    if individual_data is not None and individual_data["Collection_Date"].notna().any():
        dates = individual_data["Collection_Date"].dropna()
        date_ranges.append(("Individual", dates.min(), dates.max()))
    if pooled is not None and pooled["samples"]["Collection_Date"].notna().any():
        dates = pooled["samples"]["Collection_Date"].dropna()
        date_ranges.append(("Pooled", dates.min(), dates.max()))
    if date_ranges:
        text = "; ".join(f"{name} {first.date()} to {last.date()}" for name, first, last in date_ranges)
        if len(date_ranges) == 2:
            overlap = (min(date_ranges[0][2], date_ranges[1][2])
                       - max(date_ranges[0][1], date_ranges[1][1])).days
            text += (f"; the periods overlap for {overlap + 1} days" if overlap >= 0
                     else f"; the periods do not overlap ({-overlap} days apart)")
        facts.append(("Collection dates", text))
    return facts


def geographic_summary(individual_data, pooled_samples):
    """Group samples by documented location with their geographic precision."""

    rows = []
    sources = []
    if individual_data is not None:
        sources.append(("Individual sequence", individual_data))
    if pooled_samples is not None:
        sources.append(("Pooled sample", pooled_samples))
    for sample_type, data in sources:
        for location, group in data.groupby("Location_Name", dropna=False, sort=True):
            precision = group.get("Geographic_Precision", pd.Series(dtype=str)).astype(str)
            has_point = group["Latitude"].notna() & group["Longitude"].notna()
            usable = has_point & ~precision.str.contains("too coarse") & ~precision.str.contains("no coordinates")
            dates = group["Collection_Date"].dropna()
            rows.append({
                "Location_Name": location if pd.notna(location) and str(location) else "(no location)",
                "Sample_Type": sample_type,
                "Samples": len(group),
                "Geographic_Precision": "; ".join(sorted(set(precision))),
                "Samples_With_Coordinates": int(has_point.sum()),
                "Usable_For_Kilometer_Screening": int(usable.sum()),
                "First_Date": dates.min() if not dates.empty else pd.NaT,
                "Last_Date": dates.max() if not dates.empty else pd.NaT,
            })
    return pd.DataFrame(rows, columns=[
        "Location_Name", "Sample_Type", "Samples", "Geographic_Precision",
        "Samples_With_Coordinates", "Usable_For_Kilometer_Screening", "First_Date", "Last_Date",
    ])


def build_summary(result, individual_data):
    summary = {}
    if individual_data is not None:
        summary["individual_samples"] = len(individual_data)
    individual = result["individual"]
    if individual is not None:
        samples = individual["samples"]
        summary["individual_mapped"] = int(samples["Variant_Profile_Status"].eq("Mapped to reference").sum())
        summary["individual_variants"] = int(individual["variants"]["Variant_Key"].nunique())
    pooled = result["pooled"]
    if pooled is not None:
        summary["pooled_samples"] = len(pooled["samples"])
        summary["pooled_samples_included"] = int(pooled["samples"]["Analysis_Included"].sum())
        summary["pooled_records"] = len(pooled["variants"])
        summary["pooled_records_included"] = int(pooled["variants"]["Analysis_Included"].sum())
        summary["pooled_variants"] = int(
            pooled["variants"].loc[pooled["variants"]["Analysis_Included"], "Variant_Key"].nunique())
        dates = pd.concat([pooled["samples"]["Collection_Date"],
                           individual["samples"]["Collection_Date"] if individual else pd.Series(dtype="datetime64[ns]")])
        dates = dates.dropna()
        if not dates.empty:
            summary["date_start"] = dates.min().date().isoformat()
            summary["date_end"] = dates.max().date().isoformat()
        summary["pooled_coverage_known"] = int(pooled["samples"]["Genome_Coverage_Percent"].notna().sum())
    if result["reference"] is not None:
        summary["reference"] = (f"{result['reference']['Reference_ID']} "
                                f"({result['reference']['Length']} bases)")
    if result["comparison"] is not None:
        counts = result["comparison"]["counts"].set_index("Category")["Variants"]
        summary["category_counts"] = counts.to_dict()
        summary["shared_variants"] = int(counts.iloc[0])
    population = result["population"]
    if population is not None:
        summary["population_method"] = population["method"] if population["available"] else "Not available"
    return summary


def provenance(result, settings, labels, mafft_path):
    """A record of data sources, versions, parameters, and reproduction status."""

    rows = []

    def add(category, field, value):
        rows.append({"Category": category, "Field": field, "Value": value})

    for key, name in labels.items():
        add("Data source", key, name)
    if result["pooled"] is not None:
        add("Input format", "Pooled SNV table", result["pooled"]["input_format"])
        table = result["pooled"]["position_table"]
        add("Input format", "Per-position table",
            "Not supplied" if table is None else
            ("Depth with A/C/G/T counts" if table.attrs.get("has_base_counts") else "Depth only"))
    if result["individual"] is not None or labels.get("fasta"):
        add("Input format", "Individual sequences", "FASTA + metadata (consensus base calls)")
    reference = result["reference"]
    if reference is not None:
        add("Reference", "Identifier", reference["Reference_ID"])
        add("Reference", "Length", reference["Length"])
        add("Reference", "SHA-256", reference["Checksum"])
        add("Reference", "Source file", reference["Source_File"])
    else:
        add("Reference", "Identifier", "No reference FASTA supplied")

    add("Analysis mode", "Quality-control mode", settings["qc_mode"])
    add("Software", "Python", platform.python_version())
    for package in ["pandas", "numpy", "matplotlib", "biopython", "Flask"]:
        add("Software", package, package_version(package))
    if mafft_path:
        from system_diagnostics import read_mafft_version
        add("Software", "MAFFT", read_mafft_version(mafft_path) or f"version not reported ({mafft_path})")
    for name in ["min_variant_depth", "coverage_depth", "min_genome_coverage_percent",
                 "max_fisher_p", "min_base_quality", "date_window_days"]:
        add("Parameter", name, settings[name])
    add("Parameter", "minimum_callable_positions_for_mapping", MINIMUM_CALLABLE_POSITIONS)
    add("Parameter", "minimum_reference_identity_for_mapping", MINIMUM_MAPPING_IDENTITY)

    if result["pooled"] is not None:
        variants = result["pooled"]["variants"]
        samples = result["pooled"]["samples"]
        add("Exclusions", "Pooled records excluded", int((~variants["Analysis_Included"]).sum()))
        add("Exclusions", "Pooled samples excluded", int((~samples["Analysis_Included"]).sum()))
        for row in result["pooled"]["report"].to_dict("records"):
            add("Paper method", row["Step"], f"{row['Status']}: {row['Detail']}")
    if result["individual"] is not None:
        samples = result["individual"]["samples"]
        add("Exclusions", "Individual samples without a reference mapping",
            int((~samples["Variant_Profile_Status"].eq("Mapped to reference")).sum()))
        add("Exclusions", "Individual insertion/deletion events (not compared)",
            len(result["individual"]["indels"]))
    population = result["population"]
    if population is not None:
        add("Paper method", "Population structure",
            population["method"] if population["available"] else population["message"])
    return pd.DataFrame(rows, columns=["Category", "Field", "Value"])
