"""What one finished analysis can show: the single source of truth for the results page and report.

The results page, the PDF report, and the tests all ask this module what ran, so
the rules are written once. Every answer is read from the *outputs* of the run
(tables, figures, and the analysis-status record), never from the analysis mode
alone: a pooled upload may support allele-frequency analysis but not the
paper's PCoA, and an individual upload without a reference has no
reference-position variants.

Nothing here calculates science. It counts rows that already exist, decides
which figures are worth showing first, and writes short plain-language labels.
The result must stay JSON-serializable because the PDF report reads it back
from ``analysis_summary.json``.
"""

COVERAGE_MISSING = ("Genome-wide coverage was not supplied, so coverage percentages could not be calculated. "
                    "Read depth at the reported variants is still shown.")
NO_COORDINATES = ("Every located sample has a place or site name but no coordinates, so no map or "
                  "coordinate plot is drawn and no kilometer distances are calculated.")
NO_LOCATIONS = "No sample has a collection location, so no geographic chart is drawn."
PRECISE = ["Source-provided coordinates", "Approximate (GeoNames place)"]

# Figures shown first in each category (the rest stay one click away under "Additional figures").
PRIMARY = {
    "cross_source": "cross-source", "genetic_distances": "relationships", "relationship_network": "relationships",
    "genome_positions": "variation", "mutation_annotation_overview": "variation",
    "vaf_time_series": "temporal", "sampling_timeline": "temporal", "geography": "temporal",
    "depth_distribution": "quality", "vaf_depth": "quality", "quality": "quality", "pcoa": "advanced",
}
# Where every other figure lives (secondary).
SECONDARY = {
    "completeness": "overview", "evidence_matrix": "cross-source", "cross_source_positions": "cross-source",
    "cross_source_timeline": "cross-source", "trajectories": "temporal", "frequency_changes": "temporal",
    "observation_status": "temporal", "heatmap": "temporal", "qc_outcomes": "quality",
    "variants_per_sample": "quality", "depth_vs_snvs": "quality", "callable_positions": "quality",
    "depth_dependence": "quality", "uncertainty_models": "advanced", "replicate_concordance": "advanced",
    "mutation_profile_heatmap": "advanced", "mutation_group_heatmap": "advanced",
    "mutation_group_trajectories": "advanced", "grouping_comparison": "advanced",
}
# Order of the main figures in the PDF report (only those that are shown first on the website).
REPORT_ORDER = ["cross_source", "genetic_distances", "genome_positions", "mutation_annotation_overview",
                "vaf_time_series", "trajectories", "sampling_timeline", "geography", "depth_distribution",
                "vaf_depth", "quality", "pcoa", "relationship_network"]
# Analyses that belong to one evidence source; hidden from the short "unavailable" list when that
# source was simply not part of the run (they stay in analysis_status.csv and the provenance table).
INDIVIDUAL_ONLY = ("Sequence comparison", "Temporal screening of sequence pairs", "Geographic-distance screening",
                   "Reference-coordinate variant profiles")
POOLED_ONLY = ("Pooled variant QC", "Paper coverage filter", "Depth-versus-SNV regression",
               "Temporal frequency trends", "Pooled reference verification", "Population structure",
               "Depth-aware VAF", "Overdispersion", "Exploratory mutation-frequency groups",
               "Mutation-frequency time series", "Coverage-aware non-detection")


def count_text(value, singular, plural=None):
    return f"{value:,} {singular if value == 1 else (plural or singular + 's')}"


def mode_and_title(result):
    individual, pooled = result.get("individual") is not None, result.get("pooled") is not None
    relationship = "pairwise_distances.csv" in result["tables"]
    sprint5 = " + 5" if result.get("depth_aware") is not None else ""
    if pooled and (individual or relationship):
        if result.get("comparison") is not None:
            return "combined", f"Combined Analysis — Sprint 2 + 3{sprint5}"
        return "separate", f"Individual and Pooled Sequencing, analyzed separately — Sprint 2 + 3{sprint5}"
    if pooled:
        return "pooled", f"Pooled Pathogen Sequencing — Sprint 3{sprint5}"
    return "individual", "Individual Pathogen Sequences — Sprint 2"


def location_counts(result):
    located = result["tables"].get("sample_locations.csv")
    if located is None or located.empty:
        return {}
    return {name: int(count) for name, count in located["Precision_Category"].value_counts().items()}


def figure_reasons(result, available):
    """Figures that were drawn but would mislead as a main chart: {name: reason}. Their PNGs stay downloadable."""

    reasons = {}
    figures = result["figures"]
    if "quality" in figures and not available["coverage"] and result.get("individual") is None:
        reasons["quality"] = COVERAGE_MISSING
    if "geography" in figures and not available["coordinates"]:
        reasons["geography"] = NO_COORDINATES if available["locations"] else NO_LOCATIONS
    return reasons


def figure_roles(result, available, suppressed):
    """{figure name: (category, 'primary' or 'secondary')} for every figure that is shown."""

    roles = {}
    for name in result["figures"]:
        if name in suppressed:
            continue
        if name in PRIMARY:
            roles[name] = (PRIMARY[name], "primary")
        else:
            roles[name] = (SECONDARY.get(name, "quality"), "secondary")
    # One frequency trend first: the uncertainty-aware series when it exists, else the plain trajectories.
    if "trajectories" in roles and "vaf_time_series" not in roles:
        roles["trajectories"] = ("temporal", "primary")
    # With pooled data, the coverage bars lead only when coverage was measured (individual callable
    # positions alone, as in an individual-only run, stay a main chart).
    if "quality" in roles and not available["coverage"] and available["pooled"]:
        roles["quality"] = ("quality", "secondary")
    return roles


def key_metrics(result, available):
    """A few headline facts, each written from a calculated value."""

    summary = result["summary"]
    metrics = []
    if "individual_samples" in summary:
        mapped = (f"; {summary['individual_mapped']:,} placed on the reference genome"
                  if "individual_mapped" in summary and available["individual_variants"] else "")
        metrics.append(count_text(summary["individual_samples"], "individual sequence") + " analyzed" + mapped)
    if available["relationships"]:
        pairs = len(result["tables"]["pairwise_distances.csv"])
        metrics.append(count_text(pairs, "sequence pair") + " compared after MAFFT alignment")
    if "pooled_samples" in summary:
        metrics.append(count_text(summary["pooled_samples"], "pooled sample") + " analyzed "
                       f"({summary['pooled_samples_included']:,} kept after quality control)")
        metrics.append(count_text(summary["pooled_records"], "reported SNV observation") + " "
                       f"({summary['pooled_records_included']:,} variant records kept after quality control)")
    if summary.get("reference"):
        status = available.get("reference_status")
        metrics.append(f"Reference genome: {summary['reference']}" + (f" ({status.lower()})" if status else ""))
    elif result.get("pooled") is not None and result["pooled"].get("reference_status") is not None:
        checks = result["pooled"]["reference_status"]
        by_status = checks.groupby("Status", sort=False)["Reference_ID"].agg(lambda ids: ", ".join(map(str, ids)))
        metrics.append("Reference genome: " + "; ".join(f"{ids} ({str(status).lower()})"
                                                        for status, ids in by_status.items())
                       + " — no reference FASTA uploaded")
    if available["cross_source"]:
        metrics.append(count_text(summary.get("shared_variants", 0), "variant") + " observed in both sources")
    for flag, text in [("annotation", "Mutation annotations available"),
                       ("temporal", "Temporal patterns available"),
                       ("population", "Population structure calculated")]:
        if available[flag]:
            metrics.append(text + (" (exploratory method)" if flag == "population"
                                   and not available["population_paper_method"] else ""))
    if "date_start" in summary:
        metrics.append(f"Collection dates {summary['date_start']} to {summary['date_end']}")
    return metrics


def headline_limitations(result, available, mode):
    """The few limitations a reader must see first; each section keeps its full list."""

    items = []
    if mode != "individual" and not available["coverage"]:
        items.append("Complete genome-wide coverage unavailable: coverage percentages and the paper's >90% "
                     "coverage filter could not be checked.")
    if mode != "individual" and not available["position_depth"]:
        items.append("No per-position depth: a mutation missing from a pooled sample is unknown, not zero.")
    aware = result.get("depth_aware")
    if aware is not None:
        if not available["calibrated"]:
            items.append("Empirical overdispersion not calibrated (no technical replicates): intervals reflect "
                         "read sampling only; rho values are hypothetical scenarios.")
        if aware["reconstructed_counts"]:
            items.append("Alternative read counts were reconstructed from frequency x depth, so intervals are "
                         "approximate.")
    if available["grouping"]:
        items.append("Mutation-frequency grouping remains experimental and did not meet its validation criteria.")
    if available["population"] and not available["population_paper_method"]:
        items.append("Population structure uses an exploratory Jaccard PCoA, not the paper's Yue–Clayton method.")
    distance = next((row for row in result.get("analyses", []) if row["Analysis"] == "Geographic-distance screening"), None)
    if distance is not None and distance["Status"] != "Completed":
        items.append(f"Kilometer-distance screening was limited: {distance['Detail']}.")
    if individual_upload_failed(result):
        items.insert(0, "The uploaded individual sequences could not be compared, so no individual-sequence or "
                        "cross-source results are shown.")
    items.append("Results are screening signals; they do not confirm transmission, an outbreak, or a lineage.")
    return items


def individual_upload_failed(result):
    """Individual sequences were uploaded but produced no result (for example, MAFFT was unavailable)."""

    labels = result.get("input_labels") or {}
    return "fasta" in labels and result.get("individual") is None and "pairwise_distances.csv" not in result["tables"]


def unavailable_analyses(result, mode, suppressed):
    """Short list of analyses that were skipped or limited, with the recorded reason."""

    rows = []
    for row in result.get("analyses", []):
        if row["Status"] in {"Completed"}:
            continue
        name = row["Analysis"]
        if mode == "pooled" and name.startswith(INDIVIDUAL_ONLY):
            continue
        if mode == "individual" and name.startswith(POOLED_ONLY):
            continue
        if mode in {"pooled", "individual"} and name.startswith("Cross-source"):
            continue  # needs both sources; covered by the "not part of this analysis" line below
        rows.append({"Analysis": name, "Status": row["Status"], "Reason": row["Detail"]})
    if mode == "pooled" and individual_upload_failed(result):
        rows.append({"Analysis": "Individual-sequence and cross-source analyses (Sprint 2)", "Status": "Could not run",
                     "Reason": "Individual sequences were uploaded, but the sequence comparison did not complete "
                               "(see the message at the top of the page); only the pooled results are shown."})
    elif mode == "pooled":
        rows.append({"Analysis": "Individual-sequence and cross-source analyses (Sprint 2)",
                     "Status": "Not part of this analysis", "Reason": "No individual sequences were uploaded."})
    elif mode == "individual":
        rows.append({"Analysis": "Pooled-sequencing and cross-source analyses (Sprint 3 + 5)",
                     "Status": "Not part of this analysis", "Reason": "No processed pooled SNV table was uploaded."})
    titles = {name: figure["title"] for name, figure in result["figures"].items()}
    for name, reason in suppressed.items():
        rows.append({"Analysis": titles.get(name, name), "Status": "Chart not shown", "Reason": reason})
    return rows


def denominators(result):
    """Each count the page uses, with what it counts (the numbers come from different processing stages)."""

    summary = result["summary"]
    rows = []

    def add(label, value, definition):
        if value is not None:
            rows.append({"label": label, "value": int(value), "definition": definition})

    if "pooled_records" in summary:
        add("SNV observations in the upload", summary["pooled_records"],
            "Every pooled SNV record as uploaded (one row per sample and mutation).")
        add("Variant records after QC", summary["pooled_records_included"],
            "Records that passed the quality-control preset; excluded records stay in the downloads with their reason.")
    aware = result.get("depth_aware")
    if aware is not None:
        repeated = aware["observations"] - aware["observations_counted"]
        add("Unique mutation observations after duplicate handling", aware["observations_counted"],
            f"Records with {count_text(repeated, 'identical repeat')} counted once; the read-sampling interval is "
            "calculated for each of them, including records excluded by QC (their QC status is kept).")
    if "pooled_variants" in summary:
        add("Distinct mutations after QC", summary["pooled_variants"],
            "Different position:Ref>Alt changes among the records that passed QC (Genomic Variation, grouping "
            "eligibility).")
    annotation = result.get("annotation")
    if annotation and annotation.get("available"):
        add("Distinct mutations evaluated for annotation", annotation["evaluated"],
            "Every distinct reported change from each source, including changes seen only in records excluded by "
            "QC, so each reported mutation has a predicted consequence.")
    if "individual_variants" in summary:
        add("Distinct mutations in individual sequences", summary["individual_variants"],
            "Different reference-position changes carried by the individual sequences placed on the reference.")
    return rows


def distinct_findings(findings):
    """Main findings without the sample counts and date range that the headline counts already state."""

    return [text for text in findings
            if " records were analyzed;" not in text and not text.startswith("Collection dates run from")]


def analysis_availability(result):
    """The availability record for one result (JSON-serializable)."""

    tables, figures, summary = result["tables"], result["figures"], result["summary"]
    pooled, individual = result.get("pooled"), result.get("individual")
    mode, title = mode_and_title(result)
    locations = location_counts(result)
    aware = result.get("depth_aware")
    population = result.get("population")
    annotation = result.get("annotation")
    grouping = (aware or {}).get("grouping") or {}
    reference_status = None  # the pooled check of the shared reference named in the summary
    if pooled is not None and pooled.get("reference_status") is not None and summary.get("reference"):
        main = str(summary["reference"]).split(" (")[0]
        rows = pooled["reference_status"][pooled["reference_status"]["Reference_ID"].astype(str).eq(main)]
        reference_status = str(rows["Status"].iloc[0]) if len(rows) else None
    available = {
        "individual": individual is not None or "pairwise_distances.csv" in tables,
        "individual_variants": individual is not None,
        "pooled": pooled is not None,
        "relationships": "pairwise_distances.csv" in tables,
        "cross_source": result.get("comparison") is not None,
        "variation": "variant_positions.csv" in tables,
        "annotation": bool(annotation and annotation.get("available")),
        "dates": "date_start" in summary,
        "temporal": bool(aware and aware.get("eligible_series")) or "trajectories" in figures,
        "locations": sum(count for name, count in locations.items() if name != "No location") > 0,
        "coordinates": any(locations.get(name, 0) for name in PRECISE),
        "coverage": bool(summary.get("pooled_coverage_known")),
        "position_depth": pooled is not None and pooled.get("position_table") is not None,
        "population": bool(population and population.get("available")),
        "population_paper_method": bool(population and population.get("available")
                                        and not str(population.get("method", "")).startswith("Exploratory")),
        "depth_aware": aware is not None,
        "calibrated": bool(aware and (aware.get("uncertainty") or {}).get("calibrated_specimens")),
        "grouping": bool(grouping) and grouping.get("enabled", True),
        "grouping_experimental": bool(grouping.get("experimental")),
        "reference_status": reference_status,
    }
    suppressed = figure_reasons(result, available)
    roles = figure_roles(result, available, suppressed)
    primary = [name for name in REPORT_ORDER if roles.get(name, ("", ""))[1] == "primary"]
    return {
        "mode": mode, "title": title,
        "status": "Completed with limitations" if any(row["Status"] != "Completed"
                                                      for row in result.get("analyses", [])) else "Completed",
        "available": available,
        "locations": locations,
        "key_metrics": key_metrics(result, available),
        "limitations": headline_limitations(result, available, mode),
        "unavailable": unavailable_analyses(result, mode, suppressed),
        "suppressed_figures": suppressed,
        "figure_roles": {name: list(role) for name, role in roles.items()},
        "report_figures": primary,
        "denominators": denominators(result),
    }
