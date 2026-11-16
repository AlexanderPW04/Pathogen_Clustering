"""Turn integrated analysis results into simple values for the dashboard page.

Nothing here calculates science; it only arranges the tables, figures, and
notes that `genomics.integrated_analysis` already produced into the eight
dashboard sections. Every section has the same parts: purpose, summary
statistics, figures with captions, tables, limitations, and downloads.
Integrated (both-source) sections come before source-specific ones.
"""

from base64 import b64encode

import pandas as pd

from dashboard_guide import FIGURE_GUIDE, interpretation
from genomics.report_exports import PDF_NAME, SUMMARY_NAME, XLSX_NAME, summary_json
from genomics.spreadsheet_safety import escape_csv_bytes, safe_csv


MISSING = "—"  # an em dash marks a value that is not available
NOT_AVAILABLE = "Not available"  # a statistic the inputs cannot provide (not zero)
WIDE_FIGURE_PIXELS = 1400
TABLE_ROWS = 25  # rows shown at first; up to EXPANDED_ROWS more behind "Show all"; CSVs keep everything
EXPANDED_ROWS = 200


def format_value(value, digits=4):
    if value is None:
        return MISSING
    if isinstance(value, pd.Timestamp):
        return value.strftime("%Y-%m-%d") if pd.notna(value) else MISSING
    if isinstance(value, (bool,)):
        return "Yes" if value else "No"
    try:
        if pd.isna(value):
            return MISSING
    except (TypeError, ValueError):
        pass
    if isinstance(value, float):
        return f"{value:.{digits}g}" if abs(value) < 1e-3 and value != 0 else f"{round(value, digits)}"
    return str(value)


def table_view(table, columns=None, limit=None, title="", filename=None, expand=0):
    """Return {'title', 'columns', 'rows', 'more_rows', 'total', 'shown', 'filename'}.

    `limit` rows are shown at first; up to `expand` further rows can be opened
    on the page. Downloads always contain every row at full precision.
    """

    if table is None:
        return None
    data = table if columns is None else table[[name for name in columns if name in table.columns]]
    total = len(data)
    first = data if limit is None else data.head(limit)
    more = data.iloc[len(first):len(first) + expand] if expand else data.iloc[0:0]
    rows = [[format_value(value) for value in row] for row in first.itertuples(index=False)]
    more_rows = [[format_value(value) for value in row] for row in more.itertuples(index=False)]
    return {"title": title, "columns": list(data.columns), "rows": rows, "more_rows": more_rows,
            "total": total, "shown": len(rows) + len(more_rows), "filename": filename}


def png_width(png):
    """Pixel width from a PNG header (bytes 16-19), or 0 if it is not a PNG."""

    return int.from_bytes(png[16:20], "big") if png[:8] == b"\x89PNG\r\n\x1a\n" else 0


def figure_view(result, name, download_id):
    figure = result["figures"].get(name)
    if figure is None:
        return None
    return {
        "name": name,
        "title": figure["title"],
        "caption": figure["caption"],
        "data_uri": "data:image/png;base64," + b64encode(figure["png"]).decode("ascii"),
        "filename": f"{name}.png",
        "download_id": download_id,
        # Wider than the page: shown at full size in a scrolling box so labels stay readable.
        "wide": png_width(figure["png"]) > WIDE_FIGURE_PIXELS,
        "guide": FIGURE_GUIDE.get(name),
        "interpretation": interpretation(result, name),
    }


def stat(label, value, detail=""):
    """One summary statistic; None becomes 'Not available', which is not zero."""

    return {"label": label, "value": NOT_AVAILABLE if value is None else format_value(value),
            "available": value is not None, "detail": detail}


def short_category(name):
    """'Unresolved: pooled coverage ...' -> 'Unresolved (pooled)' so the two differ."""

    label, _, detail = name.partition(":")
    if label == "Unresolved":
        return "Unresolved (pooled)" if "pooled" in detail else "Unresolved (individual)"
    return label


SECTION_SPRINTS = {
    "overview": "Sprint 2 + Sprint 3", "cross-source": "Sprint 2 + Sprint 3 (integration)",
    "variation": "Sprint 2 + Sprint 3", "temporal": "Sprint 3 (with Sprint 2 dates)",
    "geography": "Sprint 2 + Sprint 3", "quality": "Sprint 2 + Sprint 3",
    "population": "Sprint 3 (Fontenele et al. 2021)", "relationships": "Sprint 2",
    "depth-aware": "Sprint 5 (optional; builds on Sprint 3)",
}

DEPTH_AWARE_LIMITATIONS = [
    "Higher sequencing depth generally makes a frequency more precise; low depth means more sampling "
    "uncertainty. Depth alone does not show that a variant call is correct.",
    "The 95% Wilson interval reflects read sampling only (a binomial model). It does not include "
    "sequencing or amplification errors, primer bias, mapping errors, or differences between wastewater "
    "samples, so the real uncertainty is usually larger.",
    "A high frequency does not show that a mutation spreads more easily.",
    "Similar frequencies do not show that two mutations are on the same viral genome: short reads rarely "
    "connect distant positions.",
    "Observations without quality information (base quality, strand balance, caller test) remain limited.",
    "Changes and slopes are descriptive, not statistical tests. They do not identify lineages, outbreaks, "
    "or transmission. Mutation-frequency grouping (Phase 3) is exploratory and kept in its own panel.",
]

OVERDISPERSION_LIMITATIONS = [
    "Overdispersion (rho) is estimated only from technical replicates of the same specimen. Different "
    "dates, sites, or mutations are never used instead, because real biological differences would be "
    "mistaken for sequencing noise.",
    "Calibrated beta-binomial intervals cover read sampling plus library-to-library variation. They do "
    "not correct primer bias shared by every library, mapping or reference errors, contamination, or "
    "biological differences between samples, and one rho is assumed for every site of a library.",
    "Hypothetical rho values are sensitivity scenarios: they show how the uncertainty would change, not "
    "what it is. A calibration on synthetic data does not show that real wastewater libraries are calibrated.",
]


def format_rho(value):
    return "" if value is None or pd.isna(value) else f"{value:.3g}"


def uncertainty_view(summary):
    """The Phase 2 conclusion panel: one sentence, the model table, a few examples, the reasons."""

    rows = []
    for row in summary["comparison"].to_dict("records"):
        rows.append([row["Model"], format_rho(row["Rho"]) or row["Rho_Source"], row["Calibration"],
                     row["Observations_With_Interval"],
                     "" if pd.isna(row["Median_Width"]) else f"{row['Median_Width']:.3f}",
                     "" if pd.isna(row["Median_Width_Ratio_To_Wilson"])
                     else f"{row['Median_Width_Ratio_To_Wilson']:.2f}x"])
    return {"explanation": summary["explanation"], "headline": summary["headline"],
            "calibrated": bool(summary["calibrated_specimens"]),
            "comparison": {"columns": ["Model", "rho", "Calibration", "Observations with interval",
                                       "Median 95% interval width", "Width vs basic"], "rows": rows},
            "examples": summary["examples"], "reasons": summary["reasons"]}


GROUPING_FIGURES = ["mutation_group_heatmap", "mutation_profile_heatmap", "mutation_group_trajectories",
                    "grouping_comparison"]
GROUPING_DOWNLOADS = ["sprint5_mutation_groups.csv", "sprint5_mutation_group_summary.csv",
                      "sprint5_mutation_group_trajectories.csv", "sprint5_mutation_eligibility.csv",
                      "sprint5_mutation_profiles.csv", "sprint5_grouping_parameters.csv",
                      "sprint5_grouping_comparison.csv", "sprint5_group_individual_overlap.csv"]


def grouping_view(grouping, figures, tables, available):
    """The Phase 3 panel: what was evaluated, whether groups are shown, and why."""

    if not grouping:
        return None
    mode = grouping.get("mode", "profiles")
    badge = ("Off" if not grouping.get("enabled", True) else
             "Experimental (not validated)" if grouping.get("experimental") else "Profiles only: no groups shown")
    stats = []
    if grouping.get("enabled", True) and "units" in grouping:
        stats = [
            stat("Mutations evaluated", grouping["mutations"], f"across {grouping['units']} observation unit(s)"),
            stat("Eligible for grouping", grouping.get("eligible_ordinary", 0),
                 f"known in at least {grouping['minimum_shared']} units; {grouping['eligible']} with read counts "
                 "(needed by the depth-aware method)"),
            stat("Excluded because of missing information", grouping["excluded"],
                 "reasons in sprint5_mutation_eligibility.csv"),
            stat("Calibration status", grouping["calibration"].split(":")[0], grouping["calibration"]),
        ]
        ran = [item for item in grouping["methods"] if item.get("ran")]
        if grouping.get("experimental"):
            depth = next((item for item in ran if item["method"].startswith("Depth-aware")), None)
            stats.insert(3, stat("Exploratory groups found (experimental)", depth["groups"] if depth else "None",
                                 f"{depth['grouped']} of {depth['analyzed']} analyzed mutations grouped"
                                 if depth else "the depth-aware method did not run"))
            stats.insert(4, stat("Grouping method", "Depth-aware (MAGOS-inspired)" if depth else "None",
                                 "ordinary clustering shown only for comparison"))
    methods = {"columns": ["Method", "Measurement model", "Mutations analyzed", "Groups", "Mutations grouped",
                           "Note"],
               "rows": [[item["method"], item["model"], item["analyzed"], item["groups"], item["grouped"],
                         item.get("reason") or item.get("limit_note") or ""] for item in grouping["methods"]]}
    comparison = {"columns": ["Comparison", "Mutations compared", "Adjusted Rand index", "Groups (first)",
                              "Groups (second)", "Note"],
                  "rows": [[row["Comparison"], row["Mutations_Compared"], f"{row['Adjusted_Rand_Index']:.2f}",
                            row["Groups_First"], row["Groups_Second"], row["Note"]]
                           for row in grouping.get("comparison", [])]}
    summary_table = None
    if "sprint5_mutation_group_summary.csv" in tables and len(tables["sprint5_mutation_group_summary.csv"]):
        summary_table = table_view(tables["sprint5_mutation_group_summary.csv"], columns=[
            "Group_ID", "Method", "Size", "Position_Range", "Units_With_Two_Or_More_Members", "Median_Depth",
            "Homogeneity_P", "Group_VAF_Min", "Group_VAF_Max", "Group_Stability", "Members"],
            limit=12, title="Exploratory groups (experimental; not lineages)",
            filename="sprint5_mutation_group_summary.csv")
    overlap = None
    if "sprint5_group_individual_overlap.csv" in tables:
        overlap = table_view(tables["sprint5_group_individual_overlap.csv"], limit=12,
                             title="Individual sequences compared with each group (descriptive only)",
                             filename="sprint5_group_individual_overlap.csv")
    return {"title": grouping["title"], "explanation": grouping["explanation"], "badge": badge, "mode": mode,
            "experimental": bool(grouping.get("experimental")), "stats": stats, "reasons": grouping["reasons"],
            "validation": grouping.get("validation", ""), "methods": methods if methods["rows"] else None,
            "comparison": comparison if comparison["rows"] else None, "summary_table": summary_table,
            "overlap": overlap, "figures": figures, "limitations": grouping["limitations"],
            "downloads": [name for name in GROUPING_DOWNLOADS if name in tables
                          and (available is None or name in available)]}


ANNOTATION_DOWNLOADS = ["sprint5_annotated_mutations.csv", "sprint5_mutation_annotations.csv",
                        "sprint5_annotated_observations.csv", "sprint5_annotation_sources.csv",
                        "sprint5_annotation_comparison.csv"]
ANNOTATION_TABLE_HEADINGS = {
    "Variant_Key": "Mutation", "Genes": "Gene", "Mature_Peptides": "Protein region (mature peptide)",
    "Nucleotide_Change": "Nucleotide change", "Amino_Acid_Changes": "Predicted amino-acid change",
    "Consequence_Category": "Predicted consequence", "Reference_Check": "Reference verification",
    "Pooled_Samples_Observed": "Pooled samples reporting it", "Median_VAF": "Median VAF", "Max_VAF": "Highest VAF",
    "Max_VAF_Wilson_95": "95% Wilson interval at highest VAF", "Min_Depth": "Lowest depth",
    "Max_Depth": "Highest depth", "Count_Basis": "Read-count basis",
    "Series_With_Two_Or_More_Dates": "Site series with 2+ dates",
    "Individual_Sequences_Carrying": "Individual sequences carrying it",
    "Exploratory_Group": "Exploratory group (not validated)",
}


def annotation_view(annotation, figures, tables, available, messages=()):
    """The Phase 4 panel: predicted coding consequences with their sequencing evidence."""

    if not annotation:
        return None
    view = {"title": annotation["title"], "explanation": annotation["explanation"],
            "available": annotation["available"], "limitations": annotation["limitations"],
            "message": annotation.get("message", ""), "notes": list(messages)}
    if not annotation["available"]:
        return view
    other = ", ".join(f"{name.lower()} {count}" for name, count in annotation["other_detail"].items())
    view["stats"] = [
        stat("Distinct SNVs evaluated", annotation["evaluated"]),
        stat("Successfully annotated", annotation["annotated"],
             "reference identity and reference base verified; every overlapping coding region interpreted"),
        stat("Synonymous", annotation["synonymous"], "same amino acid; not necessarily without effect"),
        stat("Missense", annotation["missense"], "different amino acid predicted"),
        stat("Other coding consequences", annotation["other_coding"], other or "stop or start codon changes"),
        stat("Noncoding", annotation["noncoding"], "outside annotated coding regions (never counted as synonymous)"),
        stat("Unverified or unannotatable", annotation["evaluated"] - annotation["annotated"],
             f"{annotation['not_annotated']} not annotated, {annotation['partly']} partly annotated; reasons in the "
             "downloads"),
    ]
    if annotation["overlap_differs"]:
        view["stats"].insert(5, stat("Different consequences in overlapping coding regions",
                                     annotation["overlap_differs"], "every consequence is kept"))
    if annotation["shared_codons"]:
        view["stats"].append(stat("SNVs sharing a codon with another reported SNV", annotation["shared_codons"],
                                  "combined change not predicted without linkage evidence"))
    references = "; ".join(f"{item['Reference_ID']} ({item['Organism']}"
                           + (f", segment {item['Segment']}" if item["Segment"] else "") + ")"
                           for item in annotation["references"])
    view["source"] = (f"{annotation['format']} {annotation['source_file']}: {references}. "
                      f"{annotation['cds_supported']} of {annotation['cds']} coding regions passed the structure "
                      "and translation checks.")
    if annotation["agreement"]:
        view["agreement"] = ("Independent check against the codons and amino acids supplied with the input ("
                             + annotation["comparison_source"] + "): "
                             + ", ".join(f"{name.lower()} {count}" for name, count in annotation["agreement"].items())
                             + " (sprint5_annotation_comparison.csv).")
    table = annotation["table"]
    shown = table[annotation["table_columns"]].rename(columns=ANNOTATION_TABLE_HEADINGS)
    view["table"] = table_view(shown, limit=TABLE_ROWS, expand=75,
                               title="Observed mutations: predicted consequence and sequencing evidence",
                               filename="sprint5_annotated_mutations.csv")
    detail_columns = ["Variant_Key", "Gene", "CDS_ID", "Product", "Strand", "Coding_Change", "Codon_Number",
                      "Codon_Position", "Ref_Codon", "Alt_Codon", "Ref_Amino_Acid", "Alt_Amino_Acid",
                      "Amino_Acid_Change", "Mature_Peptide", "Mature_Peptide_Position", "Consequence",
                      "Consequence_Detail", "Linkage_Note"]
    view["details"] = table_view(tables.get("sprint5_mutation_annotations.csv"), columns=detail_columns,
                                 limit=25, title="Technical details: one row per mutation per coding region",
                                 filename="sprint5_mutation_annotations.csv")
    view["sources"] = table_view(tables.get("sprint5_annotation_sources.csv"), columns=[
        "Reference_ID", "Segment", "Gene", "CDS_ID", "Product", "Strand", "Parts_1_Based", "Codons", "Structure",
        "Status", "Status_Detail", "Verification"], limit=40,
        title="Annotated coding regions and their checks", filename="sprint5_annotation_sources.csv")
    view["figures"] = figures
    view["downloads"] = [name for name in ANNOTATION_DOWNLOADS if name in tables
                         and (available is None or name in available)]
    return view


def section(identifier, number, title, purpose):
    return {"id": identifier, "number": number, "title": title, "purpose": purpose,
            "sprint": SECTION_SPRINTS.get(identifier, ""),
            "stats": [], "figures": [], "tables": [], "limitations": [], "downloads": [],
            "unavailable": None}


def analysis_mode(result):
    """Plain description of which kinds of evidence this run had."""

    individual, pooled = result["individual"] is not None, result["pooled"] is not None
    if result["comparison"] is not None:
        return "Combined analysis (Sprint 2 + Sprint 3) on one verified reference"
    if individual and pooled:
        return "Individual and pooled evidence, analyzed separately (no shared verified reference)"
    if pooled:
        return "Pooled sequencing analysis (Sprint 3)"
    return "Individual sequence analysis (Sprint 2)"


def major_findings(result):
    """Short plain-language findings built from the actual summary values."""

    summary = result["summary"]
    findings = []
    if "individual_samples" in summary:
        findings.append(f"{summary['individual_samples']} individual sequence records were analyzed; "
                        f"{summary.get('individual_mapped', 0)} could be placed on the reference genome.")
    if "pooled_samples" in summary:
        findings.append(f"{summary['pooled_samples']} pooled sample records were analyzed; "
                        f"{summary.get('pooled_samples_included', 0)} were kept after quality control.")
    counts = summary.get("category_counts")
    if result["comparison"] is not None and counts:
        others = ", ".join(f"{value} {short_category(name).lower()}" for name, value in counts.items()
                           if not name.startswith("Shared"))
        findings.append(f"{summary.get('shared_variants', 0)} variant(s) were observed in both sources"
                        + (f"; the others: {others}." if others else "."))
    if summary.get("population_method"):
        findings.append(f"Population structure: {summary['population_method']}.")
    aware = result.get("depth_aware")
    if aware is not None:
        findings.append(f"Sprint 5: {aware['with_interval']} of {aware['observations_counted']} pooled "
                        "observations have a read-sampling uncertainty interval; "
                        f"{aware['eligible_series']} mutation-site series have known frequencies on two or "
                        "more dates.")
        phase2 = aware.get("uncertainty")
        if phase2:
            findings.append("Sprint 5 Phase 2: " + (
                f"overdispersion was estimated from technical replicates for {phase2['calibrated_specimens']} "
                f"specimen(s) ({phase2['calibrated_observations']} calibrated intervals)."
                if phase2["calibrated_specimens"] else
                "technical overdispersion could not be estimated (no usable technical replicates); "
                "binomial intervals are the baseline."))
    if "date_start" in summary:
        findings.append(f"Collection dates run from {summary['date_start']} to {summary['date_end']}.")
    return findings


def shared_variant_story(result):
    """Follow one shared variant from the reference to both sources (built from the results)."""

    comparison = result["comparison"]
    if comparison is None:
        return None
    table = comparison["table"]
    shared = table[table["Category"].str.startswith("Shared")]
    if shared.empty:
        return None
    row = shared.sort_values(["Pooled_Reported", "Individual_Carriers"], ascending=False).iloc[0]
    key = row["Variant_Key"]
    people = result["tables"].get("individual_variants.csv")
    pools = result["tables"].get("pooled_variants_with_qc.csv")
    people = people[people["Variant_Key"].eq(key)] if people is not None else pd.DataFrame()
    if pools is not None:
        pools = pools[pools["Variant_Key"].eq(key) & pools["Analysis_Included"].astype(bool)]
    else:
        pools = pd.DataFrame()

    def listing(frame):
        return "; ".join(f"{item['Sample_ID']} ({format_value(item['Collection_Date'])}, "
                         f"{format_value(item['Location_Name'])})" for item in frame.to_dict("records"))

    frequencies = pools["Alt_Frequency"].dropna() if not pools.empty else pd.Series(dtype=float)
    days = row["Days_Pooled_Minus_Individual"]
    km = row["Nearest_Km_Between_Carriers"]
    return {
        "key": key, "reference": row["Reference_ID"], "position": int(row["Position"]),
        "ref": row["Ref"], "alt": row["Alt"],
        "individual": listing(people) or MISSING, "individual_count": len(people),
        "pooled": listing(pools) or MISSING, "pooled_count": len(pools),
        "frequency_range": (f"{frequencies.min():.3g} to {frequencies.max():.3g} "
                            f"(median {frequencies.median():.3g})") if not frequencies.empty else MISSING,
        "category": row["Category"],
        "first_individual": format_value(row["First_Individual_Date"]),
        "first_pooled": format_value(row["First_Pooled_Date"]),
        "days": None if pd.isna(days) else int(days),
        "window_pairs": format_value(row["Date_Pairs_Within_Window"]),
        "individual_places": format_value(row["Individual_Carrier_Locations"]),
        "pooled_places": format_value(row["Pooled_Reporting_Locations"]),
        "km": None if pd.isna(km) else round(float(km), 1),
    }


def build_dashboard(result, download_id, available=None):
    """Organize results into the dashboard's eight sections.

    `available` is the set of file names actually stored for download; files
    that were too large to keep are not linked.
    """

    tables = result["tables"]
    notes = result["notes"]
    summary = result["summary"]

    def figures_for(*names):
        return [view for view in (figure_view(result, name, download_id) for name in names) if view]

    def table(name, title, **options):
        options.setdefault("limit", TABLE_ROWS)
        options.setdefault("expand", EXPANDED_ROWS)
        return table_view(tables.get(name), title=title, filename=name if name in tables else None,
                          **options)

    def finish(item, figure_names, table_items, note_keys, download_names):
        item["figures"] = figures_for(*figure_names)
        item["tables"] = [entry for entry in table_items if entry]
        for key in note_keys:
            item["limitations"] += notes.get(key, [])
        item["downloads"] = [name for name in download_names
                             if (name in tables or name[:-4] in result["figures"])
                             and (available is None or name in available)]
        return item

    pooled = result["pooled"]
    individual = result["individual"]
    comparison = result["comparison"]
    population = result.get("population")

    # 1. Overview -------------------------------------------------------------
    overview = section("overview", 1, "Overview",
                       "What was analyzed, how the two kinds of genomic evidence were connected, "
                       "and which analyses ran.")
    overview["stats"] = [
        stat("Individual sequences", summary.get("individual_samples"),
             f"{summary['individual_mapped']} mapped to the reference"
             if "individual_mapped" in summary else ""),
        stat("Pooled samples", summary.get("pooled_samples"),
             f"{summary['pooled_samples_included']} kept after QC" if "pooled_samples_included" in summary else ""),
        stat("Reference FASTA", summary.get("reference")),
        stat("Variants observed in both sources", summary.get("shared_variants"),
             "needs both sources on one verified reference" if comparison is None else
             (f"{summary['incomparable_variants']} more could not be compared: reference not verified"
              if summary.get("incomparable_variants") else "")),
        stat("Collection dates", f"{summary['date_start']} to {summary['date_end']}"
             if "date_start" in summary else None),
    ]
    overview["integration"] = result.get("integration", [])
    overview["analyses"] = result.get("analyses", [])
    overview["mode"] = analysis_mode(result)
    overview["findings"] = major_findings(result)
    finish(overview, ["completeness"], [
        table_view(pd.DataFrame(result.get("completeness", [])), title="Metadata completeness by source")
        if result.get("completeness") else None,
    ], [], ["analysis_status.csv", "completeness.png"])

    # 2. Cross-source evidence ------------------------------------------------
    cross = section("cross-source", 2, "Cross-source evidence",
                    "Which reference-based variants appear in both individual sequences and pooled "
                    "samples, which are unique to one source, and which cannot be resolved because "
                    "coverage is missing. This is the integration of the two surveillance approaches.")
    if comparison is None:
        cross["unavailable"] = ("Cross-source comparison needs individual sequences, a reference "
                                "FASTA, and a pooled SNV table on the same verified reference.")
    else:
        counts = comparison["counts"].set_index("Category")["Variants"]
        cross["stats"] = [stat(short_category(name), int(value), name.split(":", 1)[1].strip())
                          for name, value in counts.items()]
    finish(cross, ["cross_source", "evidence_matrix", "cross_source_positions", "cross_source_timeline"], [
        table_view(comparison["counts"], title="Variants per evidence category") if comparison else None,
        table("cross_source_comparison.csv", "Variant-by-variant comparison", limit=60, columns=[
            "Variant_Key", "Category", "Individual_Carriers", "Individual_Reference_Calls",
            "Individual_Not_Callable", "Pooled_Reported", "Pooled_Adequate_Depth_Not_Reported",
            "Pooled_Low_Or_No_Coverage", "Pooled_Coverage_Unknown_Or_QC_Excluded",
            "Pooled_Frequency_Median", "First_Individual_Date", "First_Pooled_Date",
            "Days_Pooled_Minus_Individual", "Date_Pairs_Within_Window",
            "Individual_Carrier_Locations", "Pooled_Reporting_Locations",
            "Nearest_Km_Between_Carriers", "Note"]),
        table("variant_signature_overlap.csv",
              "Overlap with documented variant profiles (not lineage assignment)"),
    ], ["cross_source"], ["cross_source_comparison.csv", "cross_source_evidence_by_sample.csv",
                          "variant_signature_overlap.csv", "evidence_matrix.png"])
    cross["limitations"].append("Short-read pooled data cannot show which mutations occur together "
                                "in one genome, so profile overlap does not establish that a lineage "
                                "was present. Dates are context, not proof of which source saw a "
                                "change first.")

    # 3. Genomic variation ----------------------------------------------------
    variation = section("variation", 3, "Genomic variation",
                        "Which reference positions vary, in how many samples of each source, with "
                        "which alternative nucleotide, and in which annotated regions (inspired by "
                        "Fontenele et al. Fig 3).")
    positions = tables.get("variant_positions.csv")
    if positions is not None:
        positions = positions.sort_values(["Pooled_Samples", "Individual_Carriers", "Position"],
                                          ascending=[False, False, True])
        variation["stats"] = [
            stat("Distinct variants on the reference", len(positions)),
            stat("Positions with several pooled alternative bases",
                 int(positions.loc[positions["Multiple_Alt_At_Position"], "Position"].nunique())),
            stat("Variants carried by individual sequences",
                 int(positions["Individual_Carriers"].gt(0).sum()) if individual is not None else None),
        ]
    elif individual is None:
        variation["unavailable"] = "Reference-based variation needs a reference FASTA or a pooled SNV table."
    finish(variation, ["genome_positions"], [
        table_view(positions, limit=40, title="Variants by number of samples (first 40)",
                   filename="variant_positions.csv"),
        table("individual_indels.csv", "Insertions and deletions in individual sequences (not compared)"),
    ], ["variation"], ["variant_positions.csv", "individual_variants.csv", "pooled_variants_with_qc.csv",
                       "genome_positions.png"])
    variation["limitations"].append("A published novelty flag is shown only when it comes from the "
                                    "source table. A variant absent from an uploaded comparison set "
                                    "is never labeled newly discovered.")

    # 4. Temporal patterns ----------------------------------------------------
    temporal = section("temporal", 4, "Temporal patterns",
                       "How reported variant frequencies change across collection dates and sites, "
                       "when each source was sampled, and where observations are missing.")
    trends = tables.get("temporal_trends.csv")
    if trends is not None:
        trends = trends.sort_values(["Reported", "Variant_Key"], ascending=[False, True])
        temporal["stats"] = [
            stat("Variant-site series with 2+ reports", int(trends["Reported"].ge(2).sum())),
            stat("Largest change between first and last known frequency",
                 float(trends["Change_In_Callable_Frequency"].abs().max())
                 if trends["Change_In_Callable_Frequency"].notna().any() else None,
                 "reported values and base-count non-detections only"),
        ]
    elif pooled is None and "sampling_timeline" not in result["figures"]:
        temporal["unavailable"] = "Temporal analysis needs collection dates."
    finish(temporal, ["sampling_timeline", "trajectories", "frequency_changes", "observation_status",
                      "heatmap"], [
        table_view(trends, limit=30, title="Reported-frequency trends by site (first 30)",
                   filename="temporal_trends.csv"),
        table("site_date_summary.csv", "Pooled samples by site and date"),
    ], ["temporal"], ["temporal_trends.csv", "site_date_summary.csv", "pooled_observation_matrix.csv"])
    temporal["limitations"].append("Lines connect sample-level observations at their real, irregular "
                                   "dates; they are not estimates of community prevalence, and "
                                   "slopes are descriptive.")

    # 5. Geographic context ---------------------------------------------------
    geography = section("geography", 5, "Geographic context",
                        "Where samples were collected and how precisely each location is known.")
    located = tables.get("sample_locations.csv")
    if located is not None and not located.empty:
        counts = located["Precision_Category"].value_counts()
        geography["stats"] = [stat(name, int(counts.get(name, 0))) for name in [
            "Source-provided coordinates", "Approximate (GeoNames place)", "Broad administrative area",
            "Named place or site, no coordinates", "No location"]]
    finish(geography, ["geography"], [
        table("geographic_context.csv", "Samples by documented location"),
    ], ["geography"], ["geographic_context.csv", "sample_locations.csv"])

    # 6. Sequencing quality ---------------------------------------------------
    quality = section("quality", 6, "Sequencing quality",
                      "How much of each sample can be trusted, how deep the sequencing was, which "
                      "records were excluded, and which of the paper's quality steps could actually "
                      "be checked from these files.")
    if pooled is not None:
        variants = pooled["variants"]
        quality["stats"] = [
            stat("Pooled records retained", int(variants["Analysis_Included"].sum()),
                 f"of {len(variants)}"),
            stat("Pooled samples with known genome coverage",
                 int(pooled["samples"]["Genome_Coverage_Percent"].notna().sum()),
                 f"of {len(pooled['samples'])}"),
            stat("Pooled reference", "; ".join(sorted(set(pooled["reference_status"]["Status"])))
                 if pooled.get("reference_status") is not None else None,
                 "positions and Ref bases checked against a supplied FASTA only when verified"),
            stat("Median depth at reported variants",
                 float(pd.to_numeric(variants["Read_Depth"], errors="coerce").median())
                 if pd.to_numeric(variants["Read_Depth"], errors="coerce").notna().any() else None,
                 "reads"),
        ]
    if individual is not None:
        mapped = individual["samples"]["Variant_Profile_Status"].eq("Mapped to reference")
        quality["stats"].append(stat("Individual sequences mapped to the reference", int(mapped.sum()),
                                     f"of {len(mapped)}"))
    if pooled is None and individual is None:
        quality["unavailable"] = "Quality summaries need a reference FASTA or a pooled SNV table."
    finish(quality, ["quality", "qc_outcomes", "depth_distribution", "callable_positions",
                     "variants_per_sample", "depth_vs_snvs"], [
        table("paper_method_reproduction.csv", "Which Fontenele et al. (2021) methods were reproduced"),
        table("reference_verification.csv", "Pooled reference verification"),
        table("variants_per_sample.csv", "Pooled samples: SNVs, coverage, and QC",
              columns=["Sample_ID", "Collection_Date", "Location_Name", "SNV_Records",
                       "Distinct_SNVs", "Distinct_SNVs_After_QC", "Genome_Coverage_Percent",
                       "Mean_Depth", "Quality_Status"]),
        table("individual_samples.csv", "Individual sequences: reference mapping",
              columns=["Sample_ID", "Collection_Date", "Location_Name", "Sequenced_Span",
                       "Reference_Identity_Percent", "Genome_Coverage_Percent",
                       "Quality_Status", "Variant_Profile_Status"]),
    ], ["quality", "figures"], ["paper_method_reproduction.csv", "reference_verification.csv",
                                "pooled_variants_with_qc.csv",
                                "pooled_samples_with_qc.csv", "variants_per_sample.csv",
                                "individual_samples.csv"])
    fit = result.get("regression")
    quality["regression"] = fit

    # 7. Population structure -------------------------------------------------
    structure = section("population", 7, "Population structure",
                        "Do pooled samples group by place or time? The paper used Yue–Clayton "
                        "dissimilarities of A/C/G/T frequencies summed over all reference positions, "
                        "followed by principal coordinates analysis (PCoA; Fig 5).")
    structure["method"] = None
    if population is None:
        structure["unavailable"] = ("Population structure needs a pooled SNV table; the paper's exact "
                                    "PCoA also needs per-position A/C/G/T counts.")
    elif not population["available"]:
        structure["unavailable"] = population["message"]
    else:
        pcoa = population["pcoa"]
        structure["method"] = population["method"]
        structure["exploratory"] = population["method"].startswith("Exploratory")
        explained = pcoa["explained"]
        structure["stats"] = [
            stat("Samples ordinated", len(population["coordinates"])),
            stat("Axis 1 (% of positive eigenvalues)",
                 round(explained[0] * 100, 1) if len(explained) > 0 else None),
            stat("Axis 2 (% of positive eigenvalues)",
                 round(explained[1] * 100, 1) if len(explained) > 1 else None),
            stat("Negative eigenvalues", pcoa["negative_eigenvalues"]),
        ]
    finish(structure, ["pcoa"], [
        table("pcoa_coordinates.csv", "PCoA coordinates"),
        table("nucleotide_diversity.csv",
              "Mean per-position nucleotide diversity (1 - sum of squared A/C/G/T frequencies)"),
        table("yue_clayton_pairs.csv", "Shared covered positions for each sample pair (paper method)"),
    ], ["population"], ["pcoa_coordinates.csv", "dissimilarity_matrix.csv", "nucleotide_diversity.csv",
                        "yue_clayton_pairs.csv"])
    if structure.get("exploratory"):
        for item in structure["figures"]:
            if item["name"] == "pcoa" and item["guide"]:
                item["guide"] = {**item["guide"],
                                 "origin": "Exploratory alternative (not the paper's Yue–Clayton method)"}

    # 8. Individual sequence relationships (Sprint 2) -------------------------
    relationships = section("relationships", 8, "Individual sequence relationships",
                            "The original Sprint 2 screening: MAFFT-aligned genetic distances among "
                            "individual sequences combined with collection dates and locations.")
    finish(relationships, ["genetic_distances", "relationship_network"], [
        table("pairwise_distances.csv", "Pairwise distances (first 50 pairs)", limit=50, columns=[
            "Sample_ID_1", "Sample_ID_2", "Pathogen", "Gene_Segment", "P_Distance",
            "Differences", "Comparable_Positions", "Percent_Identity", "Days_Apart", "Km_Apart", "Same_Location",
            "Geographic_Precision_1", "Geographic_Precision_2"]),
    ], ["relationships"], ["pairwise_distances.csv"])

    # 9. Sprint 5: depth-aware variant analysis ------------------------------
    depth = section("depth-aware", 9, "Sprint 5 — Depth-Aware Variant Analysis",
                    "This analysis examines how strongly sequencing data support observed pathogen mutations "
                    "and how their frequencies change over time.")
    aware = result.get("depth_aware")
    depth["messages"] = []
    depth["annotation"] = annotation_view(result.get("annotation"), figures_for("mutation_annotation_overview"),
                                          tables, available, notes.get("annotation", []))
    if aware is None:
        depth["unavailable"] = ("This optional analysis needs a processed pooled SNV table: sample, position, "
                                "reference and alternative base, and frequency. Read depth (and ideally read "
                                "counts) is needed for uncertainty, and collection dates for time series.")
    else:
        depth["messages"] = aware["messages"]
        counts_detail = (f"{aware['supplied_counts']} supplied or from base counts, "
                         f"{aware['reconstructed_counts']} reconstructed from frequency x depth")
        depth["stats"] = [
            stat("Samples evaluated", aware["samples"]),
            stat("SNV observations evaluated", aware["observations_counted"],
                 f"{aware['observations'] - aware['observations_counted']} repeated record(s) counted once"
                 if aware["observations"] != aware["observations_counted"] else ""),
            stat("Observations with usable depth", aware["with_depth"]),
            stat("Observations with verified or reconstructed alternative read counts",
                 aware["supplied_counts"] + aware["reconstructed_counts"], counts_detail),
            stat("Observations with calculated uncertainty", aware["with_interval"], "95% Wilson interval"),
            stat("Observations missing required information", aware["missing_required"],
                 "no interval; see the reasons table"),
            stat("Variants eligible for time-series analysis", aware["eligible_variants"],
                 f"known values on 2+ dates at one site; {aware['eligible_series']} mutation-site series"),
        ]
        depth["grouping"] = grouping_view(aware.get("grouping"), figures_for(*GROUPING_FIGURES), tables, available)
        phase2 = aware.get("uncertainty")
        if phase2:
            depth["uncertainty"] = uncertainty_view(phase2)
            depth["stats"].append(stat(
                "Overdispersion (rho) estimated", f"{phase2['calibrated_specimens']} specimen(s)"
                if phase2["calibrated_specimens"] else "Not estimated",
                f"{phase2['calibrated_observations']} calibrated beta-binomial intervals"
                if phase2["calibrated_specimens"] else "needs documented technical replicates (Specimen_ID + distinct Library_ID)"))
    finish(depth, ["vaf_time_series", "uncertainty_models", "replicate_concordance", "vaf_depth",
                   "depth_dependence"], [
        table("sprint5_overdispersion_calibration.csv", "Overdispersion estimated from technical replicates",
              columns=["Specimen_ID", "Replicate_1", "Replicate_2", "Status", "Sites_Used", "Rho_Estimate",
                       "Rho_Lower_95", "Rho_Upper_95", "Rho_Used_For_Intervals", "Plugin_Likelihood_Rho",
                       "Status_Detail"]),
        table("sprint5_missing_information.csv", "What the data lack and what that prevents"),
        table("sprint5_temporal_summary.csv", "Mutation frequencies over time by site (first 30)", limit=30,
              columns=["Variant_Key", "Location_Name", "Samples", "Detections", "Supported_Non_Detections",
                       "Unknown_Or_Uncallable", "Known_Dated_Time_Points", "Earliest_Date", "Earliest_VAF",
                       "Latest_Date", "Latest_VAF", "Absolute_Change", "Slope_Per_Week", "Series_Status",
                       "Interpretation_Note"]),
        table("sprint5_vaf_observations.csv", "VAF, read support, and uncertainty per observation (first 25)",
              columns=["Sample_ID", "Collection_Date", "Variant_Key", "Reported_Frequency", "Alt_Read_Count",
                       "Alt_Count_Source", "Total_Depth", "VAF", "Wilson_Lower", "Wilson_Upper",
                       "Interval_Status", "QC_Status", "Validation_Warnings"]),
        table("sprint5_depth_dependence.csv", "Same VAF at different depths (illustration)"),
    ], ["depth_aware"], ["sprint5_vaf_observations.csv", "sprint5_temporal_summary.csv",
                         "sprint5_missing_information.csv", "sprint5_depth_dependence.csv",
                         "sprint5_uncertainty_models.csv", "sprint5_overdispersion_sensitivity.csv",
                         "sprint5_overdispersion_calibration.csv", "sprint5_calibration_sites.csv"])
    if aware is not None:
        depth["limitations"] += DEPTH_AWARE_LIMITATIONS + OVERDISPERSION_LIMITATIONS

    cross["story"] = shared_variant_story(result)
    view = {"summary": summary, "notes": notes, "download_id": download_id,
            "sections": [overview, cross, variation, temporal, geography, quality, structure,
                         relationships, depth]}
    for item in view["sections"]:
        view[item["id"].replace("-", "_")] = item
    view["provenance"] = table_view(tables.get("provenance.csv"), title="Provenance record")
    view["downloads"] = sorted(tables) + sorted(f"{name}.png" for name in result["figures"])
    # Phase 4A: consolidated downloads, the pre-run readiness, and the import record.
    view["report_downloads"] = ([XLSX_NAME, PDF_NAME] if available is None or SUMMARY_NAME in available else [])
    view["export_notes"] = result.get("export_notes", [])
    view["findings"] = major_findings(result)
    view["readiness"] = (result.get("import_details") or {}).get("readiness") or []
    view["import_provenance"] = (table_view(tables["import_provenance.csv"], limit=TABLE_ROWS, title="Import provenance",
                                            filename="import_provenance.csv", expand=EXPANDED_ROWS)
                                 if len(tables.get("import_provenance.csv", [])) else None)
    view["matching"] = (table_view(tables["metadata_matching.csv"], title="Sample metadata matching",
                                   filename="metadata_matching.csv")
                        if "metadata_matching.csv" in tables else None)
    view["canonical"] = result.get("canonical")
    return view


def download_files(result, extra_files=None, size_limit=8 * 1024 * 1024, reports=True):
    """Encode tables as CSV and figures as PNG for the temporary download cache.

    Phase 4A also stores a small analysis_summary.json (findings, notes,
    readiness, figure titles). The XLSX workbook and the PDF report are built
    from these stored files only when someone downloads them
    (genomics.report_exports), so the analysis request is not slowed down.
    """

    files = {}
    skipped = []
    for name, table in result["tables"].items():
        content = safe_csv(table)  # formula-like text escaped (audit P1-2)
        if len(content) > size_limit:
            skipped.append(name)
            continue
        files[name] = (content, "text/csv")
    for name, figure in result["figures"].items():
        files[f"{name}.png"] = (figure["png"], "image/png")
    for name, content in (extra_files or {}).items():
        files[name] = (escape_csv_bytes(content), "text/csv")
    result["export_notes"] = []
    if reports and (result.get("pooled") is not None or result.get("individual") is not None):
        files[SUMMARY_NAME] = (summary_json(result, major_findings(result)), "application/json")
        large = [name for name in skipped if name.endswith(".csv")]
        if large:
            result["export_notes"].append("Tables too large for the download cache are also missing from the "
                                          "workbook: " + ", ".join(large) + ".")
    return files, skipped
