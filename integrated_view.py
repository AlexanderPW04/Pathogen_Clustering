"""Turn integrated analysis results into simple values for the dashboard page.

Nothing here calculates science; it only arranges the tables, figures, and
notes that `genomics.integrated_analysis` already produced into the eight
dashboard sections. Every section has the same parts: purpose, summary
statistics, figures with captions, tables, limitations, and downloads.
Integrated (both-source) sections come before source-specific ones.
"""

from base64 import b64encode

import pandas as pd


MISSING = "—"  # an em dash marks a value that is not available
NOT_AVAILABLE = "Not available"  # a statistic the inputs cannot provide (not zero)
WIDE_FIGURE_PIXELS = 1400


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


def table_view(table, columns=None, limit=None, title="", filename=None):
    """Return {'title', 'columns', 'rows', 'total', 'shown', 'filename'} for the template."""

    if table is None:
        return None
    data = table if columns is None else table[[name for name in columns if name in table.columns]]
    total = len(data)
    if limit is not None:
        data = data.head(limit)
    rows = [[format_value(value) for value in row] for row in data.itertuples(index=False)]
    return {"title": title, "columns": list(data.columns), "rows": rows,
            "total": total, "shown": len(rows), "filename": filename}


def png_width(png):
    """Pixel width from a PNG header (bytes 16-19), or 0 if it is not a PNG."""

    return int.from_bytes(png[16:20], "big") if png[:8] == b"\x89PNG\r\n\x1a\n" else 0


def figure_view(result, name, download_id):
    figure = result["figures"].get(name)
    if figure is None:
        return None
    return {
        "title": figure["title"],
        "caption": figure["caption"],
        "data_uri": "data:image/png;base64," + b64encode(figure["png"]).decode("ascii"),
        "filename": f"{name}.png",
        "download_id": download_id,
        # Wider than the page: shown at full size in a scrolling box so labels stay readable.
        "wide": png_width(figure["png"]) > WIDE_FIGURE_PIXELS,
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


def section(identifier, number, title, purpose):
    return {"id": identifier, "number": number, "title": title, "purpose": purpose,
            "stats": [], "figures": [], "tables": [], "limitations": [], "downloads": [],
            "unavailable": None}


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
            stat("Largest change between first and last report",
                 float(trends["Change_In_Reported_Frequency"].abs().max())
                 if trends["Change_In_Reported_Frequency"].notna().any() else None),
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
        table("variants_per_sample.csv", "Pooled samples: SNVs, coverage, and QC",
              columns=["Sample_ID", "Collection_Date", "Location_Name", "SNV_Records",
                       "Distinct_SNVs", "Distinct_SNVs_After_QC", "Genome_Coverage_Percent",
                       "Mean_Depth", "Quality_Status"]),
        table("individual_samples.csv", "Individual sequences: reference mapping",
              columns=["Sample_ID", "Collection_Date", "Location_Name", "Sequenced_Span",
                       "Reference_Identity_Percent", "Genome_Coverage_Percent",
                       "Quality_Status", "Variant_Profile_Status"]),
    ], ["quality", "figures"], ["paper_method_reproduction.csv", "pooled_variants_with_qc.csv",
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
    ], ["population"], ["pcoa_coordinates.csv", "dissimilarity_matrix.csv", "nucleotide_diversity.csv"])

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

    view = {"summary": summary, "notes": notes, "download_id": download_id,
            "sections": [overview, cross, variation, temporal, geography, quality, structure,
                         relationships]}
    for item in view["sections"]:
        view[item["id"].replace("-", "_")] = item
    view["provenance"] = table_view(tables.get("provenance.csv"), title="Provenance record")
    view["downloads"] = sorted(tables) + sorted(f"{name}.png" for name in result["figures"])
    return view


def download_files(result, extra_files=None, size_limit=8 * 1024 * 1024):
    """Encode tables as CSV and figures as PNG for the temporary download cache."""

    files = {}
    skipped = []
    for name, table in result["tables"].items():
        content = table.to_csv(index=False).encode("utf-8")
        if len(content) > size_limit:
            skipped.append(name)
            continue
        files[name] = (content, "text/csv")
    for name, figure in result["figures"].items():
        files[f"{name}.png"] = (figure["png"], "image/png")
    for name, content in (extra_files or {}).items():
        files[name] = (content, "text/csv")
    return files, skipped
