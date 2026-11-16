"""Arrange a finished dashboard into a few result categories that fit what actually ran.

`integrated_view.build_dashboard` still builds the detailed sections (stats,
figures, tables, limitations) exactly as before. This module only decides
where each piece is shown, using `genomics.result_availability`:

- which categories appear for this run (an empty category is never rendered);
- which figures are shown first, which are one click away under
  "Additional figures", and which are not shown because the data cannot
  support them (their PNGs stay downloadable);
- one grouped list of downloads instead of a list under every section.

No values are calculated here.
"""

CATEGORIES = {
    "overview": ("Overview", ""),
    "cross-source": ("Cross-Source Variant Evidence",
                     "Variants compared between individual sequences and pooled samples on one verified reference "
                     "genome: seen in both, in one source only, or unresolved because coverage is missing."),
    "relationships": ("Individual Sequence Relationships",
                      "MAFFT-aligned genetic distances among individual sequences, combined with collection dates "
                      "and locations to screen for related pairs and potential clusters."),
    "variation": ("Genomic Variation",
                  "Observed SNVs, where they fall on the reference genome, and their predicted gene and protein "
                  "consequences."),
    "temporal": ("Temporal and Geographic Patterns",
                 "When and where samples were collected, and how known mutation frequencies change over time."),
    "quality": ("Sequencing Quality and Uncertainty",
                "Read depth, quality control, variant allele frequency (VAF), and how precisely each frequency is "
                "known."),
    "advanced": ("Advanced Analysis",
                 "Population structure, experimental mutation-frequency profiles, and detailed statistical methods. "
                 "These results are more method-dependent; read their labels and limitations."),
    "downloads": ("Downloads and Provenance", ""),
}
ORDER = {
    "pooled": ["overview", "variation", "temporal", "quality", "advanced", "downloads"],
    "individual": ["overview", "relationships", "variation", "temporal", "quality", "advanced", "downloads"],
    "combined": ["overview", "cross-source", "relationships", "variation", "temporal", "quality", "advanced",
                 "downloads"],
    "separate": ["overview", "relationships", "variation", "temporal", "quality", "advanced", "downloads"],
}
# Figures drawn inside their own panels (annotation, experimental grouping), not in the figure lists.
PANEL_FIGURES = {"mutation_annotation_overview", "mutation_profile_heatmap", "mutation_group_heatmap",
                 "mutation_group_trajectories", "grouping_comparison"}
# Which category each section's download names belong to.
DOWNLOAD_CATEGORY = {"overview": "downloads", "cross-source": "cross-source", "variation": "variation",
                     "temporal": "temporal", "geography": "temporal", "quality": "quality",
                     "population": "advanced", "relationships": "relationships", "depth-aware": "quality"}
DOWNLOAD_OVERRIDES = {"sprint5_temporal_summary.csv": "temporal", "sprint5_uncertainty_models.csv": "advanced",
                      "sprint5_overdispersion_sensitivity.csv": "advanced",
                      "sprint5_overdispersion_calibration.csv": "advanced", "sprint5_calibration_sites.csv": "advanced"}
GROUP_TITLES = {"downloads": "Analysis record, input mapping, and provenance"}


NAV_LABELS = {"overview": "Overview", "cross-source": "Cross-source evidence", "relationships": "Sequence relationships",
              "variation": "Genomic variation", "temporal": "Time and place", "quality": "Quality and uncertainty",
              "advanced": "Advanced", "downloads": "Downloads"}


def category(identifier, **content):
    title, purpose = CATEGORIES[identifier]
    item = {"id": identifier, "title": title, "nav": NAV_LABELS[identifier], "purpose": purpose, "stats": [],
            "statements": [], "primary": [],
            "secondary": [], "tables": [], "limitations": [], "notes": []}
    item.update(content)
    return item


def figures_in(view, roles, name_list, role):
    return [view["figure_views"][name] for name in name_list
            if name in view["figure_views"] and roles.get(name, [None, None])[1] == role]


def placed(roles, identifier, role):
    return [name for name, (where, kind) in roles.items() if where == identifier and kind == role
            and name not in PANEL_FIGURES]


def without(stats, *labels):
    return [item for item in stats if item["label"] not in labels]


def build_categories(view, availability):
    """The ordered categories for this run; each holds only content that exists."""

    roles = availability["figure_roles"]
    available = availability["available"]
    sections = {item["id"]: item for item in view["sections"]}
    order = ORDER[availability["mode"]]

    def figures(identifier):
        return (figures_in(view, roles, placed(roles, identifier, "primary"), "primary"),
                figures_in(view, roles, placed(roles, identifier, "secondary"), "secondary"))

    built = {}
    overview = sections["overview"]
    built["overview"] = category("overview", stats=[], secondary=figures("overview")[1],
                                 tables=overview["tables"])

    if available["cross_source"]:
        cross = sections["cross-source"]
        primary, secondary = figures("cross-source")
        built["cross-source"] = category("cross-source", stats=cross["stats"], primary=primary, secondary=secondary,
                                         tables=cross["tables"], limitations=cross["limitations"],
                                         story=cross.get("story"))

    if available["relationships"] or available["individual"]:
        relationships = sections["relationships"]
        primary, secondary = figures("relationships")
        built["relationships"] = category("relationships", primary=primary, secondary=secondary,
                                          tables=relationships["tables"], limitations=relationships["limitations"],
                                          sprint2=True)

    variation = sections["variation"]
    annotation = sections["depth-aware"].get("annotation")
    if available["variation"] or available["annotation"]:
        primary, secondary = figures("variation")
        built["variation"] = category("variation", stats=variation["stats"], primary=primary, secondary=secondary,
                                      tables=variation["tables"], limitations=variation["limitations"],
                                      annotation=annotation)
    elif annotation and availability["mode"] != "individual":
        # A pooled run without a usable annotation still says so once, under Genomic Variation.
        built["variation"] = category("variation", annotation=annotation)

    temporal, geography = sections["temporal"], sections["geography"]
    depth = sections["depth-aware"]
    temporal_tables = temporal["tables"] + [table for table in depth["tables"]
                                            if table["filename"] == "sprint5_temporal_summary.csv"]
    primary, secondary = figures("temporal")
    statements = [availability["suppressed_figures"]["geography"]] if "geography" in availability[
        "suppressed_figures"] else []
    if primary or secondary or temporal["stats"] or available["locations"]:
        # Location precision: only the categories that occur (the full table stays under "Detailed tables").
        located = [item for item in geography["stats"] if item["value"] != "0"]
        built["temporal"] = category("temporal", stats=temporal["stats"] + located, primary=primary,
                                     secondary=secondary, statements=statements,
                                     tables=temporal_tables + geography["tables"],
                                     limitations=temporal["limitations"] + geography["limitations"])

    quality = sections["quality"]
    primary, secondary = figures("quality")
    coverage_statement = availability["suppressed_figures"].get("quality") or (
        "Genome-wide coverage was not supplied, so coverage percentages could not be calculated."
        if available["pooled"] and not available["coverage"] else None)
    quality_stats = without(quality["stats"], "Pooled samples with known genome coverage") if (
        available["pooled"] and not available["coverage"]) else quality["stats"]
    uncertainty_stats = [item for item in depth["stats"] if item["label"] in {
        "Observations with calculated uncertainty", "Overdispersion (rho) estimated"}]
    quality_tables = quality["tables"] + [table for table in depth["tables"] if table["filename"] in {
        "sprint5_vaf_observations.csv", "sprint5_missing_information.csv", "sprint5_depth_dependence.csv"}]
    if primary or secondary or quality["stats"] or depth["stats"]:
        built["quality"] = category("quality", stats=quality_stats + uncertainty_stats, primary=primary,
                                    secondary=secondary,
                                    statements=[coverage_statement] if coverage_statement else [],
                                    tables=quality_tables, regression=quality.get("regression"),
                                    uncertainty=depth.get("uncertainty"), depth_stats=depth["stats"],
                                    depth_messages=depth.get("messages", []),
                                    limitations=quality["limitations"] + depth["limitations"])

    population = sections["population"]
    primary, secondary = figures("advanced")
    advanced_tables = (population["tables"] if available["population"] else []) + [
        table for table in depth["tables"] if table["filename"] == "sprint5_overdispersion_calibration.csv"]
    grouping = depth.get("grouping")
    if primary or secondary or grouping or available["population"]:
        built["advanced"] = category(
            "advanced", stats=population["stats"] if available["population"] else [], primary=primary,
            secondary=secondary, tables=advanced_tables, grouping=grouping,
            method=population.get("method"), exploratory=population.get("exploratory"),
            limitations=population["limitations"] if available["population"] else [])

    built["downloads"] = category("downloads")
    return [built[identifier] for identifier in order if identifier in built]


def download_groups(view, availability, stored):
    """[(title, [file names])] in category order; every stored file appears exactly once."""

    roles = availability["figure_roles"]
    where = {}
    for item in view["sections"]:
        for name in item["downloads"]:
            where.setdefault(name, DOWNLOAD_OVERRIDES.get(name, DOWNLOAD_CATEGORY.get(item["id"], "downloads")))
    depth = next(item for item in view["sections"] if item["id"] == "depth-aware")
    for panel, target in [(depth.get("annotation"), "variation"), (depth.get("grouping"), "advanced")]:
        for name in (panel or {}).get("downloads", []):
            where.setdefault(name, target)
    for name, (target, _) in roles.items():
        where.setdefault(f"{name}.png", target)
    for name in availability["suppressed_figures"]:
        where.setdefault(f"{name}.png", "quality" if name == "quality" else "temporal")
    order = ORDER[availability["mode"]][1:] + ["downloads"]
    groups = {identifier: [] for identifier in order}
    for name in sorted(stored):
        if name.endswith((".json", ".xlsx", ".pdf")):
            continue
        target = where.get(name, "downloads")
        groups.setdefault(target, []).append(name)
    return [(GROUP_TITLES.get(identifier, CATEGORIES[identifier][0]), names)
            for identifier, names in groups.items() if names]
