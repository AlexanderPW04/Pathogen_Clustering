"""Plain-language reading guides for the dashboard (no calculations).

FIGURE_GUIDE says, for every figure, what it measures, what the axes and
colors mean, where the method comes from, and its main limitation, so the
same long explanation does not have to be repeated under each figure.
`interpretation` writes one short sentence from the actual results.
GLOSSARY explains the terms a newcomer meets on the page.
"""

import pandas as pd

PAPER = "Fontenele et al. (2021)"

FIGURE_GUIDE = {
    "completeness": {
        "measures": "How many samples in each source have collection dates, locations, coordinates, and coverage.",
        "axes": "Bars show the percent of samples that have each kind of information.",
        "colors": "Green = individual sequences; blue = pooled samples.",
        "origin": "Project addition",
        "limitation": "Shows whether information is present, not whether it is correct.",
    },
    "cross_source": {
        "measures": "How many reference-based variants fall into each evidence category.",
        "axes": "Each bar is one category; its length is the number of distinct variants.",
        "colors": "One color per category (shared, one source only, unresolved, incomparable).",
        "origin": "Project addition (connects Sprint 2 and Sprint 3)",
        "limitation": "'Only' describes the uploaded evidence; it does not mean a variant is absent from the population.",
    },
    "evidence_matrix": {
        "measures": "For each variant (row) and each sample (column), what was actually observed.",
        "axes": "Rows: variants written position:Ref>Alt. Columns: individual sequences left of the black line, pooled samples right of it.",
        "colors": "Green = variant base in an individual sequence; blue shades = pooled frequency (darker is higher); white = reference base, or adequate depth without the variant; hatching = not callable, low coverage, unknown coverage, or excluded by QC.",
        "origin": "Project addition (connects Sprint 2 and Sprint 3)",
        "limitation": "Unknown coverage is hatched, never drawn as 0%.",
    },
    "cross_source_positions": {
        "measures": "Where the comparable variants lie on the reference genome, by evidence category.",
        "axes": "Horizontal axis: reference position (1-based). One lane per category.",
        "colors": "Each point is one variant; lanes separate the categories.",
        "origin": "Project addition",
        "limitation": "Positions only; it does not show frequencies or sample counts.",
    },
    "cross_source_timeline": {
        "measures": "When each shared or source-specific variant was seen in each source.",
        "axes": "Horizontal axis: collection date. One row per variant.",
        "colors": "Green triangles = individual sequences; blue circles = pooled samples (larger = higher frequency).",
        "origin": f"Exploratory (related to {PAPER} Fig 4)",
        "limitation": "Dates are context; they do not show which source detected a change first or why.",
    },
    "genome_positions": {
        "measures": "How many samples carry each variant, along the reference genome.",
        "axes": "Horizontal axis: reference position, with annotated regions underneath. Vertical axis: number of samples with the variant.",
        "colors": "Color = alternative base; filled circle = pooled samples reporting it; open diamond = individual sequences carrying it; ring = several pooled alternative bases at one position.",
        "origin": f"Inspired by {PAPER} Fig 3",
        "limitation": "Counts samples, not frequencies; the paper drew only novel variants.",
    },
    "sampling_timeline": {
        "measures": "When and where samples were collected in each source.",
        "axes": "Horizontal axis: collection date. One row per site.",
        "colors": "Triangles = individual sequences; circles = pooled samples; open marker = excluded by QC.",
        "origin": f"Partial analogue of {PAPER} Fig 1A (no map)",
        "limitation": "Samples without a date cannot be drawn and are counted in the title.",
    },
    "trajectories": {
        "measures": "Reported variant frequencies across collection dates at each site.",
        "axes": "Horizontal axis: collection date. Vertical axis: reported alternative-allele frequency (0 to 1).",
        "colors": "One color per variant; lines join reports from the same site.",
        "origin": f"Inspired by the {PAPER} temporal analysis",
        "limitation": "Gaps are missing reports, not zeros; lines are not estimates of community prevalence.",
    },
    "frequency_changes": {
        "measures": "Change between the first and last known frequency of each variant at each site (reports, plus base-count non-detections when supplied).",
        "axes": "Horizontal axis: frequency; open marker = first known value, filled marker = last known value.",
        "colors": "Color shows whether the frequency was higher or lower at the last report.",
        "origin": "Project addition",
        "limitation": "Uses only two reports per series; descriptive, not a trend test.",
    },
    "observation_status": {
        "measures": "For each pooled sample, the share of analyzed variants that were reported, absent with adequate depth, poorly covered, or of unknown coverage.",
        "axes": "Vertical axis: percent of analyzed variants. One bar per sample.",
        "colors": "Solid = reported or adequately covered; hatched = low or unknown coverage.",
        "origin": "Project addition",
        "limitation": "'Coverage unknown' is not 0%; it means the files cannot say.",
    },
    "heatmap": {
        "measures": "Reported frequency of each variant in each pooled sample.",
        "axes": "Rows: pooled samples (site, date). Columns: variants written position:Ref>Alt.",
        "colors": "Blue shade = reported frequency; plain = adequate depth but not reported; hatched = low coverage, unknown coverage, or excluded by QC.",
        "origin": "Project addition (sample-by-variant view of the paper's SNV table)",
        "limitation": "Unknown coverage is hatched, never drawn as a zero frequency.",
    },
    "geography": {
        "measures": "How precisely each sample's location is known, and where samples with coordinates lie.",
        "axes": "Bars: number of samples per precision class. Map panel (when available): longitude and latitude.",
        "colors": "One color or pattern per precision class.",
        "origin": f"Partial analogue of {PAPER} Fig 1A",
        "limitation": "State-level labels and their centroids are too coarse for kilometer distances.",
    },
    "quality": {
        "measures": "Genome coverage of each pooled sample, and callable reference positions of each individual sequence.",
        "axes": "Vertical axes: percent of genome positions at ≥10 reads (pooled) or callable (individual).",
        "colors": "Blue = pooled samples; green = individual sequences; hatched = excluded or not mapped.",
        "origin": f"Faithful to {PAPER} Fig 1C when coverage is supplied",
        "limitation": "Without coverage metadata the pooled panel shows 'unknown', not 0%.",
    },
    "qc_outcomes": {
        "measures": "Which quality checks each pooled variant record passed, failed, or could not be evaluated.",
        "axes": "One bar per check; length = percent of pooled variant records.",
        "colors": "Green = pass; red = fail; gray = not evaluable from these files.",
        "origin": f"Checks from {PAPER} methods",
        "limitation": "Base-quality filtering cannot be verified from processed records.",
    },
    "depth_distribution": {
        "measures": "How many reads covered each reported variant position.",
        "axes": "Horizontal axis: read depth (log scale). Vertical axis: number of variant records.",
        "colors": "The dashed line marks the depth threshold used for QC.",
        "origin": f"Threshold from {PAPER} methods",
        "limitation": "Depth at unreported positions is known only when a per-position table is supplied.",
    },
    "callable_positions": {
        "measures": "For each individual sequence, the share of reference positions with a usable base.",
        "axes": "Vertical axis: percent of reference positions. One bar per sequence.",
        "colors": "Solid = callable; hatched = gaps, ambiguous bases, or outside the sequenced span.",
        "origin": "Project addition",
        "limitation": "Uncallable positions are unknown, not reference matches.",
    },
    "variants_per_sample": {
        "measures": "Number of distinct SNVs in each pooled sample after QC.",
        "axes": "Vertical axis: distinct SNVs. One bar per sample.",
        "colors": "One color per sampling site (extra sites share gray).",
        "origin": f"Faithful to {PAPER} Fig 2A",
        "limitation": "SNV counts depend on sequencing depth and coverage.",
    },
    "depth_vs_snvs": {
        "measures": "Whether samples sequenced more deeply report more SNVs.",
        "axes": "Horizontal axis: mean read depth. Vertical axis: distinct SNVs after QC.",
        "colors": "Points = samples; line = least-squares fit with a bootstrap 95% band.",
        "origin": f"Partial reproduction of {PAPER} Fig 2B (model details not stated in the paper)",
        "limitation": "Descriptive only with fewer than 10 samples.",
    },
    "vaf_time_series": {
        "measures": "How the frequency of a mutation changes over time at one site, with how precisely each value is measured.",
        "axes": "Horizontal axis: collection date. Vertical axis: variant allele frequency (VAF). One panel per mutation and site.",
        "colors": "Filled = detected; open = supported non-detection or uncalled reads from base counts; bars = 95% Wilson intervals (dashed when counts were reconstructed); gray crosses = samples with no known value.",
        "origin": "Sprint 5 (project addition)",
        "limitation": "Intervals cover read sampling only. Unknown samples are left out, not set to 0; lines are descriptive, not tests of change.",
    },
    "vaf_depth": {
        "measures": "How the uncertainty of a measured frequency depends on sequencing depth.",
        "axes": "Left: depth (log) against interval width, with theoretical curves. Right: depth (log) against observed VAF.",
        "colors": "Left: gray = observations, colored lines = theory for three VAFs. Right: darker = narrower interval.",
        "origin": "Sprint 5 (project addition; binomial Wilson interval)",
        "limitation": "Depth reflects precision, not correctness; sequencing and amplification errors are not modeled.",
    },
    "depth_dependence": {
        "measures": "The same frequency (0.20) measured with 20 to 1,000 reads.",
        "axes": "Each row: a count such as 4 of 20 reads; the bar is its 95% interval.",
        "colors": "Dot = VAF; bar = interval.",
        "origin": "Sprint 5 illustration (exact counts, not uploaded data)",
        "limitation": "Illustrates read-sampling uncertainty only.",
    },
    "uncertainty_models": {
        "measures": "How much wider frequency uncertainty becomes when sequencing libraries vary more than read sampling alone (beta-binomial overdispersion rho).",
        "axes": "Left: depth (log) against the 95% interval width at VAF 0.20 for each model. Right: each observation's binomial width against its overdispersion-aware width.",
        "colors": "Black = binomial (Wilson). Solid colors = rho calibrated from technical replicates; dashed = hypothetical rho values. Right: blue = calibrated, gray = hypothetical scenario.",
        "origin": "Sprint 5 Phase 2 (project addition; beta-binomial model as in Tonkin-Hill et al. 2021)",
        "limitation": "Hypothetical curves are scenarios, not estimates. Even calibrated rho covers library-to-library variation only, not primer, mapping, or contamination errors.",
    },
    "replicate_concordance": {
        "measures": "Agreement between two separately prepared libraries of the same specimen, used to estimate overdispersion.",
        "axes": "Replicate 1 VAF against replicate 2 VAF for each informative allele.",
        "colors": "One color and marker per specimen; the legend gives its rho estimate and bootstrap range.",
        "origin": "Sprint 5 Phase 2 (site criteria from Tonkin-Hill et al. 2021)",
        "limitation": "One rho per library assumes every site varies alike; sites in one amplicon may be correlated.",
    },
    "mutation_profile_heatmap": {
        "measures": "Each mutation's known frequency in every observation unit, without grouping.",
        "axes": "Rows: mutations in genome order (those with the most known values). Columns: observation units by date.",
        "colors": "Darker blue = higher frequency. Hatched gray = unknown (no coverage information), never 0.",
        "origin": "Sprint 5 Phase 3 (mutation-frequency profiles)",
        "limitation": "Similar-looking rows are not evidence of a shared genome or lineage.",
    },
    "mutation_annotation_overview": {
        "measures": "Where each distinct observed SNV lies on the reference and its predicted coding consequence.",
        "axes": "Top: reference position (1-based), one row per predicted category, over the annotated coding "
                "regions. Bottom: number of SNVs in each coding region, split by predicted consequence.",
        "colors": "One fixed color per category (blue synonymous, orange missense, others for stop and start "
                  "changes, gray noncoding); hatched boxes are coding regions that could not be interpreted.",
        "origin": "Sprint 5 Phase 4 (reference annotation; Sequence Ontology consequence terms)",
        "limitation": "A predicted amino-acid change is not a demonstrated biological effect, and pooled data do not "
                      "show which SNVs occur on the same genome.",
    },
    "mutation_group_heatmap": {
        "measures": "Frequencies of mutations placed in experimental groups, followed by ungrouped ones.",
        "axes": "Rows: mutations, labeled with their group (dash = not grouped). Columns: observation units by date.",
        "colors": "Darker blue = higher frequency; hatched gray = unknown, never 0; black lines separate groups.",
        "origin": "Sprint 5 Phase 3 (experimental; MAGOS-inspired grouping)",
        "limitation": "The grouping did not meet its validation criteria; a group is similar frequencies, not a lineage.",
    },
    "mutation_group_trajectories": {
        "measures": "How the mean frequency of each experimental group's mutations changes over time.",
        "axes": "Collection date against the mean frequency of the group's known members; bars span the lowest to highest member.",
        "colors": "One color and marker per group; line style = sampling site.",
        "origin": "Sprint 5 Phase 3 (experimental)",
        "limitation": "Descriptive only; unknown members are left out, so a point may rest on few mutations.",
    },
    "grouping_comparison": {
        "measures": "How the depth-aware groups and the ordinary (depth-blind) groups overlap.",
        "axes": "Rows: depth-aware groups; columns: ordinary groups; each cell counts shared mutations.",
        "colors": "Darker blue = more shared mutations; '—' = not grouped by that method.",
        "origin": "Sprint 5 Phase 3 (experimental comparison)",
        "limitation": "Agreement between two methods is not accuracy; ordinary clustering grouped unrelated mutations often in simulation.",
    },
    "pcoa": {
        "measures": "Whether pooled samples group by place or time, from differences in their nucleotide frequencies.",
        "axes": "Axes 1 and 2 are the first two principal coordinates; the percentages show how much of the positive eigenvalue total each explains. The small panel lists the eigenvalues.",
        "colors": "Color = collection date; marker shape = site (or state when there are many sites).",
        "origin": f"{PAPER} Fig 5 (paper method only with per-position A/C/G/T counts)",
        "limitation": "Closeness in the plot means similar frequency profiles, not transmission.",
    },
    "genetic_distances": {
        "measures": "Genetic p-distances between pairs of individual sequences (Sprint 2).",
        "axes": "Left: histogram of p-distances (pairs per bin). Right: p-distance against days between collection dates.",
        "colors": "Dashed lines mark the screening thresholds.",
        "origin": "Sprint 2 method",
        "limitation": "Thresholds are prototype settings, not outbreak definitions.",
    },
    "relationship_network": {
        "measures": "Which individual sequences are linked by the Sprint 2 screening.",
        "axes": "Samples sit on a circle; position has no meaning.",
        "colors": "Solid green line = meets all three cutoffs; dashed gray line = genetic and date only, geography unresolved.",
        "origin": "Sprint 2 method",
        "limitation": "A line is a screening signal, never confirmed transmission.",
    },
}

GLOSSARY = [
    ("FASTA", "A text file of DNA sequences: a '>' line with the sequence name, then the letters A, C, G, T."),
    ("Processed SNV table", "A table of single-nucleotide variants already called from sequencing reads: sample, "
                            "reference position, reference and alternative base, frequency, and read depth."),
    ("Reference genome", "One agreed sequence that gives every position a number, so variants from different "
                         "files can be compared position by position."),
    ("Genetic p-distance", "The fraction of comparable positions at which two aligned sequences differ "
                           "(0.01 = 1 difference per 100 positions)."),
    ("Read depth", "How many sequencing reads covered a position; more reads mean a more reliable measurement."),
    ("Allele frequency", "In a pooled sample, the share of reads at a position that carry the alternative base."),
    ("PCoA", "Principal coordinates analysis: places samples on a plot so that samples with similar nucleotide "
             "frequencies sit close together."),
    ("Geographic precision", "How exactly a location is known: coordinates, a city, or only a state or country. "
                             "Only precise enough locations are used for kilometer distances."),
    ("Shared variant", "The same reference position and nucleotide change seen in both an individual sequence "
                       "and a pooled sample."),
    ("Unresolved evidence", "A variant that cannot be compared in one source because that source has no usable "
                            "coverage at the position. Unknown is not the same as absent."),
]


def count(value):
    return int(value) if value is not None and not pd.isna(value) else 0


def interpretation(result, name):
    """One short sentence about the actual results shown in a figure, or ''."""

    tables, summary = result["tables"], result["summary"]
    if name in {"cross_source", "evidence_matrix"} and result.get("comparison") is not None:
        counts = summary.get("category_counts", {})
        total = sum(counts.values())
        return (f"{summary.get('shared_variants', 0)} of {total} variants were observed in both sources; "
                "the rest are source-specific, unresolved, or incomparable.")
    if name == "genome_positions" and "variant_positions.csv" in tables:
        positions = tables["variant_positions.csv"]
        if positions.empty:
            return ""
        top = positions.sort_values(["Pooled_Samples", "Individual_Carriers"], ascending=False).iloc[0]
        return (f"{len(positions)} distinct variants; the most widespread is {top['Ref']}>{top['Alt']} at "
                f"position {top['Position']} ({count(top['Pooled_Samples'])} pooled samples, "
                f"{count(top['Individual_Carriers'])} individual sequences).")
    if name == "trajectories" and "temporal_trends.csv" in tables:
        trends = tables["temporal_trends.csv"]
        changes = trends["Change_In_Callable_Frequency"].dropna()
        if changes.empty:
            return ""
        row = trends.loc[changes.abs().idxmax()]
        return (f"The largest change between the first and last known frequency is "
                f"{row['Change_In_Callable_Frequency']:+.2f} for {row['Variant_Key']} at {row['Location_Name']}"
                + ("." if row["Complete_Series"] else "; some samples in that series have no known frequency."))
    if name in {"quality", "qc_outcomes"} and "pooled_samples_included" in summary:
        return (f"{summary['pooled_samples_included']} of {summary['pooled_samples']} pooled samples and "
                f"{summary.get('pooled_records_included', '—')} of {summary.get('pooled_records', '—')} "
                "pooled variant records were kept after quality control.")
    if name == "callable_positions" and "individual_mapped" in summary:
        return (f"{summary['individual_mapped']} of {summary['individual_samples']} individual sequences "
                "were mapped to the reference.")
    if name == "pcoa" and result.get("population") and result["population"].get("available"):
        explained = result["population"]["pcoa"]["explained"]
        if len(explained) >= 2:
            return (f"Axes 1 and 2 summarize {explained[0] * 100:.1f}% and {explained[1] * 100:.1f}% of the "
                    f"positive eigenvalues for {len(result['population']['coordinates'])} pooled samples.")
    if name == "genetic_distances" and "pairwise_distances.csv" in tables:
        pairs = tables["pairwise_distances.csv"]["P_Distance"].dropna()
        if not pairs.empty:
            return f"{len(pairs)} sequence pairs; p-distances range from {pairs.min():.3f} to {pairs.max():.3f}."
    if name == "variants_per_sample" and "variants_per_sample.csv" in tables:
        per_sample = tables["variants_per_sample.csv"]
        kept = per_sample[~per_sample["Quality_Status"].astype(str).str.startswith("Excluded")]
        values = kept["Distinct_SNVs_After_QC"].dropna()
        if not values.empty:
            return (f"Pooled samples kept after QC carry {int(values.min())} to {int(values.max())} "
                    "distinct SNVs.")
    if name in {"heatmap", "observation_status"} and "pooled_observation_matrix.csv" in tables:
        status = tables["pooled_observation_matrix.csv"]["Status"].value_counts()
        total = int(status.sum())
        if total:
            unknown = total - int(status.get("Variant observed", 0)) - int(
                status.get("Adequate depth; variant not reported", 0))
            return (f"Of {total} sample-variant cells, {int(status.get('Variant observed', 0))} have a reported "
                    f"frequency and {unknown} are low, unknown, or excluded coverage (not zeros).")
    annotation = result.get("annotation") or {}
    if name == "mutation_annotation_overview" and annotation.get("available"):
        return (f"Of {annotation['evaluated']} distinct SNVs, {annotation['missense']} are predicted missense, "
                f"{annotation['synonymous']} synonymous, {annotation['other_coding']} other coding changes, and "
                f"{annotation['noncoding']} noncoding; {annotation['evaluated'] - annotation['annotated']} could not "
                "be fully annotated.")
    return ""
