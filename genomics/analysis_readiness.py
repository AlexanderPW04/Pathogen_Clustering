"""Which Sprint 3 + 5 analyses the uploaded data actually support, and why.

Sprint 5 Phase 4A. Before the run, each method is given one status from the
mapped, validated inputs (never from a file or paper name):

- **Available**: the required data and validation checks are satisfied.
- **Partially available**: some scientifically meaningful calculation is
  possible, with a stated limitation.
- **Unavailable**: required information is missing or invalid.
- **Experimental**: the method exists but has not met its validation
  criteria (the Phase 3 mutation-frequency grouping).

The run itself still checks every record; this assessment only predicts
which methods have the inputs they need. An unavailable method never stops
the others.
"""

import pandas as pd

from genomics.table_layouts import BASE_COUNT_LAYOUT, WIDE_LAYOUT

AVAILABLE = "Available"
PARTIAL = "Partially available"
UNAVAILABLE = "Unavailable"
EXPERIMENTAL = "Experimental"
STATUSES = [AVAILABLE, PARTIAL, UNAVAILABLE, EXPERIMENTAL]

# What the documented formats carry (their own documentation, not guesses).
FORMAT_FIELDS = {
    "Fontenele et al. (2021) Supplementary Table 2": {
        "depth": "all", "alt": "reconstructable", "ref": False, "frequency": True, "dates": True},
    "iVar variants table": {"depth": "all", "alt": "measured", "ref": True, "frequency": True, "dates": False},
    "VCF": {"depth": "vcf", "alt": "vcf", "ref": "vcf", "frequency": True, "dates": False},
}


def filled_share(entry, field):
    """(non-blank values, rows) of a mapped field in a table entry; (0, rows) when absent."""

    table = entry.get("table")
    column = entry.get("columns", {}).get(field)
    if table is None:
        return 0, 0
    if column is None:
        return 0, len(table)
    values = table[column].fillna("").astype(str).str.strip()
    return int(values.ne("").sum()), len(table)


def variant_capabilities(review):
    """What the pooled variant files contain: depth, counts, frequencies (summed over files)."""

    entries = review["by_role"].get("snv", [])
    capability = {"files": len(entries), "depth": 0, "rows": 0, "alt": "none", "ref": False,
                  "frequency": True, "derived_frequency": False, "vcf": False, "formats": []}
    alt_kinds = set()
    for entry in entries:
        info = entry["info"]
        known = FORMAT_FIELDS.get(info.get("format")) or (FORMAT_FIELDS["VCF"] if info["kind"] == "VCF variant calls"
                                                          else None)
        if known is not None:
            capability["formats"].append(info.get("format") or info["kind"])
            capability["rows"] += info.get("records", 0)
            if known["depth"] == "all":
                capability["depth"] += info.get("records", 0)
            else:
                capability["vcf"] = True
            alt_kinds.add(known["alt"])
            capability["ref"] = capability["ref"] or known["ref"] is True
            continue
        depth, rows = filled_share(entry, "Read_Depth")
        capability["depth"] += depth
        capability["rows"] += rows
        alt, _ = filled_share(entry, "Alt_Depth")
        if info.get("layout") == BASE_COUNT_LAYOUT:
            alt_kinds.add("measured")
            capability["ref"] = True
            capability["derived_frequency"] = True
        elif alt:
            alt_kinds.add("measured" if alt == rows else "partly measured")
        elif depth:
            alt_kinds.add("reconstructable")
        if filled_share(entry, "Ref_Depth")[0]:
            capability["ref"] = True
    for kind in ("measured", "vcf", "partly measured", "reconstructable"):
        if kind in alt_kinds:
            capability["alt"] = kind
            break
    if capability["derived_frequency"] and len(entries) == sum(
            1 for entry in entries if entry["info"].get("layout") == BASE_COUNT_LAYOUT):
        capability["frequency"] = False
    return capability


def replicate_candidates(review):
    """Number of specimens with several samples and distinct documented Library_IDs (from the matching report)."""

    for row in review.get("matching", []):
        if row["Check"] == "Specimens with several sequenced samples":
            text = row["Detail"]
            return int(text.split(" have ")[0]) if text.split(" have ")[0].isdigit() else 0, row["Count"]
    return 0, 0


def analysis_readiness(review, options=None):
    """Rows {Analysis, Group, Status, Reason, Needs} for every Sprint 3 + 5 method (and Sprint 2 when relevant)."""

    options = options or {}
    roles = review["by_role"]
    rows = []

    def add(group, name, status, reason, needs):
        rows.append({"Group": group, "Analysis": name, "Status": status, "Reason": reason, "Needs": needs})

    has_snv = "snv" in roles
    reference = review.get("reference")
    pooled = next((row for row in review["completeness"] if row["Source"] == "Pooled samples"), {})
    samples = pooled.get("Samples", 0)
    capability = variant_capabilities(review)
    positions = roles.get("positions", [])
    counts_table = any(all(entry.get("columns", {}).get(f"{base}_Count") for base in "ACGT") for entry in positions)
    depth_table = bool(positions)
    dates = pooled.get("With_Collection_Date", 0)
    dates = dates if isinstance(dates, int) else 0

    if not has_snv:
        add("Pooled sequencing (Sprint 3 + 5)", "All pooled-sample analyses", UNAVAILABLE,
            "No pooled variant table, base-count table, or VCF was assigned.",
            "A processed SNV table (long or wide), iVar output, VCF, or A/C/G/T table")
    else:
        group = "Pooled sequencing (Sprint 3 + 5)"
        layouts = {entry["info"].get("layout") for entry in roles["snv"]}
        layout_note = (" (rewritten from a wide sample-by-mutation matrix)" if WIDE_LAYOUT in layouts else
                       " (derived from A/C/G/T base counts)" if BASE_COUNT_LAYOUT in layouts else "")
        add(group, "Basic SNV summaries", AVAILABLE,
            f"{samples} sample(s) with variant records{layout_note}.", "Sample, position, Ref, Alt")
        rows_total, depth = capability["rows"], capability["depth"]
        if capability["vcf"] and not depth:
            add(group, "Sequencing-depth QC", PARTIAL, "Depth comes from the VCF DP/AD fields and is checked "
                "record by record during the run.", "A depth value for each record")
        elif depth and depth == rows_total:
            add(group, "Sequencing-depth QC", AVAILABLE, "Every record has a read depth.", "A depth value per record")
        elif depth:
            add(group, "Sequencing-depth QC", PARTIAL, f"{depth:,} of {rows_total:,} records have a read depth; "
                "the others are 'not evaluable' (never passed).", "A depth value per record")
        else:
            add(group, "Sequencing-depth QC", UNAVAILABLE, "No read depth was provided, so the depth filter cannot "
                "be evaluated (records are 'not evaluable', never passed).", "A depth column")
        if capability["frequency"]:
            add(group, "Reported VAF summaries", AVAILABLE, "Reported alternative-allele frequencies are present.",
                "A frequency column (fraction or confirmed percentage)")
        else:
            add(group, "Reported VAF summaries", PARTIAL, "The source reported no frequencies; the summaries use "
                "frequencies calculated from the base counts and are labelled as calculated.", "Reported frequencies")
        alt = capability["alt"]
        if alt == "measured":
            add(group, "Count-based VAF calculations", AVAILABLE, "Alternative read counts and depth were provided "
                "(or derived from measured A/C/G/T counts).", "Alternative reads and depth")
            uncertainty = (AVAILABLE, "Measured counts give exact binomial (Wilson) intervals.")
        elif alt in {"vcf", "partly measured"}:
            add(group, "Count-based VAF calculations", PARTIAL, "Read counts are present for some records only "
                "(checked during the run); the others use reported frequencies.", "Alternative reads and depth")
            uncertainty = (PARTIAL, "Intervals only for records with counts.")
        elif alt == "reconstructable":
            add(group, "Count-based VAF calculations", PARTIAL, "Depth is available, but alternative read counts were "
                "not directly provided. Counts are reconstructed only where the reported frequency's precision "
                "makes them unique, and are labelled approximate.", "Alternative read counts")
            uncertainty = (PARTIAL, "Intervals use reconstructed (approximate) counts where reconstruction is unique.")
        else:
            add(group, "Count-based VAF calculations", UNAVAILABLE, "Neither read counts nor depth were provided, so "
                "no frequency can be calculated from reads.", "Alternative reads and depth")
            uncertainty = (UNAVAILABLE, "Without counts or depth there is no read-sampling uncertainty to estimate.")
        add(group, "Wilson uncertainty intervals", uncertainty[0], uncertainty[1], "Alternative reads and depth")
        sensitivity = (UNAVAILABLE if uncertainty[0] == UNAVAILABLE else uncertainty[0])
        add(group, "Beta-binomial sensitivity analysis", sensitivity,
            ("Hypothetical overdispersion scenarios (Advanced settings) widen the read-sampling intervals; they are "
             "scenarios, not estimates." if sensitivity != UNAVAILABLE else uncertainty[1]),
            "Counts and depth; scenario values")
        documented, shared = replicate_candidates(review)
        if documented and uncertainty[0] != UNAVAILABLE:
            add(group, "Empirical beta-binomial calibration", PARTIAL if uncertainty[0] == PARTIAL else AVAILABLE,
                f"{documented} specimen(s) have several samples with distinct documented Library_IDs; pairs are "
                "verified with the Phase 3 rules (same date and site, distinct libraries) during the run.",
                "Technical replicates: shared Specimen_ID, distinct Library_ID, counts")
        else:
            add(group, "Empirical beta-binomial calibration", UNAVAILABLE,
                ("Technical replicates were not verified, so empirical overdispersion cannot be estimated"
                 + (f" ({shared} specimen(s) share samples but lack distinct documented Library_IDs)." if shared
                    else " (no shared Specimen_ID with distinct Library_IDs).")),
                "Technical replicates: shared Specimen_ID, distinct Library_ID, counts")
        if dates >= 2:
            add(group, "Temporal mutation-frequency analysis", AVAILABLE if dates == samples else PARTIAL,
                f"{dates} of {samples} samples have a usable collection date"
                + ("" if dates == samples else "; undated samples are left out of time series") + "."
                + ("" if counts_table else " Without A/C/G/T counts, time series use detections only (an unreported "
                   "mutation is unknown, not 0%)."), "Collection dates for two or more samples")
        else:
            add(group, "Temporal mutation-frequency analysis", UNAVAILABLE,
                f"{dates} of {samples} samples have a usable collection date (missing, partial, or ambiguous dates "
                "are not guessed).", "Collection dates for two or more samples")
        if counts_table:
            add(group, "Full-count Yue-Clayton dissimilarity", AVAILABLE, "A/C/G/T counts were supplied in a "
                "per-position table.", "Per-position A/C/G/T counts for every sample")
        elif BASE_COUNT_LAYOUT in layouts:
            add(group, "Full-count Yue-Clayton dissimilarity", UNAVAILABLE, "Base counts exist only at the listed "
                "variable positions; full nucleotide coverage is unavailable, so exact full-count Yue-Clayton "
                "analysis is unsupported.", "Per-position A/C/G/T counts across the genome")
        else:
            add(group, "Full-count Yue-Clayton dissimilarity", UNAVAILABLE, "Full nucleotide coverage is unavailable, "
                "so exact full-count Yue-Clayton analysis is unsupported.", "Per-position A/C/G/T counts")
        if samples < 3:
            add(group, "PCoA", UNAVAILABLE, f"{samples} sample(s): an ordination needs at least three.",
                "Three or more samples")
        elif counts_table:
            add(group, "PCoA", AVAILABLE, "PCoA of Yue-Clayton dissimilarities (the 2021 study's method).",
                "A/C/G/T counts")
        else:
            add(group, "PCoA", PARTIAL, "Only an exploratory PCoA of Jaccard distances between reported SNV sets "
                "(not the study's method; an unreported SNV is not evidence of absence).", "A/C/G/T counts")
        add(group, "Coverage-aware non-detection", AVAILABLE if depth_table else UNAVAILABLE,
            "A per-position table classifies unreported variants as covered or not." if depth_table else
            "No per-position depth: unreported variants stay 'coverage unknown' (never 0%).",
            "Per-position depth or A/C/G/T counts")
        declared = [entry for entry in roles["snv"] if entry.get("declared_reference")]
        if reference:
            add(group, "Reference-based comparisons", AVAILABLE,
                f"Ref bases are checked against {reference['id']} ({reference['length']:,} bases)"
                + ("; the reference label was stated by you" if declared else "") + ".",
                "A reference FASTA with the same identifier")
        else:
            add(group, "Reference-based comparisons", PARTIAL,
                "Positions are compared by reference label only; without the reference FASTA, Ref bases are not "
                "verified.", "A reference FASTA")
        mode = options.get("grouping_mode", "profiles")
        add(group, "Experimental mutation-frequency grouping", EXPERIMENTAL,
            ("Mutation grouping is experimental and disabled by default; profiles are shown without groups."
             if mode != "experimental" else
             "You enabled the experimental view: groups are shown, labelled 'not validated'.")
            + (" Counts are available, so the depth-aware variant can run." if alt == "measured" else
               " Without measured counts only ordinary clustering of frequencies is possible."
               if mode == "experimental" else ""),
            "Advanced settings: 'Also show experimental groups'")
    annotation = roles.get("annotation")
    variant_refs = set()
    for entry in roles.get("snv", []):
        variant_refs |= set(entry["info"].get("reference_ids", []))
    if annotation:
        matched = bool(set(annotation[0]["info"]["reference_ids"]) & (variant_refs | (
            {reference["id"]} if reference else set())))
        add("Annotation (Sprint 5 Phase 4)", "Gene and protein annotation", AVAILABLE if matched else UNAVAILABLE,
            "The annotation describes the variants' reference (same accession.version)." if matched else
            "The annotation describes a different reference; annotations are never transferred.",
            "GenBank or GFF3 of the same accession.version")
    else:
        add("Annotation (Sprint 5 Phase 4)", "Gene and protein annotation", UNAVAILABLE,
            "Optional: no GenBank or GFF3 annotation was uploaded.", "GenBank or GFF3 of the same accession.version")
    plan = {row["Analysis"]: row for row in review.get("analyses", [])}
    cross = plan.get("Cross-source comparison")
    if cross:
        add("Individual + pooled (Sprint 2 + 3)", "Sprint 2 + 3 cross-source comparisons",
            AVAILABLE if cross["Status"] == "Will run" else UNAVAILABLE, cross["Reason"], cross["Needs"])
    if "fasta" in roles:
        translate = {"Will run": AVAILABLE, "Limited": PARTIAL, "Unavailable": UNAVAILABLE}
        for name in ("Sequence comparison (MAFFT, p-distance)", "Temporal screening of sequence pairs",
                     "Geographic-distance screening", "Reference-coordinate variant profiles"):
            if name in plan:
                row = plan[name]
                add("Individual sequences (Sprint 2)", name, translate[row["Status"]], row["Reason"], row["Needs"])
    return rows


def readiness_counts(rows):
    return {status: sum(1 for row in rows if row["Status"] == status) for status in STATUSES}


def readiness_table(rows):
    return pd.DataFrame(rows, columns=["Group", "Analysis", "Status", "Reason", "Needs"])
