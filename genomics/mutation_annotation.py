"""Sprint 5, Phase 4: mutation annotation joined to the existing sequencing evidence.

This step runs after the Sprint 3 and Sprint 5 Phase 1-3 analyses. It reads
their tables and never changes them:

- each distinct SNV (pooled and individual) is annotated once
  (`snv_consequences.annotate_variants`);
- the existing Phase 1 evidence (`sprint5_vaf_observations.csv`) and temporal
  summary (`sprint5_temporal_summary.csv`) are summarized per Variant_Key; VAF,
  depth, and intervals are copied or summarized, never recalculated;
- any annotation published with the input (the Fontenele et al. table's codons
  or iVar's GFF columns) is compared with the prediction as an independent check.

Outputs:

- `sprint5_mutation_annotations.csv`: one row per SNV per overlapping CDS;
- `sprint5_annotated_mutations.csv`: one row per SNV with its predicted consequence
  and its sequencing evidence;
- `sprint5_annotated_observations.csv`: the observed pooled records with their
  annotation (sample, date, VAF, depth, count basis, intervals, QC);
- `sprint5_annotation_sources.csv`: every coding feature, its checks, and status;
- `sprint5_annotation_comparison.csv`: agreement with annotation supplied in the input.

Predicted consequences describe sequence changes only. They do not show effects
on transmissibility, severity, immune escape, drug resistance, or fitness.
"""

import pandas as pd

from genomics import figures
from genomics.genome_annotation import SUPPORTED, load_annotation, source_table
from genomics.snv_consequences import (
    ANNOTATED,
    CATEGORY_ORDER,
    MISSENSE,
    NONCODING,
    NOT_ANNOTATED,
    OTHER_CODING,
    OVERLAP_DIFFERS,
    PARTLY_ANNOTATED,
    SYNONYMOUS,
    UNSUPPORTED_CONSEQUENCE,
    annotate_variants,
)


SECTION_TITLE = "Mutation Annotation and Protein Consequences"
EXPLANATION = ("This panel shows where each observed mutation lies in the pathogen genome, which annotated gene or "
               "coding region it affects, and whether that single nucleotide change is predicted to change the "
               "encoded amino acid. Predictions come from the uploaded reference annotation, and each one is "
               "shown with the strength of its sequencing evidence.")
MISSING_ANNOTATION = ("No genome annotation (GenBank or GFF3) was uploaded, so coding consequences cannot be "
                      "determined. Every other analysis ran normally. To add this panel, upload the annotation "
                      "that matches the reference used for variant calling (for SARS-CoV-2, NCBI NC_045512.2).")
LIMITATIONS = [
    "Predicted consequences describe the change in the encoded amino acid only. A synonymous change can still "
    "have biological effects (for example on RNA structure or codon usage), and a missense change does not by "
    "itself increase transmissibility, severity, immune or vaccine escape, drug resistance, host range, "
    "replication, or outbreak risk. Those claims need independent experimental or epidemiological evidence.",
    "Each SNV is annotated alone against the reference. Pooled sequencing reads do not show which mutations "
    "occur together on the same viral genome, so combined effects of nearby mutations are not predicted; codons "
    "changed by more than one observed SNV are flagged instead.",
    "Annotations are applied only when the variant's reference is the annotated accession.version and the "
    "reference base agrees. Records on other references are not annotated, never transferred.",
    "Coding structures that cannot be interpreted reliably (partial genes, RNA editing, translation exceptions, "
    "parts out of order, or a translation that disagrees with the record) are marked unsupported. Bases read "
    "twice because of programmed ribosomal frameshifting get no single consequence.",
    "Insertions, deletions, and multi-base changes are not annotated in this phase.",
    "Annotation does not identify lineages and does not validate the exploratory mutation-frequency groups.",
]
TABLE_ROWS = 25
OBSERVATION_COLUMNS = [
    "Sample_ID", "Collection_Date", "Location_Name", "Variant_Key", "Reported_Frequency",
    "Count_Derived_Frequency", "VAF", "Alt_Read_Count", "Total_Depth", "Alt_Count_Source",
    "Wilson_Lower", "Wilson_Upper", "Interval_Status", "BetaBinomial_Lower", "BetaBinomial_Upper",
    "Reference_Status", "QC_Status", "Analysis_Included", "Counted_Once",
]
COMPARISON_COLUMNS = ["Variant_Key", "Source_Of_Input_Annotation", "Input_Gene", "Input_Ref_Codon",
                      "Input_Alt_Codon", "Input_Amino_Acids", "Input_Syn_NonSyn", "Predicted_Genes",
                      "Predicted_Codons", "Predicted_Amino_Acids", "Predicted_Category", "Agreement", "Note"]


def distinct_variants(pooled, individual):
    """Distinct SNVs from both sources with where each was seen."""

    frames = []
    for name, data in [("Pooled", pooled["variants"] if pooled else None),
                       ("Individual", individual["variants"] if individual else None)]:
        if data is not None and len(data):
            frame = data[["Variant_Key", "Reference_ID", "Position", "Ref", "Alt"]].drop_duplicates("Variant_Key")
            frames.append(frame.assign(Source=name))
    if not frames:
        return pd.DataFrame(columns=["Variant_Key", "Reference_ID", "Position", "Ref", "Alt", "Evidence_Sources"])
    combined = pd.concat(frames, ignore_index=True)
    sources = combined.groupby("Variant_Key", sort=False)["Source"].agg(lambda values: " + ".join(sorted(set(values))))
    unique = combined.drop_duplicates("Variant_Key").drop(columns="Source")
    unique["Evidence_Sources"] = unique["Variant_Key"].map(sources)
    return unique.sort_values(["Reference_ID", "Position", "Variant_Key"]).reset_index(drop=True)


def joined_counts(frame, column):
    """Per Variant_Key: "value (count); ..." with the most frequent value first (one pass, no per-group lambdas)."""

    counts = frame.groupby(["Variant_Key", column], sort=False).size().reset_index(name="_n")
    counts = counts.sort_values(["Variant_Key", "_n"], ascending=[True, False], kind="stable")
    text = {}
    for key, value, count in counts.itertuples(index=False):
        text.setdefault(key, []).append(f"{value} ({count})")
    return pd.Series({key: "; ".join(items) for key, items in text.items()}, dtype=object)


def evidence_summary(evidence):
    """Per Variant_Key: how strongly the pooled reads support it (summaries of the Phase 1-2 values)."""

    used = evidence[evidence["Counted_Once"] & evidence["Analysis_Included"]]
    grouped = used.groupby("Variant_Key", sort=False)
    summary = pd.DataFrame({
        "Pooled_Samples_Observed": grouped["Sample_ID"].nunique(),
        "First_Collection_Date": grouped["Collection_Date"].min(),
        "Last_Collection_Date": grouped["Collection_Date"].max(),
        "Median_VAF": grouped["VAF"].median(),
        "Max_VAF": grouped["VAF"].max(),
        "Median_Depth": grouped["Total_Depth"].median(),
        "Min_Depth": grouped["Total_Depth"].min(),
        "Max_Depth": grouped["Total_Depth"].max(),
        "Observations_With_Wilson_Interval": grouped["Wilson_Lower"].count(),
        # Phase 2 fills beta-binomial bounds only for replicate-calibrated observations.
        "Observations_With_Calibrated_BetaBinomial": (grouped["BetaBinomial_Lower"].count()
                                                      if "BetaBinomial_Lower" in used else 0),
    })
    summary["Count_Basis"] = joined_counts(used, "Alt_Count_Source")
    statuses = used[["Variant_Key", "Reference_Status"]].drop_duplicates()
    summary["Reference_Status"] = statuses.groupby("Variant_Key", sort=False)["Reference_Status"].agg(
        lambda values: "; ".join(sorted(map(str, values))))
    best = used.dropna(subset=["VAF"]).sort_values("VAF", ascending=False).drop_duplicates("Variant_Key")
    best = best.set_index("Variant_Key")
    interval = best["Wilson_Lower"].map("{:.3g}".format) + "-" + best["Wilson_Upper"].map("{:.3g}".format)
    summary["Max_VAF_Wilson_95"] = interval.where(best["Wilson_Lower"].notna(), "").reindex(summary.index).fillna("")
    excluded = evidence[~evidence["Analysis_Included"]].groupby("Variant_Key")["Sample_ID"].size()
    summary = summary.reindex(summary.index.union(excluded.index))
    summary["Observations_Excluded_By_QC"] = excluded.reindex(summary.index).fillna(0).astype(int)
    summary["Pooled_Samples_Observed"] = summary["Pooled_Samples_Observed"].fillna(0).astype(int)
    return summary


def temporal_summary(series):
    """Per Variant_Key: how many site series exist and how many have known values on 2+ dates."""

    if series is None or series.empty:
        return pd.DataFrame(index=pd.Index([], name="Variant_Key"))
    grouped = series.groupby("Variant_Key", sort=False)
    largest = (series.assign(_size=series["Absolute_Change"].abs()).sort_values("_size", ascending=False)
               .drop_duplicates("Variant_Key").set_index("Variant_Key")["Absolute_Change"])
    sizes = grouped.size()
    return pd.DataFrame({
        "Site_Series": sizes,
        "Series_With_Two_Or_More_Dates": (series["Known_Dated_Time_Points"] >= 2).groupby(
            series["Variant_Key"], sort=False).sum().astype(int),
        "Largest_Absolute_Change": largest.reindex(sizes.index),
        "Series_Status": joined_counts(series, "Series_Status").reindex(sizes.index),
    })


def input_annotation_comparison(variants, consequences, summary):
    """Agreement with codons/amino acids supplied in the input (published or caller-generated), if any."""

    sources = [("Published columns", "Published_Annotation", "Published_Ref_Codon",
                "Published_Alt_Codon", "Published_Ref_AA", "Published_Alt_AA", "Published_Syn_NonSyn"),
               ("iVar GFF columns", "Source_Gene_Feature", "Source_Ref_Codon",
                "Source_Alt_Codon", "Source_Ref_AA", "Source_Alt_AA", None)]
    source = next((item for item in sources if item[2] in variants.columns and variants[item[2]].notna().any()),
                  None)
    if source is None:
        return None
    label, gene_column, ref_codon, alt_codon, ref_aa, alt_aa, syn = source
    columns = ["Variant_Key"] + [name for name in [gene_column, ref_codon, alt_codon, ref_aa, alt_aa, syn]
                                 if name and name in variants.columns]
    first = {row["Variant_Key"]: row for row in
             variants[columns].drop_duplicates("Variant_Key").to_dict("records")}
    coding = {}
    for row in consequences.loc[consequences["Codon_Number"].notna(),
                                ["Variant_Key", "Gene", "CDS_ID", "Ref_Codon", "Alt_Codon",
                                 "Amino_Acid_Change"]].to_dict("records"):
        coding.setdefault(row["Variant_Key"], []).append(row)
    category = dict(zip(summary["Variant_Key"], summary["Consequence_Category"]))

    def text(item, column):
        value = item.get(column, "") if column else ""
        return "" if pd.isna(value) or str(value).strip().upper() in {"", "NA", "NAN"} else str(value).strip()

    rows = []
    for key in summary["Variant_Key"]:
        item = first.get(key)
        if item is None:
            continue
        gene, codon_1, codon_2 = text(item, gene_column), text(item, ref_codon).upper(), text(item, alt_codon).upper()
        amino = f"{text(item, ref_aa)}>{text(item, alt_aa)}" if text(item, ref_aa) else ""
        group = coding.get(key)
        note = ""
        if not codon_1:
            agreement = ("Agrees (no codon in either)" if group is None else
                         "Differs (input gives no codon; a coding consequence is predicted)")
            if syn and text(item, syn).casefold().startswith("syn") and category.get(key) == NONCODING:
                note = ("The input labels this noncoding SNV 'Syn'; this application reports it as noncoding, "
                        "not synonymous.")
        elif group is None:
            agreement = "Differs (input gives a codon; no supported coding consequence is predicted)"
        else:
            names = {row["Gene"].casefold() for row in group} | {row["CDS_ID"].casefold() for row in group}
            gene_ok = not gene or gene.casefold() in names or any(
                row["CDS_ID"] and row["CDS_ID"].casefold() in gene.casefold() for row in group)
            match = any(row["Ref_Codon"] == codon_1 and row["Alt_Codon"] == codon_2 and (
                not amino or (row["Amino_Acid_Change"] and
                              f"{row['Amino_Acid_Change'][0]}>{row['Amino_Acid_Change'][-1]}" == amino))
                        for row in group)
            agreement = "Agrees" if match and gene_ok else "Differs"
            if agreement == "Differs":
                note = "Codon, amino acid, or gene differs from the input; check the reference and annotation."
        rows.append({"Variant_Key": key, "Source_Of_Input_Annotation": label, "Input_Gene": gene,
                     "Input_Ref_Codon": codon_1, "Input_Alt_Codon": codon_2, "Input_Amino_Acids": amino,
                     "Input_Syn_NonSyn": text(item, syn) if syn else "",
                     "Predicted_Genes": "" if group is None else "; ".join(dict.fromkeys(r["Gene"] for r in group)),
                     "Predicted_Codons": "" if group is None else "; ".join(
                         f"{r['Ref_Codon']}>{r['Alt_Codon']}" for r in group),
                     "Predicted_Amino_Acids": "" if group is None else "; ".join(
                         f"{r['Amino_Acid_Change'][0]}>{r['Amino_Acid_Change'][-1]}" for r in group
                         if r["Amino_Acid_Change"]),
                     "Predicted_Category": category.get(key, ""), "Agreement": agreement, "Note": note})
    return pd.DataFrame(rows, columns=COMPARISON_COLUMNS)


def category_counts(summary):
    counts = summary["Consequence_Category"].value_counts()
    return {name: int(counts.get(name, 0)) for name in CATEGORY_ORDER}


def annotated_mutation_table(variants, summary, evidence, series, individual, groups):
    """One row per SNV: predicted consequence plus its sequencing evidence (existing values only)."""

    table = variants[["Variant_Key", "Evidence_Sources"]].merge(summary, on="Variant_Key", how="left")
    if evidence is not None and len(evidence):
        table = table.join(evidence_summary(evidence), on="Variant_Key")
    if series is not None and len(series):
        table = table.join(temporal_summary(series), on="Variant_Key")
    if individual is not None and len(individual["variants"]):
        carriers = individual["variants"].groupby("Variant_Key")["Sample_ID"].nunique()
        table["Individual_Sequences_Carrying"] = table["Variant_Key"].map(carriers).fillna(0).astype(int)
    if groups is not None and len(groups) and "Group_ID" in groups:
        grouped = groups[groups["Group_ID"].astype(str).ne("")].drop_duplicates("Variant_Key")
        label = grouped.set_index("Variant_Key")["Group_ID"].astype(str) + " (experimental, not validated)"
        table["Exploratory_Group"] = table["Variant_Key"].map(label).fillna("")
    return table


def run_mutation_annotation(result, annotation_path, reference, pooled, individual, add_figure, note,
                            label=None, annotation=None, genetic_code=None):
    """Add the Phase 4 tables, one figure, and a view summary to `result`; returns the summary.

    `label` is the uploaded file name shown to the user (the server copy has a generic name).
    """

    if not annotation_path:
        info = {"available": False, "title": SECTION_TITLE, "explanation": EXPLANATION,
                "message": MISSING_ANNOTATION, "limitations": LIMITATIONS}
        result["annotation"] = info
        return info
    if annotation is None:  # the caller may already have loaded it (for the gene track)
        annotation = load_annotation(annotation_path, reference, genetic_code)
    if label:
        annotation["source_file"] = label
    for text in annotation["notes"]:
        note("annotation", text)
    tables = result["tables"]
    sources = source_table(annotation)
    tables["sprint5_annotation_sources.csv"] = sources
    variants = distinct_variants(pooled, individual)
    evidence = tables.get("sprint5_vaf_observations.csv")
    samples_by_variant = None
    if evidence is not None and len(evidence):
        reported = evidence[evidence["Counted_Once"] & evidence["Analysis_Included"]]
        samples_by_variant = reported.groupby("Variant_Key")["Sample_ID"].agg(set).to_dict()
    consequences, summary = annotate_variants(variants, annotation, samples_by_variant)
    tables["sprint5_mutation_annotations.csv"] = consequences
    groups = tables.get("sprint5_mutation_groups.csv") if (result.get("depth_aware") or {}).get(
        "grouping", {}).get("experimental") else None
    mutations = annotated_mutation_table(variants, summary, evidence, tables.get("sprint5_temporal_summary.csv"),
                                         individual, groups)
    tables["sprint5_annotated_mutations.csv"] = mutations
    if evidence is not None and len(evidence):
        observed = evidence[[name for name in OBSERVATION_COLUMNS if name in evidence.columns]]
        lookup = summary.set_index("Variant_Key")
        observed = observed.assign(
            Gene=observed["Variant_Key"].map(lookup["Genes"]),
            Predicted_Amino_Acid_Change=observed["Variant_Key"].map(lookup["Amino_Acid_Changes"]),
            Consequence_Category=observed["Variant_Key"].map(lookup["Consequence_Category"]),
            Annotation_Status=observed["Variant_Key"].map(lookup["Annotation_Status"]))
        tables["sprint5_annotated_observations.csv"] = observed
    comparison = input_annotation_comparison(pooled["variants"], consequences, summary) if pooled else None
    if comparison is not None:
        tables["sprint5_annotation_comparison.csv"] = comparison

    counts = category_counts(summary)
    statuses = summary["Annotation_Status"].value_counts()
    supported = sources[sources["Status"].eq(SUPPORTED)]
    shared = int(summary["Shared_Codon_With"].ne("").sum())
    annotated_references = sorted(set(summary.loc[summary["Annotation_Status"].eq(ANNOTATED), "Reference_ID"]))
    if statuses.get(NOT_ANNOTATED, 0):
        reasons = consequences[consequences["Consequence"].eq(NOT_ANNOTATED)]["Reference_Check"].value_counts()
        note("annotation", f"{int(statuses[NOT_ANNOTATED])} SNV(s) were not annotated: "
             + "; ".join(f"{name} ({count})" for name, count in reasons.items()) + ".")
    if shared:
        note("annotation", f"{shared} SNV(s) share a codon with another observed SNV; their combined amino-acid "
                           "change is not predicted because pooled reads do not show linkage.")
    if len(sources) - len(supported):
        note("annotation", f"{len(sources) - len(supported)} coding feature(s) are unsupported and give no "
                           "amino-acid prediction (see sprint5_annotation_sources.csv).")

    figure_name = None
    plotted = consequences[consequences["Reference_ID"].isin(annotated_references)] if annotated_references \
        else consequences.iloc[0:0]
    if len(plotted):
        main = plotted["Reference_ID"].value_counts().index[0]
        cds_rows = [cds for cds in annotation["cds"] if cds["Reference_ID"] == main]
        length = annotation["records"][main]["Length"]
        other = sorted(set(plotted["Reference_ID"]) - {main})
        add_figure(result, "mutation_annotation_overview", "Where the mutations fall and what they are predicted "
                   "to change",
                   f"Top: each distinct SNV at its position on {main} (one row per predicted category), over the "
                   "annotated coding regions. Bottom: SNVs per coding region by predicted consequence; a base in "
                   "two overlapping coding regions is counted in both. Predictions describe sequence changes only."
                   + (" Other annotated references (" + ", ".join(other) + ") are in the tables only."
                      if other else ""),
                   lambda: figures.mutation_annotation_figure(plotted[plotted["Reference_ID"].eq(main)],
                                                              summary[summary["Reference_ID"].eq(main)],
                                                              cds_rows, length))
        figure_name = "mutation_annotation_overview"

    columns = ["Variant_Key", "Genes", "Mature_Peptides", "Nucleotide_Change", "Amino_Acid_Changes",
               "Consequence_Category", "Reference_Check", "Pooled_Samples_Observed", "Median_VAF",
               "Max_VAF", "Max_VAF_Wilson_95", "Min_Depth", "Max_Depth", "Count_Basis",
               "Series_With_Two_Or_More_Dates", "Individual_Sequences_Carrying", "Exploratory_Group"]
    sort_columns = [name for name in ["Pooled_Samples_Observed", "Individual_Sequences_Carrying"]
                    if name in mutations]
    shown = mutations.sort_values(sort_columns + ["Position"], ascending=[False] * len(sort_columns) + [True]) \
        if sort_columns else mutations.sort_values("Position")
    agreement = comparison["Agreement"].value_counts().to_dict() if comparison is not None else {}
    info = {
        "available": True, "title": SECTION_TITLE, "explanation": EXPLANATION, "limitations": LIMITATIONS,
        "format": annotation["format"], "source_file": annotation["source_file"],
        "references": [{"Reference_ID": key, "Organism": record["Organism"], "Segment": record["Segment"],
                        "Length": record["Length"], "Sequence_Check": record["Sequence_Check"]}
                       for key, record in annotation["records"].items()],
        "cds": len(sources), "cds_supported": len(supported),
        "evaluated": len(summary), "annotated": int(statuses.get(ANNOTATED, 0)),
        "partly": int(statuses.get(PARTLY_ANNOTATED, 0)), "not_annotated": int(statuses.get(NOT_ANNOTATED, 0)),
        "counts": counts, "synonymous": counts[SYNONYMOUS], "missense": counts[MISSENSE],
        "other_coding": sum(counts[name] for name in OTHER_CODING),
        "other_detail": {name: counts[name] for name in OTHER_CODING if counts[name]},
        "overlap_differs": counts[OVERLAP_DIFFERS], "noncoding": counts[NONCODING],
        "unsupported": counts[UNSUPPORTED_CONSEQUENCE], "unannotatable": counts[UNSUPPORTED_CONSEQUENCE]
        + counts[NOT_ANNOTATED], "shared_codons": shared, "agreement": agreement,
        "comparison_source": comparison["Source_Of_Input_Annotation"].iloc[0]
        if comparison is not None and len(comparison) else "",
        "figure": figure_name, "table_columns": [name for name in columns if name in shown],
        "table": shown, "annotated_references": annotated_references,
    }
    result["annotation"] = info
    return info
