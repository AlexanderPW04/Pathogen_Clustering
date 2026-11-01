"""Turn integrated analysis results into simple values for the dashboard page.

Nothing here calculates science; it only formats tables and figures that
`genomics.integrated_analysis` already produced.
"""

from base64 import b64encode

import pandas as pd


MISSING = "—"  # an em dash marks a value that is not available


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


def table_view(table, columns=None, limit=None, title=""):
    """Return {'title', 'columns', 'rows', 'total', 'shown'} for the template."""

    if table is None:
        return None
    data = table if columns is None else table[[name for name in columns if name in table.columns]]
    total = len(data)
    if limit is not None:
        data = data.head(limit)
    rows = [[format_value(value) for value in row] for row in data.itertuples(index=False)]
    return {"title": title, "columns": list(data.columns), "rows": rows,
            "total": total, "shown": len(rows)}


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
    }


def build_dashboard(result, download_id):
    """Organize results into the dashboard's eight sections."""

    tables = result["tables"]
    notes = result["notes"]

    def figures_for(*names):
        return [view for view in (figure_view(result, name, download_id) for name in names) if view]

    view = {"summary": result["summary"], "notes": notes, "download_id": download_id}

    view["quality"] = {
        "figures": figures_for("quality", "variants_per_sample", "depth_vs_snvs"),
        "tables": [table for table in [
            table_view(tables.get("paper_method_reproduction.csv"),
                       title="Which Fontenele et al. (2021) methods were reproduced"),
            table_view(tables.get("variants_per_sample.csv"), title="Pooled samples: SNVs, coverage, and QC",
                       columns=["Sample_ID", "Collection_Date", "Location_Name", "SNV_Records",
                                "Distinct_SNVs", "Distinct_SNVs_After_QC", "Genome_Coverage_Percent",
                                "Mean_Depth", "Quality_Status"]),
            table_view(tables.get("individual_samples.csv"), title="Individual sequences: reference mapping",
                       columns=["Sample_ID", "Collection_Date", "Location_Name", "Sequenced_Span",
                                "Reference_Identity_Percent", "Genome_Coverage_Percent",
                                "Quality_Status", "Variant_Profile_Status"]),
        ] if table],
    }
    fit = result.get("regression")
    view["regression"] = fit

    positions = tables.get("variant_positions.csv")
    if positions is not None:
        positions = positions.sort_values(["Pooled_Samples", "Individual_Carriers", "Position"],
                                          ascending=[False, False, True])
    view["variation"] = {
        "figures": figures_for("genome_positions"),
        "tables": [table for table in [
            table_view(positions, limit=40, title="Variants by number of samples (first 40)"),
            table_view(tables.get("individual_indels.csv"),
                       title="Insertions and deletions in individual sequences (not compared)"),
        ] if table],
    }

    trends = tables.get("temporal_trends.csv")
    if trends is not None:
        trends = trends.sort_values(["Reported", "Variant_Key"], ascending=[False, True])
    view["temporal"] = {
        "figures": figures_for("trajectories", "heatmap"),
        "tables": [table for table in [
            table_view(trends, limit=30, title="Reported-frequency trends by site (first 30)"),
            table_view(tables.get("site_date_summary.csv"), title="Pooled samples by site and date"),
        ] if table],
    }

    population = result.get("population")
    view["population"] = None
    if population is not None:
        pcoa = population.get("pcoa")
        view["population"] = {
            "available": population["available"],
            "method": population["method"],
            "message": population["message"],
            "figures": figures_for("pcoa"),
            "eigen": None if pcoa is None else {
                "axis_1": pcoa["explained"][0] * 100 if len(pcoa["explained"]) > 0 else 0,
                "axis_2": pcoa["explained"][1] * 100 if len(pcoa["explained"]) > 1 else 0,
                "negative": pcoa["negative_eigenvalues"],
            },
            "tables": [table for table in [
                table_view(tables.get("pcoa_coordinates.csv"), title="PCoA coordinates"),
                table_view(tables.get("nucleotide_diversity.csv"),
                           title="Mean per-position nucleotide diversity (1 - sum of squared A/C/G/T frequencies)"),
            ] if table],
        }

    view["geography"] = {"tables": [table for table in [
        table_view(tables.get("geographic_context.csv"), title="Samples by documented location"),
    ] if table]}

    comparison = result.get("comparison")
    view["cross_source"] = None
    if comparison is not None:
        view["cross_source"] = {
            "figures": figures_for("cross_source", "cross_source_timeline"),
            "tables": [table for table in [
                table_view(comparison["counts"], title="Variants per evidence category"),
                table_view(comparison["table"], limit=60, title="Variant-by-variant comparison",
                           columns=["Variant_Key", "Category", "Individual_Carriers",
                                    "Individual_Reference_Calls", "Individual_Not_Callable",
                                    "Pooled_Reported", "Pooled_Adequate_Depth_Not_Reported",
                                    "Pooled_Low_Or_No_Coverage", "Pooled_Coverage_Unknown_Or_QC_Excluded",
                                    "Pooled_Frequency_Median", "First_Individual_Date",
                                    "First_Pooled_Date", "Days_Pooled_Minus_Individual",
                                    "Date_Pairs_Within_Window", "Individual_Carrier_Locations",
                                    "Pooled_Reporting_Locations", "Nearest_Km_Between_Carriers", "Note"]),
                table_view(tables.get("variant_signature_overlap.csv"),
                           title="Overlap with documented variant profiles (not lineage assignment)"),
            ] if table],
        }

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
