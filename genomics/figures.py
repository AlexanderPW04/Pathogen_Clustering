"""Scientific figures for the integrated dashboard (PNG bytes, no Flask code).

Colors come from a colorblind-checked categorical palette, assigned in a fixed
order, and a single-hue blue ramp for frequencies. Every figure answers one
surveillance question; its underlying numbers are also offered as a CSV.
"""

from io import BytesIO

import matplotlib

matplotlib.use("Agg")

from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.dates import AutoDateLocator, ConciseDateFormatter, date2num, num2date  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from genomics.unified_schema import EXCLUDED_BY_QC, VARIANT_OBSERVED  # noqa: E402


CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
               "#6250d6", "#e34948"]
OTHER_GRAY = "#8a8985"
TEXT = "#0b0b0b"
MUTED_TEXT = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"
NOT_REPORTED_FILL = "#ffffff"
LOW_COVERAGE_FILL = "#c9c8c3"
UNKNOWN_FILL = "#ecebe8"
BLUE_RAMP = LinearSegmentedColormap.from_list(
    "blue_ramp", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
)
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "p", "h", "*"]


def new_figure(width, height):
    figure = Figure(figsize=(width, height), layout="constrained", facecolor=SURFACE)
    return figure


def style_axis(axis):
    axis.set_facecolor(SURFACE)
    axis.grid(True, color=GRID, linewidth=0.8)
    axis.set_axisbelow(True)
    for side in ["top", "right"]:
        axis.spines[side].set_visible(False)
    for side in ["left", "bottom"]:
        axis.spines[side].set_color(MUTED_TEXT)
    axis.tick_params(colors=MUTED_TEXT, labelsize=8)
    axis.title.set_color(TEXT)
    axis.xaxis.label.set_color(TEXT)
    axis.yaxis.label.set_color(TEXT)


def to_png(figure):
    image = BytesIO()
    figure.savefig(image, format="png", dpi=110, facecolor=SURFACE)
    figure.clear()  # figures hold reference cycles; empty this one before the next is drawn
    return image.getvalue()


def location_colors(locations):
    """Fixed-order colors by sorted location name; extra sites fold into gray."""

    colors = {}
    for index, name in enumerate(sorted({str(value) for value in locations})):
        colors[name] = CATEGORICAL[index] if index < len(CATEGORICAL) else OTHER_GRAY
    return colors


def date_axis(axis):
    locator = AutoDateLocator()
    axis.xaxis.set_major_locator(locator)
    axis.xaxis.set_major_formatter(ConciseDateFormatter(locator))


def variants_per_sample_figure(per_sample):
    """Paper Fig 2A analogue: distinct SNVs per pooled sample, by site."""

    data = per_sample.sort_values(["Location_Name", "Collection_Date", "Sample_ID"]).reset_index(drop=True)
    colors = location_colors(data["Location_Name"])
    figure = new_figure(max(6, 0.35 * len(data) + 2), 3.8)
    axis = figure.subplots()
    style_axis(axis)
    for index, row in data.iterrows():
        included = bool(row["Analysis_Included"])
        height = row["Distinct_SNVs_After_QC"] if included else row["Distinct_SNVs"]
        axis.bar(index, height, width=0.7, color=colors[str(row["Location_Name"])],
                 alpha=1.0 if included else 0.35, hatch=None if included else "///",
                 edgecolor=SURFACE if included else MUTED_TEXT, linewidth=0.8)
    axis.set_xticks(range(len(data)))
    axis.set_xticklabels(data["Sample_ID"], rotation=60, ha="right", fontsize=7)
    axis.set_ylabel("Distinct SNVs (after QC)")
    axis.set_title("Distinct SNVs per pooled sample, colored by sampling site", fontsize=10)
    handles = [Patch(color=color, label=name) for name, color in colors.items()]
    if (~data["Analysis_Included"].astype(bool)).any():
        handles.append(Patch(facecolor="none", edgecolor=MUTED_TEXT, hatch="///",
                             label="Excluded by QC (all records shown)"))
    axis.legend(handles=handles, fontsize=7, frameon=False, loc="upper left",
                bbox_to_anchor=(1.0, 1.0))
    return to_png(figure)


def depth_regression_figure(per_sample, fit):
    """Paper Fig 2B analogue: mean depth versus distinct SNVs, with a 95% band."""

    data = per_sample[per_sample["Analysis_Included"].astype(bool) & per_sample["Mean_Depth"].notna()]
    colors = location_colors(data["Location_Name"])
    figure = new_figure(6.4, 4.0)
    axis = figure.subplots()
    style_axis(axis)
    if fit is not None and fit["lower"] is not None:
        axis.fill_between(fit["grid"], fit["lower"], fit["upper"], color="#cde2fb",
                          label="95% bootstrap band of the fitted line")
        axis.plot(fit["grid"], fit["slope"] * fit["grid"] + fit["intercept"],
                  color=CATEGORICAL[0], linewidth=2, label="Least-squares line")
    for name, group in data.groupby("Location_Name"):
        axis.scatter(group["Mean_Depth"], group["Distinct_SNVs_After_QC"], s=64,
                     color=colors[str(name)], edgecolor=SURFACE, linewidth=1.5,
                     label=str(name), zorder=3)
    axis.set_xlabel("Mean read depth across the genome (reads per position)")
    axis.set_ylabel("Distinct SNVs (after QC)")
    title = "Sequencing depth versus detected SNVs"
    if fit is not None:
        note = (f"n = {fit['n']}, slope = {fit['slope']:.3f} SNVs per unit of mean depth, "
                f"R² = {fit['r_squared']:.2f}, permutation p = {fit['permutation_p']:.3f}")
        if fit["small_sample"]:
            note += "\nFewer than 10 samples: descriptive only, not reliable inference"
        axis.text(0.01, 0.99, note, transform=axis.transAxes, va="top", fontsize=7,
                  color=MUTED_TEXT)
    axis.set_title(title, fontsize=10)
    axis.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    return to_png(figure)


BASE_COLORS = {"A": CATEGORICAL[0], "C": CATEGORICAL[1], "G": CATEGORICAL[2], "T": CATEGORICAL[3]}


def genome_position_figure(position_summary, genes):
    """Paper Fig 3 analogue: samples carrying each variant along the genome.

    Color shows the alternative base; marker shape shows the evidence source
    (filled circle = pooled samples reporting, open diamond = individual
    sequences carrying). A dark ring marks positions where pooled samples
    report more than one alternative base (red points in the paper).
    """

    figure = new_figure(9.5, 4.6)
    axis, gene_axis = figure.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [6, 1]})
    style_axis(axis)
    pooled = position_summary[position_summary["Pooled_Samples"] > 0]
    individual = position_summary[position_summary["Individual_Carriers"] > 0]
    crowded = len(position_summary) > 300  # thousands of SNVs need smaller marks
    size, ring, star = (16, 55, 8) if crowded else (56, 150, 26)
    axis.vlines(pooled["Position"], 0, pooled["Pooled_Samples"], color=GRID,
                linewidth=0.6 if crowded else 1.2)
    for base, color in BASE_COLORS.items():
        rows = pooled[pooled["Alt"].eq(base)]
        axis.scatter(rows["Position"], rows["Pooled_Samples"], s=size, color=color,
                     edgecolor=SURFACE, linewidth=0.5 if crowded else 1.2, zorder=3)
        carriers = individual[individual["Alt"].eq(base)]
        axis.scatter(carriers["Position"], carriers["Individual_Carriers"], s=70, marker="D",
                     facecolor="none", edgecolor=color, linewidth=1.8, zorder=4)
    multiple = pooled[pooled["Multiple_Alt_At_Position"]]
    if not multiple.empty:
        axis.scatter(multiple["Position"], multiple["Pooled_Samples"], s=ring, facecolor="none",
                     edgecolor=TEXT, linewidth=0.7 if crowded else 1.2, zorder=5)
    flagged = position_summary[position_summary["Published_Novelty"].astype(str).str.casefold().isin(
        ["yes", "true", "novel"])] if "Published_Novelty" in position_summary else pd.DataFrame()
    if not flagged.empty:
        axis.scatter(flagged["Position"], flagged["Pooled_Samples"], s=star, marker="*",
                     color=TEXT, linewidth=0, zorder=6)
    axis.set_ylabel("Samples with the variant")
    axis.set_title("Where variants occur along the reference genome", fontsize=10)
    axis.yaxis.get_major_locator().set_params(integer=True)
    axis.set_ylim(bottom=0)
    handles = [Patch(color=color, label=f"Alternative base {base}") for base, color in BASE_COLORS.items()]
    handles += [
        Line2D([], [], marker="o", linestyle="", color=MUTED_TEXT, markersize=7,
               label="Pooled samples reporting (filled circle)"),
        Line2D([], [], marker="D", linestyle="", markerfacecolor="none", color=MUTED_TEXT,
               markersize=7, label="Individual sequences carrying (open diamond)"),
    ]
    if not multiple.empty:
        handles.append(Line2D([], [], marker="o", linestyle="", markerfacecolor="none", color=TEXT,
                              markersize=10, label="Position with more than one pooled alternative base"))
    if not flagged.empty:
        handles.append(Line2D([], [], marker="*", linestyle="", color=TEXT, markersize=7,
                              label="Published novelty flag (source annotation)"))
    axis.legend(handles=handles, fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))

    gene_axis.set_facecolor(SURFACE)
    gene_axis.set_yticks([])
    for side in ["top", "right", "left"]:
        gene_axis.spines[side].set_visible(False)
    gene_axis.tick_params(colors=MUTED_TEXT, labelsize=8)
    span = (max((end for _, _, end in genes), default=1)
            - min((start for _, start, _ in genes), default=0))
    for index, (name, start, end) in enumerate(genes):
        gene_axis.add_patch(Rectangle((start, 0.2), end - start, 0.6,
                                      color=GRID if index % 2 == 0 else "#d3d2cd"))
        if end - start >= 0.04 * span:  # narrow bands stay unlabeled rather than overlap
            gene_axis.text((start + end) / 2, 0.5, name, ha="center", va="center",
                           fontsize=7, color=TEXT, clip_on=True)
    gene_axis.set_ylim(0, 1)
    gene_axis.set_xlabel("Reference position (1-based nucleotide coordinate)"
                         + ("; narrow annotation bands are unlabeled" if genes else ""))
    if not genes:
        gene_axis.text(0.5, 0.5, "No gene annotation supplied", transform=gene_axis.transAxes,
                       ha="center", va="center", fontsize=7, color=MUTED_TEXT)
    return to_png(figure)


def trajectory_figure(matrix, variant_keys):
    """Reported frequencies over time, with unreported observations kept apart."""

    figure = new_figure(9.0, 4.6)
    axis = figure.subplots()
    style_axis(axis)
    locations = sorted(matrix["Location_Name"].astype(str).unique())
    styles = ["-", "--", ":", "-."]
    for index, key in enumerate(variant_keys):
        color = CATEGORICAL[index % len(CATEGORICAL)]
        rows = matrix[matrix["Variant_Key"].eq(key)]
        for location_index, location in enumerate(locations):
            site = rows[rows["Location_Name"].astype(str).eq(location)].sort_values("Collection_Date")
            reported = site[site["Status"].eq(VARIANT_OBSERVED)]
            if not reported.empty:
                axis.plot(reported["Collection_Date"], reported["Alt_Frequency"], color=color,
                          linestyle=styles[location_index % len(styles)], linewidth=2,
                          marker=MARKERS[location_index % len(MARKERS)], markersize=8,
                          markeredgecolor=SURFACE)
            not_reported = site[site["Status"].astype(str).str.startswith("Adequate")]
            axis.scatter(not_reported["Collection_Date"], [-0.08] * len(not_reported), s=40,
                         facecolor="none", edgecolor=color, linewidth=1.5)
            uncertain = site[~site["Status"].eq(VARIANT_OBSERVED)
                             & ~site["Status"].astype(str).str.startswith("Adequate")]
            axis.scatter(uncertain["Collection_Date"], [-0.17] * len(uncertain), s=40,
                         marker="x", color=color, linewidth=1.5)
    axis.axhspan(-0.22, -0.02, color="#f1f0ec", zorder=0)
    axis.axhline(0, color=MUTED_TEXT, linewidth=0.8)
    axis.set_ylim(-0.22, 1.02)
    axis.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axis.text(1.0, -0.08, " not reported, adequate depth", transform=axis.get_yaxis_transform(),
              va="center", fontsize=7, color=MUTED_TEXT)
    axis.text(1.0, -0.17, " low/zero/unknown coverage or QC-excluded", transform=axis.get_yaxis_transform(),
              va="center", fontsize=7, color=MUTED_TEXT)
    date_axis(axis)
    axis.set_ylabel("Reported alternative-allele frequency")
    axis.set_title("Reported variant frequencies over time (gray strip: no frequency measured)",
                   fontsize=10)
    handles = [Line2D([], [], color=CATEGORICAL[index % len(CATEGORICAL)], linewidth=2, label=key)
               for index, key in enumerate(variant_keys)]
    handles += [Line2D([], [], color=MUTED_TEXT, linestyle=styles[index % len(styles)],
                       marker=MARKERS[index % len(MARKERS)], label=location)
                for index, location in enumerate(locations)]
    axis.legend(handles=handles, fontsize=7, frameon=False, loc="upper left",
                bbox_to_anchor=(1.12, 1.0))
    return to_png(figure)


def heatmap_figure(matrix, variant_keys):
    """Sample-by-variant grid: reported frequency, or why none was measured."""

    samples = matrix.drop_duplicates("Sample_ID").sort_values(
        ["Location_Name", "Collection_Date", "Sample_ID"])
    sample_ids = list(samples["Sample_ID"])
    figure = new_figure(max(6.5, 0.45 * len(variant_keys) + 3.5), max(3.2, 0.32 * len(sample_ids) + 1.8))
    axis = figure.subplots()
    axis.set_facecolor(SURFACE)
    # Only the displayed variants: the full matrix can have hundreds of thousands of rows.
    shown = matrix[matrix["Variant_Key"].isin(variant_keys)]
    lookup = {(row["Sample_ID"], row["Variant_Key"]): row for row in shown.to_dict("records")}
    for row_index, sample_id in enumerate(sample_ids):
        for column_index, key in enumerate(variant_keys):
            cell = lookup.get((sample_id, key))
            status = cell["Status"] if cell else "Coverage unknown"
            if status == VARIANT_OBSERVED:
                color, hatch = BLUE_RAMP(float(cell["Alt_Frequency"])), None
            elif str(status).startswith("Adequate"):
                color, hatch = NOT_REPORTED_FILL, None
            elif status == "Coverage unknown":
                color, hatch = UNKNOWN_FILL, "xx"
            elif status == EXCLUDED_BY_QC:
                color, hatch = UNKNOWN_FILL, ".."
            else:
                color, hatch = LOW_COVERAGE_FILL, "///"
            axis.add_patch(Rectangle((column_index, row_index), 1, 1, facecolor=color,
                                     edgecolor=SURFACE, linewidth=2, hatch=hatch))
    axis.set_xlim(0, len(variant_keys))
    axis.set_ylim(len(sample_ids), 0)
    axis.set_xticks([index + 0.5 for index in range(len(variant_keys))])
    axis.set_xticklabels([key.split(":", 1)[1] for key in variant_keys], rotation=60,
                         ha="right", fontsize=7)
    labels = [f"{row['Sample_ID']} ({row['Location_Name']}, "
              f"{row['Collection_Date'].date() if pd.notna(row['Collection_Date']) else 'no date'})"
              for row in samples.to_dict("records")]
    axis.set_yticks([index + 0.5 for index in range(len(sample_ids))])
    axis.set_yticklabels(labels, fontsize=7)
    for side in ["top", "right", "left", "bottom"]:
        axis.spines[side].set_visible(False)
    axis.tick_params(colors=MUTED_TEXT, length=0)
    axis.set_title("Pooled sample-by-variant frequencies (position:Ref>Alt)", fontsize=10, color=TEXT)
    scale = figure.colorbar(ScalarMappable(cmap=BLUE_RAMP), ax=axis, shrink=0.6)
    scale.set_label("Reported frequency", fontsize=8)
    scale.ax.tick_params(labelsize=7)
    handles = [
        Patch(facecolor=NOT_REPORTED_FILL, edgecolor=MUTED_TEXT, label="Not reported, adequate depth"),
        Patch(facecolor=LOW_COVERAGE_FILL, hatch="///", edgecolor=MUTED_TEXT, label="Low or zero coverage"),
        Patch(facecolor=UNKNOWN_FILL, hatch="xx", edgecolor=MUTED_TEXT, label="Coverage unknown"),
        Patch(facecolor=UNKNOWN_FILL, hatch="..", edgecolor=MUTED_TEXT, label="Reported but excluded by QC"),
    ]
    figure.legend(handles=handles, fontsize=7, frameon=False, loc="outside lower center", ncol=2)
    return to_png(figure)


def pcoa_figure(coordinates, pcoa, method, exploratory):
    """PCoA scatter colored by collection date, one marker shape per site, and
    an eigenvalue panel showing how much of the variation each axis explains."""

    figure = new_figure(9.6, 4.8)
    axis, scree = figure.subplots(1, 2, gridspec_kw={"width_ratios": [3.2, 1]})
    style_axis(axis)
    dates = coordinates["Collection_Date"]
    has_dates = dates.notna().any()
    numbers = date2num(dates.dropna()) if has_dates else None
    groups = coordinates["Location_Name"].astype(str)
    group_title = "Sampling site (marker shape); open = no date"
    if groups.nunique() > len(MARKERS) and groups.str.contains(":").all():
        # Too many sites for distinct shapes: group by the area before ":" (the
        # state in the paper's "STATE: LOCATION" labels), as the paper's Fig 5 did.
        groups = groups.str.split(":").str[0].str.strip()
        group_title = "Sampling area, the label before ':' (marker shape); open = no date"
    locations = sorted(groups.unique())
    scatter = None
    for index, location in enumerate(locations):
        rows = coordinates[groups.eq(location)]
        dated = rows[rows["Collection_Date"].notna()] if has_dates else rows
        undated = rows[rows["Collection_Date"].isna()] if has_dates else rows.iloc[0:0]
        marker = MARKERS[index % len(MARKERS)]
        colors = date2num(dated["Collection_Date"]) if has_dates else CATEGORICAL[0]
        scatter = axis.scatter(dated["PCoA_1"], dated["PCoA_2"], c=colors, cmap=BLUE_RAMP,
                               vmin=numbers.min() if has_dates else None,
                               vmax=numbers.max() if has_dates else None,
                               marker=marker, s=80, edgecolor=TEXT, linewidth=0.6, zorder=3)
        if not undated.empty:
            axis.scatter(undated["PCoA_1"], undated["PCoA_2"], marker=marker, s=80,
                         facecolor="none", edgecolor=MUTED_TEXT, linewidth=1.2, zorder=3)
    if len(coordinates) <= 30:
        for row in coordinates.to_dict("records"):
            axis.annotate(row["Sample_ID"], (row["PCoA_1"], row["PCoA_2"]), fontsize=6,
                          color=MUTED_TEXT, xytext=(4, 4), textcoords="offset points")
    explained = pcoa["explained"]
    first = explained[0] * 100 if len(explained) > 0 else 0
    second = explained[1] * 100 if len(explained) > 1 else 0
    axis.set_xlabel(f"PCoA axis 1 ({first:.1f}% of positive eigenvalues)")
    axis.set_ylabel(f"PCoA axis 2 ({second:.1f}% of positive eigenvalues)")
    suffix = "\n(EXPLORATORY - not the paper's Yue-Clayton PCoA)" if exploratory else ""
    axis.set_title(method + suffix, fontsize=9)
    if has_dates and scatter is not None:
        bar = figure.colorbar(scatter, ax=axis, shrink=0.7, pad=0.01)
        ticks = np.linspace(numbers.min(), numbers.max(), 4)
        bar.set_ticks(ticks)
        bar.set_ticklabels([num2date(value).strftime("%Y-%m-%d") for value in ticks])
        bar.set_label("Collection date", fontsize=8)
        bar.ax.tick_params(labelsize=7)

    style_axis(scree)
    eigenvalues = np.asarray(pcoa.get("eigenvalues", []), dtype=float)
    positive_sum = eigenvalues[eigenvalues > 0].sum()
    shown = eigenvalues[:8]
    if positive_sum > 0 and len(shown):
        percents = shown / positive_sum * 100
        colors = [CATEGORICAL[0] if value >= 0 else CATEGORICAL[1] for value in percents]
        scree.bar(np.arange(1, len(shown) + 1), percents, color=colors, width=0.7)
        scree.axhline(0, color=MUTED_TEXT, linewidth=0.8)
        scree.set_xticks(np.arange(1, len(shown) + 1))
    scree.set_xlabel("Axis")
    scree.set_ylabel("% of positive eigenvalues")
    negatives = pcoa.get("negative_eigenvalues", 0)
    plural = "s" if negatives != 1 else ""
    scree.set_title(f"Explained variation\n({negatives} negative eigenvalue{plural})", fontsize=8)
    handles = [Line2D([], [], marker=MARKERS[index % len(MARKERS)], linestyle="", color=MUTED_TEXT,
                      markersize=8, label=location) for index, location in enumerate(locations[:12])]
    if len(locations) > 12:
        handles.append(Line2D([], [], linestyle="", label=f"+ {len(locations) - 12} more sites"))
    figure.legend(handles=handles, fontsize=7, frameon=False, loc="outside lower center",
                  ncol=min(6, len(handles)), title=group_title, title_fontsize=7)
    return to_png(figure)


CATEGORY_COLORS = [CATEGORICAL[0], CATEGORICAL[2], CATEGORICAL[1], "#b0afaa", "#9a9994", "#5f5e5a"]


def cross_source_figure(counts):
    """How many variants fall in each cross-source evidence category."""

    figure = new_figure(8.0, 3.2)
    axis = figure.subplots()
    style_axis(axis)
    labels = list(counts["Category"])
    values = list(counts["Variants"])
    positions = range(len(labels))[::-1]
    axis.barh(list(positions), values, color=CATEGORY_COLORS[: len(labels)], height=0.6)
    for position, value in zip(positions, values):
        axis.text(value, position, f" {value}", va="center", fontsize=8, color=TEXT)
    axis.set_yticks(list(positions))
    axis.set_yticklabels(labels, fontsize=7)
    axis.set_xlabel("Distinct reference-based variants")
    axis.xaxis.get_major_locator().set_params(integer=True)
    axis.set_title("Cross-source evidence for each variant", fontsize=10)
    return to_png(figure)


def shared_timeline_figure(comparison, individual_variants, pooled_variants):
    """When each shared or source-specific variant was seen in each source."""

    keys = list(comparison.loc[comparison["Category"].str.startswith(("Shared", "Individual only",
                                                                      "Pooled only")), "Variant_Key"])
    figure = new_figure(8.5, max(2.8, 0.45 * len(keys) + 1.6))
    axis = figure.subplots()
    style_axis(axis)
    for row, key in enumerate(keys):
        people = individual_variants[individual_variants["Variant_Key"].eq(key)]
        pools = pooled_variants[pooled_variants["Variant_Key"].eq(key)]
        axis.scatter(people["Collection_Date"], [row + 0.15] * len(people), marker="^", s=70,
                     color=CATEGORICAL[2], edgecolor=SURFACE, linewidth=1.2, zorder=3)
        if not pools.empty:
            axis.scatter(pools["Collection_Date"], [row - 0.15] * len(pools), s=40 + 160 * pools["Alt_Frequency"],
                         color=CATEGORICAL[0], edgecolor=SURFACE, linewidth=1.2, zorder=3)
    axis.set_yticks(range(len(keys)))
    axis.set_yticklabels([key.split(":", 1)[1] for key in keys], fontsize=7)
    axis.set_ylim(-0.7, max(len(keys) - 0.3, 0.7))
    date_axis(axis)
    axis.set_title("When each variant was observed in each source (dates are context, not causation)",
                   fontsize=9)
    handles = [
        Line2D([], [], marker="^", linestyle="", color=CATEGORICAL[2], markersize=8,
               label="Individual sequence carrying"),
        Line2D([], [], marker="o", linestyle="", color=CATEGORICAL[0], markersize=8,
               label="Pooled sample reporting (size = frequency)"),
    ]
    axis.legend(handles=handles, fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    return to_png(figure)


def quality_figure(pooled_samples, individual_samples, coverage_threshold):
    """Genome coverage of pooled samples and callable span of individual sequences."""

    panels = [frame for frame in (pooled_samples, individual_samples) if frame is not None and not frame.empty]
    figure = new_figure(9.0, 3.6)
    axes = figure.subplots(1, max(1, len(panels)), squeeze=False)[0]
    index = 0
    if pooled_samples is not None and not pooled_samples.empty:
        axis = axes[index]
        index += 1
        style_axis(axis)
        data = pooled_samples.sort_values("Sample_ID").reset_index(drop=True)
        for position, row in data.iterrows():
            value = row["Genome_Coverage_Percent"]
            if pd.isna(value):
                axis.text(position, 2, "unknown", rotation=90, ha="center", va="bottom",
                          fontsize=6, color=MUTED_TEXT)
                continue
            included = bool(row.get("Analysis_Included", True))
            axis.bar(position, value, width=0.7, color=CATEGORICAL[0], alpha=1.0 if included else 0.35,
                     hatch=None if included else "///", edgecolor=SURFACE if included else MUTED_TEXT)
        axis.axhline(coverage_threshold, color=TEXT, linestyle="--", linewidth=1)
        axis.text(-0.4, 102, f"dashed line: >{coverage_threshold:.0f}% paper threshold; hatched = excluded",
                  ha="left", va="bottom", fontsize=7, color=TEXT)
        axis.set_xticks(range(len(data)))
        axis.set_xticklabels(data["Sample_ID"], rotation=60, ha="right", fontsize=7)
        axis.set_ylim(0, 112)
        axis.set_ylabel("Genome positions at >= 10x (%)")
        axis.set_title("Pooled samples: genome coverage", fontsize=9)
    if individual_samples is not None and not individual_samples.empty:
        axis = axes[index]
        style_axis(axis)
        data = individual_samples.sort_values("Sample_ID")
        mapped = data["Variant_Profile_Status"].eq("Mapped to reference")
        axis.bar(np.arange(len(data)), data["Genome_Coverage_Percent"],
                 color=[CATEGORICAL[2] if value else "#c9c8c3" for value in mapped], width=0.7)
        axis.set_xticks(range(len(data)))
        axis.set_xticklabels(data["Sample_ID"], rotation=60, ha="right", fontsize=7)
        axis.set_ylim(0, 105)
        axis.set_ylabel("Callable reference positions (%)")
        axis.set_title("Individual sequences: callable reference positions", fontsize=9)
        axis.legend(handles=[Patch(color=CATEGORICAL[2], label="Mapped to reference"),
                             Patch(color="#c9c8c3", label="Mapping not established")],
                    fontsize=7, frameon=False, loc="lower left")
    return to_png(figure)


# ----------------------------------------------------------------------------
# Figures added for the guided dashboard. Source colors are fixed everywhere:
# individual sequences green, pooled samples blue.
# ----------------------------------------------------------------------------

SOURCE_COLORS = {"Individual sequences": CATEGORICAL[2], "Pooled samples": CATEGORICAL[0]}
CATEGORY_LANES = [
    "Shared: observed in both sources",
    "Individual only: not detected in adequately covered pooled samples",
    "Pooled only: reference base called in individual sequences",
    "Unresolved: pooled coverage missing or insufficient",
    "Unresolved: individual sequences not callable here",
    "Incomparable: reference not verified as the same",
]


def source_color(sample_type):
    return SOURCE_COLORS["Individual sequences" if str(sample_type).startswith("Individual")
                         else "Pooled samples"]


def completeness_figure(rows):
    """Share of samples per source with each kind of metadata.

    A field that the uploaded files cannot provide is labeled with the reason
    (for example "not supplied") instead of being drawn as a 0% bar.
    """

    fields = [("With_Collection_Date", "Collection date"), ("With_Location", "Location or site"),
              ("With_Coordinates", "Coordinates"), ("With_Genome_Coverage", "Genome coverage")]
    figure = new_figure(8.0, 3.2)
    axis = figure.subplots()
    style_axis(axis)
    height = 0.8 / max(1, len(rows))
    for row_index, row in enumerate(rows):
        color = SOURCE_COLORS.get(row["Source"], OTHER_GRAY)
        for field_index, (key, _) in enumerate(fields):
            y = field_index + (row_index - (len(rows) - 1) / 2) * height
            value = row.get(key)
            if isinstance(value, (int, float)) and not isinstance(value, bool) and row["Samples"]:
                percent = 100 * value / row["Samples"]
                axis.barh(y, percent, height=height * 0.9, color=color)
                axis.text(percent + 1, y, f"{int(value)} of {row['Samples']}", va="center",
                          fontsize=7, color=TEXT)
            else:
                reason = value if isinstance(value, str) else "not supplied"
                axis.text(1, y, f"{row['Source']}: {reason}", va="center", fontsize=7,
                          color=MUTED_TEXT, style="italic")
    axis.set_yticks(range(len(fields)))
    axis.set_yticklabels([label for _, label in fields], fontsize=8)
    axis.invert_yaxis()
    axis.set_xlim(0, 118)
    axis.set_xlabel("Samples with the information (%)")
    axis.set_title("Metadata completeness by evidence source", fontsize=10)
    axis.legend(handles=[Patch(color=SOURCE_COLORS.get(row["Source"], OTHER_GRAY), label=row["Source"])
                         for row in rows], fontsize=7, frameon=False, loc="upper left",
                bbox_to_anchor=(1.0, 1.0))
    return to_png(figure)


def depth_distribution_figure(variants, minimum_depth, position_table=None):
    """Read depth of pooled variant records (and of all listed positions)."""

    panels = 2 if position_table is not None and not position_table.empty else 1
    figure = new_figure(9.0 if panels == 2 else 6.8, 3.4)
    axes = figure.subplots(1, panels, squeeze=False)[0]
    depths = pd.to_numeric(variants["Read_Depth"], errors="coerce")
    known = depths.notna() & depths.gt(0)
    axis = axes[0]
    style_axis(axis)
    if known.any():
        bins = np.logspace(0, np.log10(max(depths[known].max(), minimum_depth * 2)) + 0.05, 30)
        included = variants["Analysis_Included"].astype(bool)
        axis.hist([depths[known & included], depths[known & ~included]], bins=bins, stacked=True,
                  color=[CATEGORICAL[0], "#c9c8c3"], label=["Retained records", "Excluded by QC"],
                  edgecolor=SURFACE, linewidth=0.5)
        axis.set_xscale("log")
        axis.axvline(minimum_depth, color=TEXT, linestyle="--", linewidth=1)
        axis.text(minimum_depth, axis.get_ylim()[1] * 0.95, f" {minimum_depth}x paper minimum",
                  fontsize=7, color=TEXT, va="top")
        axis.legend(fontsize=7, frameon=False, loc="upper right")
    missing = int((~known).sum())
    axis.set_xlabel("Total read depth at the variant position (reads, log scale)")
    axis.set_ylabel("Variant records")
    axis.set_title("Depth of reported pooled variants"
                   + (f" ({missing} without depth)" if missing else ""), fontsize=9)
    if panels == 2:
        axis = axes[1]
        style_axis(axis)
        values = position_table["Depth"].to_numpy(dtype=float)
        zero = int((values == 0).sum())
        positive = values[values > 0]
        if len(positive):
            axis.hist(positive, bins=np.logspace(0, np.log10(positive.max()) + 0.05, 30),
                      color=CATEGORICAL[0], edgecolor=SURFACE, linewidth=0.5)
            axis.set_xscale("log")
        axis.axvline(10, color=TEXT, linestyle="--", linewidth=1)
        axis.set_xlabel("Read depth per listed position (reads, log scale)")
        axis.set_ylabel("Sample-positions")
        axis.set_title(f"All listed positions ({zero} with 0 reads; dashed = 10x)", fontsize=9)
    return to_png(figure)


QC_CHECKS = [("QC_Variant_Depth", "Variant depth >= 20x"), ("QC_Fisher_Test", "Fisher p < 0.05"),
             ("QC_Masked_Site", "Masked-site list"), ("QC_Sample_Coverage", "Sample coverage > 90%"),
             ("QC_Reference_Base", "Ref matches reference"), ("QC_Base_Quality", "Base quality Q20")]


def qc_outcome_figure(variants):
    """For each quality check: records passing, failing, or not verifiable."""

    figure = new_figure(8.6, 3.4)
    axis = figure.subplots()
    style_axis(axis)
    total = max(1, len(variants))
    labels = []
    for index, (column, label) in enumerate(QC_CHECKS):
        values = (variants[column].astype(str) if column in variants
                  else pd.Series(["Not evaluable"] * len(variants)))
        passed = int(values.eq("Pass").sum())
        failed = int(values.str.startswith("Fail").sum())
        other = len(variants) - passed - failed
        left = 0
        for count, color, hatch in [(passed, CATEGORICAL[0], None), (failed, CATEGORICAL[1], None),
                                    (other, "#e4e3df", "///")]:
            if count:
                axis.barh(index, 100 * count / total, left=left, color=color, hatch=hatch,
                          edgecolor=SURFACE if hatch is None else MUTED_TEXT, linewidth=0.6, height=0.6)
                left += 100 * count / total
        labels.append(f"{label}\n{passed} pass, {failed} fail, {other} not verifiable")
    axis.set_yticks(range(len(QC_CHECKS)))
    axis.set_yticklabels(labels, fontsize=7)
    axis.invert_yaxis()
    axis.set_xlim(0, 100)
    axis.set_xlabel(f"Pooled variant records (%, n = {len(variants)})")
    axis.set_title("Which quality checks could be applied to the pooled records", fontsize=10)
    axis.legend(handles=[Patch(color=CATEGORICAL[0], label="Pass"), Patch(color=CATEGORICAL[1], label="Fail"),
                         Patch(facecolor="#e4e3df", hatch="///", edgecolor=MUTED_TEXT,
                               label="Not evaluable or not verifiable")],
                fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    return to_png(figure)


CALL_GROUPS = [("R", "Reference base called", "#86b6ef", None),
               ("V", "Alternative base called", CATEGORICAL[1], None),
               ("N", "N or ambiguity code", "#c9c8c3", None),
               ("D", "Deletion relative to reference", CATEGORICAL[6], None),
               (".", "Outside the sequenced span", "#ecebe8", "xx")]


def callable_positions_figure(individual_samples, calls_by_sample, reference_length):
    """Per individual sequence: how each reference position was called."""

    data = individual_samples.sort_values("Sample_ID").reset_index(drop=True)
    figure = new_figure(max(6.5, 0.4 * len(data) + 3.5), 3.6)
    axis = figure.subplots()
    style_axis(axis)
    for index, row in data.iterrows():
        calls = calls_by_sample.get(row["Sample_ID"])
        if calls is None:
            axis.text(index, 2, "mapping not\nestablished", rotation=90, ha="center", va="bottom",
                      fontsize=6, color=MUTED_TEXT)
            continue
        bottom = 0
        for code, _, color, hatch in CALL_GROUPS:
            share = 100 * calls.count(code) / reference_length
            if share:
                axis.bar(index, share, bottom=bottom, color=color, hatch=hatch, width=0.7,
                         edgecolor=SURFACE if hatch is None else MUTED_TEXT, linewidth=0.5)
                bottom += share
    axis.set_xticks(range(len(data)))
    axis.set_xticklabels(data["Sample_ID"], rotation=60, ha="right", fontsize=7)
    axis.set_ylim(0, 105)
    axis.set_ylabel("Reference positions (%)")
    axis.set_title("Individual sequences: callable and uncallable reference positions", fontsize=9)
    axis.legend(handles=[Patch(facecolor=color, hatch=hatch, edgecolor=MUTED_TEXT, label=label)
                         for _, label, color, hatch in CALL_GROUPS],
                fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    return to_png(figure)


def sampling_timeline_figure(samples):
    """Collection dates by site for both sources, at their actual (irregular) dates.

    `samples` has Sample_Type, Location_Name, Collection_Date, Analysis_Included.
    """

    dated = samples[samples["Collection_Date"].notna()].copy()
    dated["Site"] = dated["Location_Name"].fillna("").astype(str).replace("", "(no location)")
    sites = sorted(dated["Site"].unique(), key=lambda name: (name == "(no location)", name))
    figure = new_figure(9.0, max(2.6, 0.3 * len(sites) + 1.6))
    axis = figure.subplots()
    style_axis(axis)
    rows = {site: index for index, site in enumerate(sites)}
    for sample_type, marker, offset in [("Individual sequence", "^", 0.15), ("Pooled sample", "o", -0.15)]:
        group = dated[dated["Sample_Type"].eq(sample_type)]
        if group.empty:
            continue
        color = source_color(sample_type)
        kept = group["Analysis_Included"].fillna(True).astype(bool)
        y = group["Site"].map(rows) + offset
        axis.scatter(group.loc[kept, "Collection_Date"], y[kept], marker=marker, s=60, color=color,
                     edgecolor=SURFACE, linewidth=1, zorder=3)
        axis.scatter(group.loc[~kept, "Collection_Date"], y[~kept], marker=marker, s=60,
                     facecolor="none", edgecolor=color, linewidth=1.4, zorder=3)
    axis.set_yticks(range(len(sites)))
    axis.set_yticklabels(sites, fontsize=7)
    axis.set_ylim(len(sites) - 0.5, -0.5)
    date_axis(axis)
    undated = len(samples) - len(dated)
    axis.set_title("Sampling dates and sites (points sit at the actual, irregular dates)"
                   + (f"; {undated} sample(s) without a date not shown" if undated else ""), fontsize=9)
    axis.legend(handles=[
        Line2D([], [], marker="^", linestyle="", color=SOURCE_COLORS["Individual sequences"], markersize=8,
               label="Individual sequence"),
        Line2D([], [], marker="o", linestyle="", color=SOURCE_COLORS["Pooled samples"], markersize=8,
               label="Pooled sample"),
        Line2D([], [], marker="o", linestyle="", markerfacecolor="none", color=MUTED_TEXT, markersize=8,
               label="Open = excluded by QC"),
    ], fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    return to_png(figure)


def frequency_change_figure(trends, limit=15):
    """First versus last reported frequency for variants reported twice or more at one site."""

    data = trends[trends["Change_In_Reported_Frequency"].notna()].copy()
    data["Size"] = data["Change_In_Reported_Frequency"].abs()
    data = data.sort_values(["Size", "Variant_Key"], ascending=[False, True]).head(limit)
    data = data.sort_values("Change_In_Reported_Frequency").reset_index(drop=True)
    figure = new_figure(9.0, max(2.6, 0.36 * len(data) + 1.4))
    axis = figure.subplots()
    style_axis(axis)
    for index, row in data.iterrows():
        rising = row["Change_In_Reported_Frequency"] > 0
        color = CATEGORICAL[1] if rising else CATEGORICAL[0]
        axis.plot([row["First_Reported_Frequency"], row["Last_Reported_Frequency"]], [index, index],
                  color=color, linewidth=2)
        axis.scatter(row["First_Reported_Frequency"], index, s=50, facecolor=SURFACE, edgecolor=color,
                     linewidth=1.6, zorder=3)
        axis.scatter(row["Last_Reported_Frequency"], index, s=50, color=color, edgecolor=SURFACE, zorder=3)
    labels = []
    for row in data.to_dict("records"):
        days = (row["Last_Reported_Date"] - row["First_Reported_Date"]).days
        labels.append(f"{row['Variant_Key'].split(':', 1)[1]} at {row['Location_Name'] or 'no site'}"
                      f" ({row['Reported']} reports over {days} days)")
    axis.set_yticks(range(len(data)))
    axis.set_yticklabels(labels, fontsize=7)
    axis.set_xlim(-0.02, 1.02)
    axis.set_xlabel("Reported alternative-allele frequency (open = first report, filled = last report)")
    axis.set_title("Change between the first and last reported frequency at the same site", fontsize=9)
    axis.legend(handles=[Patch(color=CATEGORICAL[1], label="Higher at the last report"),
                         Patch(color=CATEGORICAL[0], label="Lower at the last report")],
                fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    return to_png(figure)


STATUS_GROUPS = [
    ("Reported", lambda status: status.eq(VARIANT_OBSERVED), CATEGORICAL[0], None),
    ("Adequate depth, not reported", lambda status: status.str.startswith("Adequate"), "#ffffff", None),
    ("Low or no coverage", lambda status: status.str.startswith("Low") | status.str.startswith("Not covered"),
     LOW_COVERAGE_FILL, "///"),
    ("Coverage unknown", lambda status: status.eq("Coverage unknown"), UNKNOWN_FILL, "xx"),
    ("Reported but excluded by QC", lambda status: status.eq(EXCLUDED_BY_QC), UNKNOWN_FILL, ".."),
]


def observation_status_figure(matrix):
    """For each pooled sample: what is known at each analyzed variant position."""

    samples = matrix.drop_duplicates("Sample_ID").sort_values(["Collection_Date", "Sample_ID"])
    figure = new_figure(max(6.5, 0.3 * len(samples) + 3.5), 3.8)
    axis = figure.subplots()
    style_axis(axis)
    statuses = matrix.groupby("Sample_ID")["Status"].apply(lambda values: values.astype(str))
    for index, sample_id in enumerate(samples["Sample_ID"]):
        status = statuses.loc[sample_id]
        bottom = 0
        for _, test, color, hatch in STATUS_GROUPS:
            share = 100 * test(status).sum() / max(1, len(status))
            if share:
                axis.bar(index, share, bottom=bottom, color=color, hatch=hatch, width=0.75,
                         edgecolor=MUTED_TEXT, linewidth=0.4)
                bottom += share
    labels = [f"{row['Sample_ID']} ({row['Collection_Date'].date() if pd.notna(row['Collection_Date']) else 'no date'})"
              for row in samples.to_dict("records")]
    axis.set_xticks(range(len(samples)))
    axis.set_xticklabels(labels, rotation=70, ha="right", fontsize=6 if len(samples) > 25 else 7)
    axis.set_ylim(0, 100)
    axis.set_ylabel("Analyzed variants (%)")
    axis.set_title("What is known for each pooled sample at the analyzed variant positions", fontsize=9)
    axis.legend(handles=[Patch(facecolor=color, hatch=hatch, edgecolor=MUTED_TEXT, label=label)
                         for label, _, color, hatch in STATUS_GROUPS],
                fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    return to_png(figure)


PRECISION_ORDER = ["Source-provided coordinates", "Approximate (GeoNames place)",
                   "Broad administrative area", "Named place or site, no coordinates", "No location"]
PRECISION_STYLE = {
    "Source-provided coordinates": (CATEGORICAL[0], None),
    "Approximate (GeoNames place)": ("#86b6ef", None),
    "Broad administrative area": ("#c9c8c3", "///"),
    "Named place or site, no coordinates": ("#ecebe8", None),
    "No location": ("#ffffff", "xx"),
}


def geography_figure(samples):
    """Location precision by source, and a plain coordinate plot (not a map).

    `samples` has Sample_Type, Precision_Category, Latitude, Longitude.
    Broad administrative areas and sites without coordinates are counted but
    never drawn as points.
    """

    points = samples[samples["Precision_Category"].isin(PRECISION_ORDER[:2])
                     & samples["Latitude"].notna() & samples["Longitude"].notna()]
    panels = 2 if not points.empty else 1
    figure = new_figure(9.6 if panels == 2 else 7.0, 3.8)
    axes = figure.subplots(1, panels, squeeze=False)[0]
    axis = axes[0]
    style_axis(axis)
    sources = list(dict.fromkeys(samples["Sample_Type"]))
    for row, source in enumerate(sources):
        group = samples[samples["Sample_Type"].eq(source)]
        left = 0
        for category in PRECISION_ORDER:
            count = int(group["Precision_Category"].eq(category).sum())
            if count:
                color, hatch = PRECISION_STYLE[category]
                axis.barh(row, count, left=left, color=color, hatch=hatch, edgecolor=MUTED_TEXT,
                          linewidth=0.5, height=0.55)
                axis.text(left + count / 2, row, str(count), ha="center", va="center", fontsize=7, color=TEXT)
                left += count
    axis.set_yticks(range(len(sources)))
    axis.set_yticklabels(sources, fontsize=8)
    axis.invert_yaxis()
    axis.set_xlabel("Samples")
    axis.set_title("How precisely each sample's location is known", fontsize=9)
    axis.legend(handles=[Patch(facecolor=PRECISION_STYLE[name][0], hatch=PRECISION_STYLE[name][1],
                               edgecolor=MUTED_TEXT, label=name) for name in PRECISION_ORDER],
                fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
    if panels == 2:
        axis = axes[1]
        style_axis(axis)
        for source in sources:
            group = points[points["Sample_Type"].eq(source)]
            color = source_color(source)
            exact = group["Precision_Category"].eq(PRECISION_ORDER[0])
            marker = "^" if source.startswith("Individual") else "o"
            axis.scatter(group.loc[exact, "Longitude"], group.loc[exact, "Latitude"], marker=marker, s=60,
                         color=color, edgecolor=SURFACE, zorder=3)
            axis.scatter(group.loc[~exact, "Longitude"], group.loc[~exact, "Latitude"], marker=marker, s=60,
                         facecolor="none", edgecolor=color, linewidth=1.4, zorder=3)
        axis.set_xlabel("Longitude (degrees)")
        axis.set_ylabel("Latitude (degrees)")
        axis.set_title("Coordinates (filled = source-provided, open = GeoNames estimate)", fontsize=8)
        axis.margins(0.2)
    return to_png(figure)


EVIDENCE_STYLE = {
    "Reference base called": ("#ffffff", None),
    "Other alternative base called": (CATEGORICAL[3], None),
    "Adequate depth; variant not reported": ("#ffffff", None),
    "Low coverage": (LOW_COVERAGE_FILL, "///"),
    "Not covered or not callable": (LOW_COVERAGE_FILL, "///"),
    "Coverage unknown": (UNKNOWN_FILL, "xx"),
    EXCLUDED_BY_QC: (UNKNOWN_FILL, ".."),
}


def evidence_matrix_figure(cells, variant_keys, individual_ids, pooled_ids, pooled_labels):
    """Variant-by-sample evidence across both sources on one verified reference."""

    columns = list(individual_ids) + list(pooled_ids)
    figure = new_figure(max(7.0, 0.36 * len(columns) + 4.0), max(3.4, 0.3 * len(variant_keys) + 2.6))
    axis = figure.subplots()
    axis.set_facecolor(SURFACE)
    lookup = {(row["Variant_Key"], row["Sample_ID"], row["Sample_Type"]): row for row in cells}
    for row_index, key in enumerate(variant_keys):
        for column_index, sample_id in enumerate(columns):
            sample_type = "Individual sequence" if column_index < len(individual_ids) else "Pooled sample"
            cell = lookup.get((key, sample_id, sample_type))
            status = cell["Status"] if cell else "Coverage unknown"
            color, hatch = EVIDENCE_STYLE.get(status, (UNKNOWN_FILL, "xx"))
            if status == VARIANT_OBSERVED:
                hatch = None
                color = (SOURCE_COLORS["Individual sequences"] if sample_type.startswith("Individual")
                         else BLUE_RAMP(float(cell["Alt_Frequency"])))
            axis.add_patch(Rectangle((column_index, row_index), 1, 1, facecolor=color, hatch=hatch,
                                     edgecolor="#d3d2cd", linewidth=0.8))
    if individual_ids and pooled_ids:
        axis.axvline(len(individual_ids), color=TEXT, linewidth=2)
    axis.set_xlim(0, len(columns))
    axis.set_ylim(len(variant_keys), 0)
    axis.set_xticks([index + 0.5 for index in range(len(columns))])
    axis.set_xticklabels(list(individual_ids) + list(pooled_labels), rotation=70, ha="right", fontsize=6.5)
    axis.set_yticks([index + 0.5 for index in range(len(variant_keys))])
    axis.set_yticklabels([key.split(":", 1)[1] for key in variant_keys], fontsize=7)
    for side in ["top", "right", "left", "bottom"]:
        axis.spines[side].set_visible(False)
    axis.tick_params(colors=MUTED_TEXT, length=0)
    axis.set_title("Evidence for each variant (left of the line: individual sequences; right: pooled samples)",
                   fontsize=9, color=TEXT)
    scale = figure.colorbar(ScalarMappable(cmap=BLUE_RAMP), ax=axis, shrink=0.5, pad=0.01)
    scale.set_label("Pooled reported frequency", fontsize=7)
    scale.ax.tick_params(labelsize=6)
    handles = [
        Patch(facecolor=SOURCE_COLORS["Individual sequences"], label="Individual: variant base called"),
        Patch(facecolor="#ffffff", edgecolor=MUTED_TEXT,
              label="Reference base called / pooled: adequate depth, not reported"),
        Patch(facecolor=CATEGORICAL[3], label="Individual: other alternative base"),
        Patch(facecolor=LOW_COVERAGE_FILL, hatch="///", edgecolor=MUTED_TEXT, label="Not callable or low coverage"),
        Patch(facecolor=UNKNOWN_FILL, hatch="xx", edgecolor=MUTED_TEXT, label="Coverage unknown (not 0%)"),
        Patch(facecolor=UNKNOWN_FILL, hatch="..", edgecolor=MUTED_TEXT, label="Pooled call excluded by QC"),
    ]
    figure.legend(handles=handles, fontsize=6.5, frameon=False, loc="outside lower center", ncol=2)
    return to_png(figure)


def cross_source_positions_figure(comparison, reference_length=None):
    """Reference positions of each comparable variant, one lane per evidence category."""

    data = comparison[~comparison["Category"].str.startswith("Incomparable")]
    lanes = [category for category in CATEGORY_LANES if category in set(data["Category"])]
    figure = new_figure(9.0, max(2.4, 0.5 * len(lanes) + 1.4))
    axis = figure.subplots()
    style_axis(axis)
    for index, category in enumerate(lanes):
        rows = data[data["Category"].eq(category)]
        color = CATEGORY_COLORS[CATEGORY_LANES.index(category)]
        axis.scatter(rows["Position"], [index] * len(rows), marker="|", s=260, color=color, linewidth=2.2)
        axis.text(1.0, index, f"  {len(rows)}", transform=axis.get_yaxis_transform(), va="center",
                  fontsize=7, color=TEXT)
    axis.set_yticks(range(len(lanes)))
    axis.set_yticklabels([lane.split(":")[0] for lane in lanes], fontsize=8)
    axis.set_ylim(len(lanes) - 0.5, -0.5)
    if reference_length:
        axis.set_xlim(0, reference_length + 1)
    axis.set_xlabel("Reference position (1-based)")
    axis.set_title("Where shared and source-specific variants lie on the reference (count at right)",
                   fontsize=9)
    return to_png(figure)


def distance_figure(pairs, thresholds):
    """Sprint 2 pairs: p-distance distribution and genetic versus temporal distance."""

    max_p, max_days, max_km = thresholds
    figure = new_figure(9.6, 3.6)
    left, right = figure.subplots(1, 2)
    for axis in (left, right):
        style_axis(axis)
    distances = pd.to_numeric(pairs["P_Distance"], errors="coerce")
    known = distances.dropna()
    if not known.empty:
        upper = max(known.max(), max_p) * 1.05 or 0.01
        left.hist(known, bins=np.linspace(0, upper, 25), color=CATEGORICAL[2], edgecolor=SURFACE)
    left.axvline(max_p, color=TEXT, linestyle="--", linewidth=1)
    left.text(max_p, left.get_ylim()[1] * 0.95, f" cutoff {max_p:g}", fontsize=7, va="top")
    left.yaxis.get_major_locator().set_params(integer=True)
    left.set_xlabel("Pairwise p-distance (differing / comparable positions)")
    left.set_ylabel("Sequence pairs")
    left.set_title(f"Genetic distances ({len(known)} pairs)", fontsize=9)

    days = pd.to_numeric(pairs["Days_Apart"], errors="coerce")
    km = pd.to_numeric(pairs["Km_Apart"], errors="coerce")
    same = pairs["Same_Location"].astype(str).str.casefold().eq("true")
    within = km.le(max_km) | (km.isna() & same)
    beyond = km.gt(max_km)
    unresolved = ~within & ~beyond
    groups = [("Within the km cutoff (or same named place)", within, CATEGORICAL[0], "o"),
              ("Beyond the km cutoff", beyond, CATEGORICAL[1], "s"),
              ("Geography unresolved", unresolved, "#9a9994", "x")]
    right.add_patch(Rectangle((0, 0), max_days, max_p, facecolor="#cde2fb", alpha=0.45, edgecolor="none",
                              zorder=1))
    for label, mask, color, marker in groups:
        mask = mask & days.notna() & distances.notna()
        right.scatter(days[mask], distances[mask], s=36, color=color, marker=marker,
                      label=f"{label} ({int(mask.sum())})", zorder=3)
    right.set_xlabel("Days between collection dates")
    right.set_ylabel("p-distance")
    undated = int((days.isna() & distances.notna()).sum())
    right.set_title("Genetic vs temporal distance (shaded = both cutoffs met)"
                    + (f"; {undated} pairs lack dates" if undated else ""), fontsize=8)
    right.legend(fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    return to_png(figure)


def relationship_network_figure(sample_ids, strict_pairs, candidate_pairs):
    """Samples on a circle; solid lines = all three cutoffs met, dashed = geography unresolved."""

    names = list(dict.fromkeys(sample_ids))
    angles = np.linspace(0, 2 * np.pi, len(names), endpoint=False)
    position = {name: (np.cos(angle), np.sin(angle)) for name, angle in zip(names, angles)}
    figure = new_figure(6.4, 5.6)
    axis = figure.subplots()
    axis.set_facecolor(SURFACE)
    axis.set_aspect("equal")
    axis.axis("off")
    for pairs, style, color, width in [(candidate_pairs, "--", "#9a9994", 1.4),
                                       (strict_pairs, "-", CATEGORICAL[2], 2.6)]:
        for first, second in pairs:
            if first in position and second in position:
                (x1, y1), (x2, y2) = position[first], position[second]
                axis.plot([x1, x2], [y1, y2], linestyle=style, color=color, linewidth=width, zorder=1)
    connected = {name for pair in list(strict_pairs) + list(candidate_pairs) for name in pair}
    for name, (x, y) in position.items():
        axis.scatter(x, y, s=90, color=CATEGORICAL[2] if name in connected else "#d3d2cd",
                     edgecolor=SURFACE, linewidth=1.5, zorder=3)
        axis.text(x * 1.15, y * 1.15, name, ha="center", va="center", fontsize=7, color=TEXT)
    axis.set_xlim(-1.45, 1.45)
    axis.set_ylim(-1.45, 1.45)
    axis.set_title("Screening relationships among individual sequences (not transmission)", fontsize=9)
    figure.legend(handles=[Line2D([], [], color=CATEGORICAL[2], linewidth=2.6,
                                  label=f"Meets all three cutoffs ({len(strict_pairs)})"),
                           Line2D([], [], color="#9a9994", linestyle="--", linewidth=1.4,
                                  label=f"Genetic + date only; geography unresolved ({len(candidate_pairs)})")],
                  fontsize=7, frameon=False, loc="outside lower center")
    return to_png(figure)
