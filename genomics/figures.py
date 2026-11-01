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
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "<"]


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


def genome_position_figure(position_summary, genes):
    """Paper Fig 3 analogue: samples carrying each variant along the genome."""

    figure = new_figure(9.0, 4.2)
    axis, gene_axis = figure.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [6, 1]})
    style_axis(axis)
    pooled = position_summary[position_summary["Pooled_Samples"] > 0]
    individual = position_summary[position_summary["Individual_Carriers"] > 0]

    single = pooled[~pooled["Multiple_Alt_At_Position"]]
    multiple = pooled[pooled["Multiple_Alt_At_Position"]]
    axis.vlines(pooled["Position"], 0, pooled["Pooled_Samples"], color=GRID, linewidth=1.5)
    axis.scatter(single["Position"], single["Pooled_Samples"], s=64, color=CATEGORICAL[0],
                 edgecolor=SURFACE, linewidth=1.5, zorder=3, label="Pooled samples reporting")
    if not multiple.empty:
        axis.scatter(multiple["Position"], multiple["Pooled_Samples"], s=64, color=CATEGORICAL[1],
                     edgecolor=SURFACE, linewidth=1.5, zorder=3,
                     label="Pooled: position with more than one alternative base")
    if not individual.empty:
        axis.scatter(individual["Position"], individual["Individual_Carriers"], s=70, marker="D",
                     color=CATEGORICAL[2], edgecolor=SURFACE, linewidth=1.5, zorder=4,
                     label="Individual sequences carrying")
    flagged = position_summary[position_summary["Published_Novelty"].astype(str).str.casefold().isin(
        ["yes", "true", "novel"])] if "Published_Novelty" in position_summary else pd.DataFrame()
    if not flagged.empty:
        axis.scatter(flagged["Position"], flagged["Pooled_Samples"], s=150, facecolor="none",
                     edgecolor=TEXT, linewidth=1.2, zorder=5,
                     label="Published novelty flag (source annotation)")
    axis.set_ylabel("Samples with the variant")
    axis.set_title("Where variants occur along the reference genome", fontsize=10)
    axis.yaxis.get_major_locator().set_params(integer=True)
    axis.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))

    gene_axis.set_facecolor(SURFACE)
    gene_axis.set_yticks([])
    for side in ["top", "right", "left"]:
        gene_axis.spines[side].set_visible(False)
    gene_axis.tick_params(colors=MUTED_TEXT, labelsize=8)
    for index, (name, start, end) in enumerate(genes):
        gene_axis.add_patch(Rectangle((start, 0.2), end - start, 0.6,
                                      color=GRID if index % 2 == 0 else "#d3d2cd"))
        gene_axis.text((start + end) / 2, 0.5, name, ha="center", va="center",
                       fontsize=7, color=TEXT, clip_on=True)
    gene_axis.set_ylim(0, 1)
    gene_axis.set_xlabel("Reference position (1-based)")
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
    lookup = {(row["Sample_ID"], row["Variant_Key"]): row for row in matrix.to_dict("records")}
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
    """PCoA scatter colored by collection date with one marker shape per site."""

    figure = new_figure(7.2, 4.6)
    axis = figure.subplots()
    style_axis(axis)
    dates = coordinates["Collection_Date"]
    has_dates = dates.notna().any()
    numbers = date2num(dates.dropna()) if has_dates else None
    locations = sorted(coordinates["Location_Name"].astype(str).unique())
    scatter = None
    for index, location in enumerate(locations):
        rows = coordinates[coordinates["Location_Name"].astype(str).eq(location)]
        colors = date2num(rows["Collection_Date"]) if has_dates else CATEGORICAL[0]
        scatter = axis.scatter(rows["PCoA_1"], rows["PCoA_2"], c=colors, cmap=BLUE_RAMP,
                               vmin=numbers.min() if has_dates else None,
                               vmax=numbers.max() if has_dates else None,
                               marker=MARKERS[index % len(MARKERS)], s=80, edgecolor=TEXT,
                               linewidth=0.6, label=location, zorder=3)
    for row in coordinates.to_dict("records"):
        axis.annotate(row["Sample_ID"], (row["PCoA_1"], row["PCoA_2"]), fontsize=6,
                      color=MUTED_TEXT, xytext=(4, 4), textcoords="offset points")
    explained = pcoa["explained"]
    first = explained[0] * 100 if len(explained) > 0 else 0
    second = explained[1] * 100 if len(explained) > 1 else 0
    axis.set_xlabel(f"PCoA 1 ({first:.1f}% of positive eigenvalues)")
    axis.set_ylabel(f"PCoA 2 ({second:.1f}% of positive eigenvalues)")
    prefix = "EXPLORATORY - not the paper's PCoA: " if exploratory else ""
    axis.set_title(prefix + method, fontsize=9)
    if has_dates and scatter is not None:
        bar = figure.colorbar(scatter, ax=axis, shrink=0.7)
        ticks = np.linspace(numbers.min(), numbers.max(), 4)
        bar.set_ticks(ticks)
        bar.set_ticklabels([num2date(value).strftime("%Y-%m-%d") for value in ticks])
        bar.set_label("Collection date", fontsize=8)
        bar.ax.tick_params(labelsize=7)
    handles = [Line2D([], [], marker=MARKERS[index % len(MARKERS)], linestyle="", color=MUTED_TEXT,
                      markersize=8, label=location) for index, location in enumerate(locations)]
    figure.legend(handles=handles, fontsize=7, frameon=False, loc="outside lower center",
                  ncol=min(4, len(handles)), title="Sampling site (marker shape)", title_fontsize=7)
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
