"""Scientific figures for the integrated dashboard (PNG bytes, no Flask code).

Colors come from a colorblind-checked categorical palette, assigned in a fixed
order, and a single-hue blue ramp for frequencies. Every figure answers one
surveillance question; its underlying numbers are also offered as a CSV.
"""

from io import BytesIO

import matplotlib

matplotlib.use("Agg")

from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.collections import PatchCollection  # noqa: E402
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
    matrix = matrix[matrix["Collection_Date"].notna()]  # undated samples cannot be placed in time
    known = (matrix["Callable_Frequency"] if "Callable_Frequency" in matrix
             else matrix["Alt_Frequency"].where(matrix["Status"].eq(VARIANT_OBSERVED)))
    for index, key in enumerate(variant_keys):
        color = CATEGORICAL[index % len(CATEGORICAL)]
        rows = matrix[matrix["Variant_Key"].eq(key)]
        for location_index, location in enumerate(locations):
            site = rows[rows["Location_Name"].astype(str).eq(location)].sort_values("Collection_Date")
            values = known[site.index]
            measured = site[values.notna()]
            if not measured.empty:
                axis.plot(measured["Collection_Date"], values[measured.index], color=color,
                          linestyle=styles[location_index % len(styles)], linewidth=2,
                          marker=MARKERS[location_index % len(MARKERS)], markersize=8,
                          markeredgecolor=SURFACE)
            counted = measured[~measured["Status"].eq(VARIANT_OBSERVED)]
            axis.scatter(counted["Collection_Date"], values[counted.index], s=60, facecolor=SURFACE,
                         edgecolor=color, linewidth=1.6, zorder=3)
            not_reported = site[site["Status"].astype(str).str.startswith("Adequate") & values.isna()]
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
    axis.text(1.0, -0.08, " not reported, adequate depth, no base counts", transform=axis.get_yaxis_transform(),
              va="center", fontsize=7, color=MUTED_TEXT)
    axis.text(1.0, -0.17, " low/zero/unknown coverage or QC-excluded", transform=axis.get_yaxis_transform(),
              va="center", fontsize=7, color=MUTED_TEXT)
    date_axis(axis)
    axis.set_ylabel("Alternative-allele frequency")
    axis.set_title("Known variant frequencies over time (gray strip: frequency not known)",
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
        # Undated observations cannot be placed on a time axis (audit S3-04).
        people = individual_variants[individual_variants["Variant_Key"].eq(key)
                                     & individual_variants["Collection_Date"].notna()]
        pools = pooled_variants[pooled_variants["Variant_Key"].eq(key) & pooled_variants["Collection_Date"].notna()]
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
    """First versus last callable frequency for variants with two or more at one site.

    Callable values are reported frequencies plus base-count values at
    adequately covered unreported positions (0 when no alternative reads).
    """

    data = trends[trends["Change_In_Callable_Frequency"].notna()].copy()
    data["Size"] = data["Change_In_Callable_Frequency"].abs()
    data = data.sort_values(["Size", "Variant_Key"], ascending=[False, True]).head(limit)
    data = data.sort_values("Change_In_Callable_Frequency").reset_index(drop=True)
    figure = new_figure(9.0, max(2.6, 0.36 * len(data) + 1.4))
    axis = figure.subplots()
    style_axis(axis)
    for index, row in data.iterrows():
        rising = row["Change_In_Callable_Frequency"] > 0
        color = CATEGORICAL[1] if rising else CATEGORICAL[0]
        axis.plot([row["First_Callable_Frequency"], row["Last_Callable_Frequency"]], [index, index],
                  color=color, linewidth=2)
        axis.scatter(row["First_Callable_Frequency"], index, s=50, facecolor=SURFACE, edgecolor=color,
                     linewidth=1.6, zorder=3)
        axis.scatter(row["Last_Callable_Frequency"], index, s=50, color=color, edgecolor=SURFACE, zorder=3)
    labels = []
    for row in data.to_dict("records"):
        days = (row["Last_Callable_Date"] - row["First_Callable_Date"]).days
        missing = "" if row["Complete_Series"] else ", some samples unknown"
        labels.append(f"{row['Variant_Key'].split(':', 1)[1]} at {row['Location_Name'] or 'no site'}"
                      f" ({row['Callable_Dated']} known values over {days} days{missing})")
    axis.set_yticks(range(len(data)))
    axis.set_yticklabels(labels, fontsize=7)
    axis.set_xlim(-0.02, 1.02)
    axis.set_xlabel("Alternative-allele frequency (open = first known value, filled = last known value)")
    axis.set_title("Change between the first and last known frequency at the same site", fontsize=9)
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


# --- Sprint 5: depth-aware variant analysis --------------------------------------------

def vaf_time_series_figure(points, unknown, chosen):
    """Small multiples: one mutation at one site per panel, VAF over time with 95% intervals.

    Filled markers are detections, open markers on 0 are supported non-detections
    (0 alternative reads in base counts); gray crosses below the axis are samples
    whose value is unknown. Dashed bars come from reconstructed counts.
    """

    columns = 2 if len(chosen) > 1 else 1
    rows = int(np.ceil(len(chosen) / columns))
    figure = new_figure(9.0, 2.5 * rows + 0.6)
    axes = np.atleast_1d(figure.subplots(rows, columns, squeeze=False)).ravel()
    for axis, (variant, location), color in zip(axes, chosen, CATEGORICAL):
        style_axis(axis)
        data = points[points["Variant_Key"].eq(variant) & points["Location_Name"].eq(location)
                      & points["Collection_Date"].notna()]
        upper_values = [data["Wilson_Upper"].to_numpy(dtype=float), data["VAF"].to_numpy(dtype=float), [0.05]]
        if "BetaBinomial_Upper" in data:
            upper_values.append(data["BetaBinomial_Upper"].to_numpy(dtype=float))
        top = float(np.nanmax(np.concatenate(upper_values)))
        top = min(1.0, top * 1.15)
        axis.plot(data["Collection_Date"], data["VAF"], color=color, linewidth=1.5, zorder=2)
        for row in data.to_dict("records"):
            detected = row["State"] == "Detected"
            if pd.notna(row.get("BetaBinomial_Lower", np.nan)):
                axis.vlines(row["Collection_Date"], row["BetaBinomial_Lower"], row["BetaBinomial_Upper"],
                            color=color, linewidth=6, alpha=0.25, zorder=1)
            if pd.notna(row["Wilson_Lower"]):
                approximate = str(row["Interval_Status"]).startswith("Approximate")
                axis.vlines(row["Collection_Date"], row["Wilson_Lower"], row["Wilson_Upper"], color=color,
                            linewidth=1.6, linestyle="--" if approximate else "-", zorder=2)
            axis.scatter([row["Collection_Date"]], [row["VAF"]], s=46, zorder=3, linewidth=1.5,
                         color=color if detected else SURFACE, edgecolor=color)
        gaps = unknown[unknown["Variant_Key"].eq(variant) & unknown["Location_Name"].eq(location)]
        axis.scatter(gaps["Collection_Date"], [-0.08 * top] * len(gaps), marker="x", s=30,
                     color=OTHER_GRAY, linewidth=1.2, zorder=3)
        axis.axhline(0, color=MUTED_TEXT, linewidth=0.8)
        axis.set_ylim(-0.14 * top, top)
        date_axis(axis)
        axis.set_title(f"{variant.split(':', 1)[1]} at {location or 'no site'}", fontsize=8)
        axis.set_ylabel("VAF", fontsize=8)
    for axis in axes[len(chosen):]:
        axis.set_visible(False)
    figure.legend(handles=[
        Line2D([], [], marker="o", linestyle="", color=MUTED_TEXT, label="Detected"),
        Line2D([], [], marker="o", linestyle="", markerfacecolor=SURFACE, color=MUTED_TEXT,
               label="Supported non-detection or uncalled reads (base counts)"),
        Line2D([], [], color=MUTED_TEXT, label="95% Wilson interval (read sampling only)"),
        Line2D([], [], color=MUTED_TEXT, linestyle="--", label="Interval from reconstructed counts"),
        Line2D([], [], marker="x", linestyle="", color=OTHER_GRAY, label="Unknown: no value (not 0)"),
    ] + ([Line2D([], [], color=MUTED_TEXT, linewidth=6, alpha=0.25,
                 label="Calibrated beta-binomial 95% interval (technical replicates)")]
         if "BetaBinomial_Lower" in points and points["BetaBinomial_Lower"].notna().any() else []),
        fontsize=7, frameon=False, loc="outside lower center", ncol=3)
    return to_png(figure)


def vaf_depth_figure(evidence, minimum_depth):
    """Why depth matters: interval width falls with depth, for the same VAF.

    Left: each observation's 95% interval width against its depth, with the
    theoretical width for three VAFs. Right: VAF against depth, shaded by width.
    """

    data = evidence[evidence["Wilson_Lower"].notna() & evidence["Counted_Once"]].copy()
    data["Depth"] = data["Total_Depth"].astype(float)
    figure = new_figure(10.0, 3.8)
    left, right = figure.subplots(1, 2)
    for axis in (left, right):
        style_axis(axis)
        axis.set_xscale("log")
        axis.axvline(minimum_depth, color=MUTED_TEXT, linewidth=0.8, linestyle=":")
    left.scatter(data["Depth"], data["Interval_Width"], s=8, color=OTHER_GRAY, alpha=0.5, linewidth=0,
                 label="Observations")
    grid = np.geomspace(max(1.0, data["Depth"].min() if len(data) else 10),
                        max(20.0, data["Depth"].max() if len(data) else 1000), 120)
    for vaf, color in zip([0.02, 0.2, 0.5], CATEGORICAL):
        left.plot(grid, expected_width(vaf, grid), color=color, linewidth=1.8, label=f"Expected width, VAF = {vaf}")
    left.set_xlabel("Total depth (reads, log scale)")
    left.set_ylabel("95% interval width")
    left.set_title("Interval width shrinks as depth grows", fontsize=9)
    left.legend(fontsize=7, frameon=False)
    widths = data["Interval_Width"]
    bins = pd.cut(widths, [-np.inf, 0.02, 0.1, np.inf], labels=["width <= 0.02", "0.02-0.1", "> 0.1"])
    for label, color in zip(bins.cat.categories, ["#0d366b", "#3987e5", "#cde2fb"]):
        chosen = data[bins.eq(label)]
        right.scatter(chosen["Depth"], chosen["VAF"], s=10, color=color, edgecolor=MUTED_TEXT,
                      linewidth=0.2, label=str(label))
    right.set_yscale("symlog", linthresh=0.01)
    right.set_ylim(0, 1.05)
    right.set_xlabel("Total depth (reads, log scale)")
    right.set_ylabel("Observed VAF")
    right.set_title(f"Observed VAF by depth (dotted line: {minimum_depth}x QC minimum)", fontsize=9)
    right.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0),
                 title="95% interval", title_fontsize=7)
    return to_png(figure)


def expected_width(vaf, depth):
    from genomics.frequency_uncertainty import expected_interval_width  # local import keeps figures light
    return expected_interval_width(vaf, depth)


def depth_dependence_figure(table):
    """The same VAF (0.20) at different depths: deeper sequencing, narrower interval."""

    figure = new_figure(7.5, 0.5 * len(table) + 1.4)
    axis = figure.subplots()
    style_axis(axis)
    for index, row in enumerate(table.to_dict("records")):
        axis.hlines(index, row["Wilson_Lower"], row["Wilson_Upper"], color=CATEGORICAL[0], linewidth=3)
        axis.scatter(row["VAF"], index, s=50, color=TEXT, zorder=3)
        axis.text(row["Wilson_Upper"] + 0.01, index, f"{row['Wilson_Lower']:.3f} to {row['Wilson_Upper']:.3f}",
                  va="center", fontsize=7, color=MUTED_TEXT)
    axis.set_yticks(range(len(table)))
    axis.set_yticklabels([f"{row['Alt_Reads']} of {row['Depth']} reads" for row in table.to_dict("records")],
                         fontsize=8)
    axis.set_xlim(0, 0.55)
    axis.set_xlabel("Variant allele frequency (dot) and 95% Wilson interval (bar)")
    axis.set_title("Same VAF, different depth: the interval narrows with more reads", fontsize=9)
    return to_png(figure)


def uncertainty_models_figure(evidence, scenarios, calibrated_rhos, reference_rho=None):
    """Binomial versus beta-binomial interval width (Sprint 5 Phase 2).

    Left: theoretical width at a frequency of 0.20 against depth for the
    binomial (Wilson) interval and each beta-binomial rho shown. Right: the
    uploaded observations, Wilson width against the overdispersion-aware
    width (calibrated where available, otherwise the reference hypothetical
    rho). Points above the diagonal are wider than read sampling alone implies.
    """

    from genomics.beta_binomial import expected_score_width, score_interval
    from genomics.frequency_uncertainty import expected_interval_width

    figure = new_figure(10.0, 3.9)
    left, right = figure.subplots(1, 2)
    for axis in (left, right):
        style_axis(axis)
    depths = np.unique(np.round(np.logspace(1, 4, 120)).astype(int))
    left.plot(depths, expected_interval_width(0.2, depths), color=TEXT, linewidth=2,
              label="Binomial (Wilson): read sampling only")
    lines = [(rho, f"Calibrated rho {rho:.3g} (upper 95% bound)", "-") for rho in calibrated_rhos[:3]]
    lines += [(rho, f"Hypothetical rho {rho:g}" if "Hypothetical" in source else f"Borrowed rho {rho:.3g}", "--")
              for rho, source in scenarios[:4]
              if not any(abs(rho - shown) < 1e-12 for shown in calibrated_rhos)]  # the same curve once
    for (rho, label, style), color in zip(lines, CATEGORICAL):
        left.plot(depths, expected_score_width(0.2, depths, rho), color=color, linewidth=1.8, linestyle=style,
                  label=label)
    left.set_xscale("log")
    left.set_xlabel("Sequencing depth (reads, log scale)")
    left.set_ylabel("Width of the 95% interval at VAF 0.20")
    left.set_title("More reads narrow the binomial interval;\nlibrary variation sets a floor", fontsize=9)
    left.legend(fontsize=7, frameon=False)

    data = evidence[evidence["Counted_Once"] & evidence["Wilson_Lower"].notna()]
    wilson_width = (data["Wilson_Upper"] - data["Wilson_Lower"]).to_numpy(dtype=float)
    calibrated = data["BetaBinomial_Width"].notna().to_numpy() if "BetaBinomial_Width" in data else np.zeros(len(data), bool)
    depth = data["Total_Depth"].astype(float).to_numpy()
    if reference_rho is not None and (~calibrated).any():
        lower, upper = score_interval(data["Alt_Read_Count"].astype(float).to_numpy()[~calibrated], depth[~calibrated],
                                      reference_rho)
        right.scatter(wilson_width[~calibrated], upper - lower, s=8, color=OTHER_GRAY, alpha=0.5, linewidth=0,
                      rasterized=True, label=f"Not calibrated: hypothetical rho {reference_rho:g}")
    if calibrated.any():
        right.scatter(wilson_width[calibrated], data["BetaBinomial_Width"].to_numpy(dtype=float)[calibrated], s=16,
                      color=CATEGORICAL[0], linewidth=0, label="Calibrated from technical replicates")
    top = max(0.05, float(np.nanmax(np.r_[wilson_width, 0.05])))
    right.plot([0, top], [0, top], color=MUTED_TEXT, linewidth=0.8, linestyle=":", label="Equal width")
    right.set_xlabel("Binomial (Wilson) interval width")
    right.set_ylabel("Overdispersion-aware interval width")
    right.set_title("Uploaded observations", fontsize=9)
    right.legend(fontsize=7, frameon=False, loc="upper left")
    return to_png(figure)


def replicate_concordance_figure(sites, calibration):
    """Frequencies of the same allele in two technical replicates of one specimen."""

    used = sites[sites["Used"]]
    figure = new_figure(5.6, 4.6)
    axis = figure.subplots()
    style_axis(axis)
    axis.plot([0, 1], [0, 1], color=MUTED_TEXT, linewidth=0.8, linestyle=":")
    results = calibration.set_index("Specimen_ID")
    for (specimen, group), color, marker in zip(used.groupby("Specimen_ID", sort=True), CATEGORICAL, MARKERS):
        row = results.loc[specimen]
        if pd.notna(row["Rho_Estimate"]):
            label = (f"{specimen}: rho {row['Rho_Estimate']:.3g} (range {row['Rho_Lower_95']:.3g} to "
                     f"{row['Rho_Upper_95']:.3g}), {len(group)} sites")
        else:
            label = f"{specimen}: not estimated ({len(group)} sites)"
        axis.scatter(group["VAF_1"].astype(float), group["VAF_2"].astype(float), s=22, color=color, marker=marker,
                     linewidth=0, alpha=0.85, label=label)
    top = float(np.nanmax(np.r_[used["VAF_1"].astype(float), used["VAF_2"].astype(float), 0.1])) * 1.08
    axis.set_xlim(0, min(1, top))
    axis.set_ylim(0, min(1, top))
    axis.set_xlabel("VAF in replicate 1")
    axis.set_ylabel("VAF in replicate 2")
    axis.set_title("Technical replicates of the same specimen\n(dotted line: perfect agreement)", fontsize=9)
    axis.legend(fontsize=6.5, frameon=False, loc="upper left")
    return to_png(figure)


# ----------------------------------------------------------------------------
# Sprint 5 Phase 3: exploratory mutation-frequency groups
# ----------------------------------------------------------------------------

def mutation_group_heatmap_figure(F, row_labels, row_groups, unit_labels):
    """Known frequencies (blue) by mutation (rows, ordered by group) and unit (columns, by date).

    Unknown cells are light gray with a hatch: they are not zeros.
    """

    rows, columns = F.shape
    figure = new_figure(min(11.5, 3.2 + 0.32 * columns), min(13, 1.6 + 0.2 * rows))
    axis = figure.subplots()
    axis.set_facecolor(UNKNOWN_FILL)
    # pcolormesh, not imshow: imshow resamples into full-canvas float buffers (~65 MB for 60 x 52 on mmc2).
    image = axis.pcolormesh(np.arange(columns + 1) - 0.5, np.arange(rows + 1) - 0.5, np.ma.masked_invalid(F),
                            cmap=BLUE_RAMP, vmin=0, vmax=1)
    axis.set_xlim(-0.5, columns - 0.5)
    axis.set_ylim(rows - 0.5, -0.5)  # first row at the top, as in a table
    # One collection for every unknown cell (thousands of separate patches cost ~20 MB on mmc2).
    unknown = [Rectangle((j - 0.5, i - 0.5), 1, 1) for i, j in zip(*np.nonzero(~np.isfinite(F)))]
    if unknown:
        axis.add_collection(PatchCollection(unknown, facecolor=UNKNOWN_FILL, edgecolor=GRID, hatch="///",
                                            linewidth=0))
    boundaries = [i for i in range(1, rows) if row_groups[i] != row_groups[i - 1]]
    for boundary in boundaries:
        axis.axhline(boundary - 0.5, color=TEXT, linewidth=1.2)
    axis.set_yticks(range(rows))
    axis.set_yticklabels([f"{group}  {label}" if group else f"—  {label}"
                          for group, label in zip(row_groups, row_labels)], fontsize=6.5)
    axis.set_xticks(range(columns))
    axis.set_xticklabels(unit_labels, rotation=60, ha="right", fontsize=6.5)
    axis.tick_params(colors=MUTED_TEXT, length=0)
    for side in axis.spines.values():
        side.set_visible(False)
    bar = figure.colorbar(image, ax=axis, shrink=0.6, pad=0.01)
    bar.set_label("Observed frequency (VAF)", fontsize=8)
    bar.ax.tick_params(labelsize=7)
    axis.set_title("Mutation frequencies by exploratory group (hatched gray = unknown, not 0)", fontsize=9,
                   color=TEXT)
    return to_png(figure)


def group_trajectory_figure(trajectories, group_ids):
    """Mean member frequency of each group over time, with the member range as a bar; one line style per site."""

    figure = new_figure(7.6, 4.4)
    axis = figure.subplots()
    style_axis(axis)
    sites = sorted(trajectories["Location_Name"].astype(str).unique())
    styles = ["-", "--", ":", "-."]
    dated = trajectories[trajectories["Collection_Date"].notna()].copy()
    dated["Collection_Date"] = pd.to_datetime(dated["Collection_Date"])
    for number, group_id in enumerate(group_ids):
        color = CATEGORICAL[number % len(CATEGORICAL)]
        rows = dated[dated["Group_ID"].eq(group_id)]
        for site_number, site in enumerate(sites):
            part = rows[rows["Location_Name"].astype(str).eq(site)].sort_values("Collection_Date")
            if part.empty:
                continue
            axis.errorbar(part["Collection_Date"], part["Mean_Member_VAF"],
                          yerr=[part["Mean_Member_VAF"] - part["Min_Member_VAF"],
                                part["Max_Member_VAF"] - part["Mean_Member_VAF"]],
                          color=color, linestyle=styles[site_number % len(styles)], marker=MARKERS[number],
                          markersize=4, linewidth=1.3, capsize=2, elinewidth=0.8)
    handles = [Line2D([], [], color=CATEGORICAL[n % len(CATEGORICAL)], marker=MARKERS[n], linewidth=1.3,
                      label=f"Group {group_id}") for n, group_id in enumerate(group_ids)]
    if len(sites) > 1:
        handles += [Line2D([], [], color=MUTED_TEXT, linestyle=styles[n % len(styles)], label=site)
                    for n, site in enumerate(sites[:len(styles)])]
    axis.legend(handles=handles, fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    locator = AutoDateLocator()
    axis.xaxis.set_major_locator(locator)
    axis.xaxis.set_major_formatter(ConciseDateFormatter(locator))
    axis.set_ylim(0, 1)
    axis.set_ylabel("Mean frequency of the group's mutations")
    axis.set_title("Exploratory groups over time (bars: lowest to highest member frequency)", fontsize=9)
    return to_png(figure)


def grouping_comparison_figure(depth_labels, ordinary_labels):
    """How many mutations each pair of (depth-aware group, ordinary group) shares; '—' = ungrouped."""

    rows = sorted(set(depth_labels), key=lambda v: (v == "—", len(v), v))
    columns = sorted(set(ordinary_labels), key=lambda v: (v == "—", len(v), v))
    table = np.zeros((len(rows), len(columns)), dtype=int)
    for a, b in zip(depth_labels, ordinary_labels):
        table[rows.index(a), columns.index(b)] += 1
    figure = new_figure(min(10, max(6.0, 2.6 + 0.45 * len(columns))), min(9, max(3.0, 1.8 + 0.35 * len(rows))))
    axis = figure.subplots()
    image = axis.imshow(np.where(table > 0, table, np.nan), cmap=BLUE_RAMP, aspect="auto")
    for i in range(len(rows)):
        for j in range(len(columns)):
            if table[i, j]:
                axis.text(j, i, str(table[i, j]), ha="center", va="center", fontsize=7,
                          color="white" if table[i, j] > table.max() / 2 else TEXT)
    axis.set_xticks(range(len(columns)))
    axis.set_xticklabels(columns, fontsize=7, rotation=45, ha="right")
    axis.set_yticks(range(len(rows)))
    axis.set_yticklabels(rows, fontsize=7)
    axis.set_xlabel("Ordinary clustering group (depth ignored)")
    axis.set_ylabel("Depth-aware group")
    for side in axis.spines.values():
        side.set_visible(False)
    axis.tick_params(length=0, colors=MUTED_TEXT)
    del image
    axis.set_title("Mutations shared between the two methods' groups\n('—' = not grouped by that method)",
                   fontsize=9)
    return to_png(figure)


ANNOTATION_COLORS = {
    "Synonymous": CATEGORICAL[0], "Missense": CATEGORICAL[1], "Stop gained": CATEGORICAL[2],
    "Stop lost": CATEGORICAL[3], "Stop retained": CATEGORICAL[4], "Start lost": CATEGORICAL[5],
    "Differs between overlapping coding regions": CATEGORICAL[6], "Noncoding": OTHER_GRAY,
    "Unsupported or ambiguous": LOW_COVERAGE_FILL, "Not annotated": UNKNOWN_FILL,
}


def cds_lanes(cds_rows):
    """Lane number per CDS so that overlapping coding regions are drawn on separate lines."""

    lanes, ends = [], []
    for cds in sorted(cds_rows, key=lambda item: min(start for start, _ in item["Parts"])):
        start = min(part[0] for part in cds["Parts"])
        end = max(part[1] for part in cds["Parts"])
        lane = next((index for index, last in enumerate(ends) if last < start), None)
        if lane is None:
            lane = len(ends)
            ends.append(end)
        else:
            ends[lane] = end
        lanes.append((cds, lane))
    return lanes


def mutation_annotation_figure(consequences, summary, cds_rows, reference_length):
    """Top: SNV positions by predicted category over the CDS map. Bottom: SNVs per CDS by consequence."""

    present = [name for name in ANNOTATION_COLORS if name in set(summary["Consequence_Category"])]
    per_cds = consequences[consequences["CDS_ID"].ne("")]
    labels = []
    genes = [cds["Gene"] or cds["CDS_ID"] for cds in cds_rows]
    for cds in cds_rows:
        name = cds["Gene"] or cds["CDS_ID"]
        labels.append(f"{name} ({cds['Product'][:28]})" if genes.count(name) > 1 and cds["Product"] else name)
    bar_rows = len(cds_rows) + 1
    figure = new_figure(10.5, 3.0 + 0.24 * len(present) + 0.26 * bar_rows)
    grid = figure.add_gridspec(2, 1, height_ratios=[1.0 + 0.12 * len(present), 0.35 + 0.16 * bar_rows])
    top, bottom = figure.add_subplot(grid[0]), figure.add_subplot(grid[1])
    style_axis(top)
    style_axis(bottom)

    lanes = cds_lanes(cds_rows)
    lane_count = max([lane for _, lane in lanes], default=-1) + 1
    for cds, lane in lanes:
        y = -1.0 - lane * 0.8
        for start, end in cds["Parts"]:
            top.add_patch(Rectangle((start, y - 0.28), end - start + 1, 0.56, facecolor="#dfe9f7",
                                    edgecolor="#1c5cab", linewidth=0.6,
                                    hatch="//" if cds["Status"] != "Supported" else None))
        first = min(part[0] for part in cds["Parts"])
        last = max(part[1] for part in cds["Parts"])
        if (last - first) > 0.025 * reference_length:
            name = cds["Gene"] or cds["CDS_ID"]
            if genes.count(name) > 1 and cds["Product"]:
                name = cds["Product"]
            top.text((first + last) / 2, y, name + (" (- strand)" if cds["Strand"] == "-" else ""),
                     ha="center", va="center", fontsize=7, color=TEXT)
    for row, name in enumerate(present):
        positions = summary.loc[summary["Consequence_Category"].eq(name), "Position"].astype(float)
        top.scatter(positions, np.full(len(positions), row), marker="|", s=60, linewidths=1.0,
                    color=ANNOTATION_COLORS[name] if name not in {"Unsupported or ambiguous", "Not annotated"}
                    else MUTED_TEXT)
    top.set_yticks(list(range(len(present))) + [-1.0 - lane * 0.8 for lane in range(lane_count)])
    top.set_yticklabels([f"{name} ({int(summary['Consequence_Category'].eq(name).sum())})" for name in present]
                        + ["Coding regions" if lane == 0 else "" for lane in range(lane_count)], fontsize=7.5)
    top.set_ylim(-1.0 - (lane_count - 1) * 0.8 - 0.6, len(present) - 0.4)
    top.set_xlim(0, reference_length)
    top.set_xlabel("Reference position (1-based)")
    top.set_title("Distinct SNVs along the reference by predicted consequence", fontsize=10, loc="left")
    top.grid(False, axis="y")

    categories = [name for name in ANNOTATION_COLORS if name not in {"Noncoding", "Differs between overlapping "
                                                                     "coding regions", "Not annotated"}]
    rows = []
    for cds in cds_rows:
        values = per_cds[per_cds["CDS_ID"].eq(cds["CDS_ID"]) & per_cds["Reference_ID"].eq(cds["Reference_ID"])]
        counts = values["Consequence"].value_counts()
        rows.append([int(counts.get(name, 0)) for name in categories])
    noncoding = int(summary["Consequence_Category"].eq("Noncoding").sum())
    labels.append("Outside coding regions")
    rows.append([0] * len(categories))
    y = np.arange(len(labels))[::-1]
    left = np.zeros(len(labels))
    for index, name in enumerate(categories):
        values = np.array([row[index] for row in rows], dtype=float)
        if values.any():
            bottom.barh(y, values, left=left, color=ANNOTATION_COLORS[name], edgecolor=SURFACE, linewidth=1.0,
                        height=0.7, label=name, hatch="//" if name == "Unsupported or ambiguous" else None)
            left += values
    if noncoding:
        bottom.barh(y[-1], noncoding, color=ANNOTATION_COLORS["Noncoding"], edgecolor=SURFACE, height=0.7,
                    label="Noncoding")
        left[-1] = noncoding
    for value, position in zip(left, y):
        bottom.text(value, position, f" {int(value)}", va="center", fontsize=7, color=MUTED_TEXT)
    bottom.set_yticks(y)
    bottom.set_yticklabels(labels, fontsize=7.5)
    bottom.set_xlim(0, max(left.max(), 1) * 1.12)
    bottom.set_xlabel("Distinct SNVs (a base in two overlapping coding regions is counted in each)")
    bottom.set_title("Predicted consequences per coding region", fontsize=10, loc="left")
    bottom.legend(fontsize=7, frameon=False, loc="lower right", ncol=2)
    bottom.grid(False, axis="y")
    return to_png(figure)
