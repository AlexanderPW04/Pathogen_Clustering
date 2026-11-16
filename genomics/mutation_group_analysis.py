"""Sprint 5 Phase 3: the "Exploratory Mutation-Frequency Groups" analysis for the website.

This runs after the Phase 1-2 depth-aware analysis and reuses its evidence
table, the Sprint 3 observation matrix, and the Phase 2 calibration. It
never changes those results. Everything it produces is EXPLORATORY: a group
is mutations with similar observed frequency patterns, not a viral genome,
haplotype, or lineage, and nothing here detects outbreaks or transmission.

Decision from the pre-registered simulation study (validation/simulate_mutation_groups.py,
docs/SPRINT5_PHASE3_IMPLEMENTATION_REPORT.md):

- The depth-aware method placed far fewer unrelated mutations in groups than
  ordinary clustering, but it did NOT meet the pre-registered limit (at most
  10% of unrelated mutations grouped) in every condition, and it recovered
  planted groups well only at high depth with little missing coverage and no
  overdispersion. By rule P6 it is therefore not offered by default.
- Default mode ("profiles"): coverage-preserving mutation profiles, the
  eligibility table, and an ungrouped frequency heatmap. No groups.
- "experimental" mode (Advanced settings): the depth-aware groups, with the
  ordinary method as a labelled comparison, both marked as not validated.
- Detection-only input (no A/C/G/T counts) never receives groups.
- "off": nothing from Phase 3 is computed.
"""

import numpy as np
import pandas as pd

from genomics import figures
from genomics import mutation_grouping as mg
from genomics.mutation_profiles import build_mutation_profiles
from genomics.unified_schema import REFERENCE_OBSERVED, VARIANT_OBSERVED

SECTION_TITLE = "Exploratory Mutation-Frequency Groups"
EXPLANATION = ("Mutations with similar frequency patterns across samples may be grouped together for comparison. "
               "These groups do not necessarily represent the same viral genome or lineage.")
LIMITATIONS = [
    "Exploratory only: a group means similar observed frequencies in the shared samples, not that the mutations "
    "sit on the same genome, form a lineage, or share ancestry. Read-level linkage is not used.",
    "Mutations are compared only over samples where both have a known value; unknown coverage is never 0.",
    "The depth-aware method assumes the stated measurement model (binomial, or beta-binomial with rho). Primer "
    "bias, amplicon dropout, and mapping errors are not modeled.",
    "Validation is on simulated mixtures with known groups only; no real wastewater data with known group "
    "structure were available.",
    "Groups are not compared with lineage definitions and must not be used to confirm outbreaks or transmission.",
]

PROFILES_MODE, EXPERIMENTAL_MODE, OFF_MODE = "profiles", "experimental", "off"
GROUPING_MODES = [PROFILES_MODE, EXPERIMENTAL_MODE, OFF_MODE]
# Measured on synthetic data (attempt 2, held-out seeds); quoted on the page.
VALIDATION_SUMMARY = (
    "Validation on simulated mixtures (29 mutations, 12 samples): with the correct noise model the depth-aware "
    "method placed 0-24% of unrelated mutations in a group, against 10-100% for ordinary clustering, and its "
    "groups were 69-100% pure. It did not meet the pre-set limit of 10% in every condition (more mutations, "
    "mixed depths, and missing coverage raised it), and it recovered planted groups well only at 100-300x with "
    "little missing coverage and no extra library variation; otherwise it mostly declined to group.")
NOT_SHOWN_REASON = ("Exploratory groups are not shown by default because the grouping method did not meet its "
                    "pre-set validation criteria (see below). The mutation profiles and eligibility table are "
                    "available; groups can be shown as an experimental view in Advanced settings.")
# Pre-registered rule P6: detection-only input failed P1 for both methods.
DEPTH_AWARE_DETECTION_ONLY_SUPPORTED = False
DETECTION_ONLY_REASON = (
    "No per-position A/C/G/T counts were supplied, so a mutation the caller did not report is unknown rather "
    "than a measured 0. Frequency patterns are then compared only over samples where both mutations were "
    "reported, which in the simulation study placed unrelated mutations in the same group too often "
    "(see the Phase 3 report). No exploratory groups are shown for this input; every other analysis is "
    "unaffected.")
MIN_MINIMUM_SHARED, MAX_MINIMUM_SHARED = mg.DEFAULT_MINIMUM_SHARED, 20
HEATMAP_ROWS = 60
TRAJECTORY_GROUPS = 6

SUMMARY_COLUMNS = [
    "Group_ID", "Method", "Model", "Assumed_Rho", "Size", "Members", "Reference_ID", "Position_Range",
    "Units_With_Two_Or_More_Members", "Median_Depth", "Count_Basis", "Homogeneity_Chi2", "Homogeneity_DF",
    "Homogeneity_P", "Max_Pairwise_Distance", "Group_VAF_Min", "Group_VAF_Max", "Group_Stability",
    "Interpretation",
]
TRAJECTORY_COLUMNS = ["Group_ID", "Method", "Unit_ID", "Collection_Date", "Location_Name", "Members",
                      "Members_Known", "Mean_Member_VAF", "Min_Member_VAF", "Max_Member_VAF", "Total_Depth"]
OVERLAP_COLUMNS = ["Group_ID", "Method", "Individual_Sample", "Group_Mutations", "Carried", "Reference_Called",
                   "Not_Callable_Or_Unknown", "Description", "Interpretation"]
INTERPRETATION = "Exploratory: similar observed frequencies only; not a genome, haplotype, or lineage."
OVERLAP_NOTE = ("Descriptive only. Carrying some or all of a group's mutations does not place this sequence in "
                "the group, and a single shared mutation is not evidence of a relationship.")


def parse_grouping_mode(text):
    value = str(text or PROFILES_MODE).strip().lower()
    if value not in GROUPING_MODES:
        raise ValueError("Choose how mutation-frequency groups are handled: profiles only, experimental groups, "
                         "or off.")
    return value


def parse_minimum_shared(text):
    """Validate the minimum-overlap setting (values below the validated minimum are refused)."""

    try:
        value = int(str(text).strip())
    except ValueError:
        raise ValueError("The minimum number of shared samples for mutation grouping must be a whole number.") from None
    if not MIN_MINIMUM_SHARED <= value <= MAX_MINIMUM_SHARED:
        raise ValueError(f"The minimum number of shared samples for mutation grouping must be from "
                         f"{MIN_MINIMUM_SHARED} (the smallest value validated in simulation) to {MAX_MINIMUM_SHARED}.")
    return value


def _info(eligibility):
    return eligibility.set_index("Variant_Key")[["Reference_ID", "Position", "Ref", "Alt"]].to_dict("index")


def _group_tables(method_label, model, rho_text, groups, keys, arrays, result, units, unit_info, stability,
                  profiles):
    """Summary and trajectory rows for named groups ({name: indices})."""

    summary, trajectory = [], []
    distance = result["pairs"]["distance"]
    for name, members in groups.items():
        member_keys = [keys[i] for i in members]
        F, N = arrays["F"][members], arrays["N"][members]
        known = np.isfinite(F)
        test = result["tests"].get(tuple(sorted(members)), (np.nan, np.nan, np.nan))
        sub = distance[np.ix_(members, members)]
        totals = np.where(known, F, 0.0).sum(axis=0)
        present_count = known.sum(axis=0)
        means = np.where(present_count > 0, totals / np.maximum(present_count, 1), np.nan)
        rows = profiles[profiles["Variant_Key"].isin(member_keys)]
        positions = rows["Position"].astype(int)
        summary.append({
            "Group_ID": name, "Method": method_label, "Model": model, "Assumed_Rho": rho_text,
            "Size": len(members), "Members": "; ".join(member_keys),
            "Reference_ID": "; ".join(sorted(set(rows["Reference_ID"].astype(str)))),
            "Position_Range": f"{positions.min()}-{positions.max()}" if len(positions) else "",
            "Units_With_Two_Or_More_Members": int((known.sum(axis=0) >= 2).sum()),
            "Median_Depth": float(np.nanmedian(N)) if np.isfinite(N).any() else np.nan,
            "Count_Basis": "; ".join(sorted(set(rows["Count_Basis"]))),
            "Homogeneity_Chi2": test[0], "Homogeneity_DF": test[1], "Homogeneity_P": test[2],
            "Max_Pairwise_Distance": float(sub.max()) if len(members) > 1 else np.nan,
            "Group_VAF_Min": float(np.nanmin(means)) if np.isfinite(means).any() else np.nan,
            "Group_VAF_Max": float(np.nanmax(means)) if np.isfinite(means).any() else np.nan,
            "Group_Stability": stability.get(name, np.nan), "Interpretation": INTERPRETATION,
        })
        for column, unit in enumerate(units):
            present = known[:, column]
            if not present.any():
                continue
            values = F[present, column]
            depth = N[present, column]
            trajectory.append({
                "Group_ID": name, "Method": method_label, "Unit_ID": unit,
                "Collection_Date": unit_info[unit]["Collection_Date"], "Location_Name": unit_info[unit]["Location_Name"],
                "Members": len(members), "Members_Known": int(present.sum()), "Mean_Member_VAF": float(values.mean()),
                "Min_Member_VAF": float(values.min()), "Max_Member_VAF": float(values.max()),
                "Total_Depth": float(np.nansum(depth)) if np.isfinite(depth).any() else np.nan})
    return summary, trajectory


def _run(kind, eligibility, profiles, units_table, unit_rho_map, rho_mode, minimum_shared, stability_allowed):
    """Run one method on its eligible (possibly size-limited) mutations."""

    column = "Eligible_Ordinary" if kind == "ordinary" else "Eligible_Depth_Aware"
    used_units = sorted(profiles["Unit_ID"].unique())
    chosen, limit_note = mg.select_mutations(eligibility, column, len(used_units))
    keys = list(chosen["Variant_Key"])
    arrays = mg.profile_arrays(profiles, keys, used_units, counts_only=(kind != "ordinary"))
    references = list(chosen["Reference_ID"].astype(str))
    if rho_mode is None:
        unit_rho = np.zeros(len(used_units))
    elif isinstance(rho_mode, float):
        unit_rho = np.full(len(used_units), rho_mode)
    else:
        unit_rho = np.array([unit_rho_map.get(unit, np.nan) for unit in used_units], dtype=float)
    if len(keys) < 2:
        result = {"groups": [], "pairs": {"distance": np.full((len(keys), len(keys)), np.inf)}, "tests": {},
                  "merges": []}
    else:
        result = mg.run_method(kind, arrays, references, unit_rho, minimum_shared)
    stability = None
    budget = len(keys) ** 2 * len(used_units) ** 2
    if stability_allowed and result["groups"] and budget <= mg.STABILITY_BUDGET:
        stability = mg.leave_one_unit_out_stability(kind, arrays, references, unit_rho, minimum_shared,
                                                    result["groups"])
    return {"keys": keys, "arrays": arrays, "result": result, "units": used_units, "unit_rho": unit_rho,
            "limit_note": limit_note, "stability": stability, "chosen": chosen}


def individual_overlap(groups_by_method, cells):
    """Descriptive table: how many of each group's mutations every individual sequence carries."""

    if cells is None or cells.empty:
        return pd.DataFrame(columns=OVERLAP_COLUMNS)
    individual = cells[cells["Sample_Type"].astype(str).eq("Individual sequence")]
    if individual.empty:
        return pd.DataFrame(columns=OVERLAP_COLUMNS)
    rows = []
    for method, groups in groups_by_method.items():
        for name, member_keys in groups.items():
            subset = individual[individual["Variant_Key"].isin(member_keys)]
            if subset.empty:
                continue
            for sample, group in subset.groupby("Sample_ID", sort=True):
                status = group["Status"].astype(str)
                carried = int(status.eq(VARIANT_OBSERVED).sum())
                reference = int(status.eq(REFERENCE_OBSERVED).sum())
                other = len(member_keys) - carried - reference
                callable_total = carried + reference
                description = ("carries every callable group mutation" if callable_total and carried == callable_total
                               else "carries some group mutations" if carried
                               else "carries none of the callable group mutations")
                rows.append({"Group_ID": name, "Method": method, "Individual_Sample": sample,
                             "Group_Mutations": len(member_keys), "Carried": carried, "Reference_Called": reference,
                             "Not_Callable_Or_Unknown": other, "Description": description,
                             "Interpretation": OVERLAP_NOTE})
    return pd.DataFrame(rows, columns=OVERLAP_COLUMNS)


def run_mutation_grouping(evidence, matrix, calibration, samples, position_table, settings):
    """Return the Phase 3 tables, figure inputs, and page summary (never raises for unsupported data)."""

    mode = settings.get("grouping_mode", PROFILES_MODE)
    enabled = mode != OFF_MODE
    experimental = mode == EXPERIMENTAL_MODE
    minimum_shared = int(settings.get("grouping_minimum_shared", mg.DEFAULT_MINIMUM_SHARED))
    hypothetical = [float(value) for value in settings.get("overdispersion_rhos", [])]
    summary = {"title": SECTION_TITLE, "explanation": EXPLANATION, "limitations": LIMITATIONS,
               "validation": VALIDATION_SUMMARY, "mode": mode, "experimental": experimental,
               "enabled": enabled, "minimum_shared": minimum_shared, "reasons": [], "methods": [],
               "groups": 0, "mutations": 0, "eligible": 0, "excluded": 0, "calibration": "",
               "primary_method": "", "comparison": []}
    tables = {}
    if not enabled:
        summary["reasons"].append("Mutation-frequency profiles and grouping were switched off in Advanced "
                                  "settings; every other analysis ran as usual.")
        return {"summary": summary, "tables": tables, "figures": {}}
    if matrix is None or matrix.empty:
        summary["reasons"].append("No quality-controlled pooled variant observations are available to group.")
        return {"summary": summary, "tables": tables, "figures": {}}

    built = build_mutation_profiles(matrix, evidence, calibration, samples, minimum_shared)
    profiles, eligibility, units_table = built["profiles"], built["eligibility"], built["units"]
    tables["sprint5_mutation_profiles.csv"] = profiles
    tables["sprint5_mutation_eligibility.csv"] = eligibility
    unit_info = units_table.set_index("Unit_ID")[["Collection_Date", "Location_Name"]].to_dict("index")
    has_base_counts = position_table is not None and position_table.attrs.get("has_base_counts", False)
    detection_only = not has_base_counts
    summary.update({"mutations": len(eligibility), "units": int(profiles["Unit_ID"].nunique()),
                    "verified_replicate_specimens": len(built["verified_pairs"]),
                    "eligible_ordinary": int(eligibility["Eligible_Ordinary"].sum()),
                    "eligible": int(eligibility["Eligible_Depth_Aware"].sum()),
                    "detection_only": detection_only})
    summary["excluded"] = summary["mutations"] - summary["eligible_ordinary"]

    used_units = set(profiles["Unit_ID"])
    calibrated = bool(used_units) and all(unit in built["rho_for_unit"] for unit in used_units)
    if calibrated:
        mode_label, rho_mode = mg.CALIBRATED_MODE, "calibrated"
        summary["calibration"] = (f"Empirically calibrated: every one of the {len(used_units)} observation units has "
                                  "a rho from its own verified technical replicates.")
    else:
        mode_label, rho_mode = mg.BINOMIAL_MODE, None
        partial = len(set(built["rho_for_unit"]) & used_units)
        summary["calibration"] = ("Binomial baseline: no empirical overdispersion estimate covers every sample, so "
                                  "independent read sampling (rho = 0) is assumed and stated."
                                  + (f" {partial} of {len(used_units)} units have a calibrated rho, which is not "
                                     "enough to call the grouping calibrated." if partial else ""))

    assignments, group_summary, trajectories, comparison = [], [], [], []
    named_groups, runs = {}, {}
    plan = [("ordinary", mg.ORDINARY_METHOD, "Frequencies only (depth ignored)", None, "O"),
            ("depth", mg.DEPTH_AWARE_METHOD, mode_label, rho_mode, "G")]
    known_counts = eligibility.set_index("Variant_Key")
    if not experimental:
        summary["reasons"].append(NOT_SHOWN_REASON)
        plan = []
    elif detection_only and not DEPTH_AWARE_DETECTION_ONLY_SUPPORTED:
        summary["reasons"].append(DETECTION_ONLY_REASON)
        plan = []
    elif eligibility["Eligible_Depth_Aware"].sum() < 2:
        # The ordinary method is shown only as a comparison for the depth-aware one.
        plan = [item for item in plan if item[0] == "depth"]
    for kind, label, model, rho_setting, prefix in plan:
        column = "Eligible_Ordinary" if kind == "ordinary" else "Eligible_Depth_Aware"
        if eligibility[column].sum() < 2:
            summary["methods"].append({"method": label, "model": model, "ran": False, "groups": 0, "grouped": 0,
                                       "analyzed": int(eligibility[column].sum()),
                                       "reason": (f"Fewer than two mutations have known values in at least "
                                                  f"{minimum_shared} observation units"
                                                  + (" with read counts" if kind != "ordinary" else "")
                                                  + ", so no pair can be compared.")})
            continue
        run = _run(kind, eligibility, profiles, units_table, built["rho_for_unit"], rho_setting, minimum_shared, True)
        runs[kind] = run
        known = known_counts["Known_Units" if kind == "ordinary" else "Known_Units_With_Counts"]
        rho_text = "calibrated per unit" if rho_setting == "calibrated" else ("0" if kind != "ordinary" else "")
        table, groups = mg.assignment_table(run["keys"], _info(eligibility), run["result"], label, model, rho_text,
                                            prefix, known, run["stability"])
        assignments.append(table)
        stability = {}
        if run["stability"] is not None:
            original = [sorted(g) for g in run["result"]["groups"]]
            stability = {name: run["stability"][original.index(sorted(g))] for name, g in groups.items()}
        rows, traj = _group_tables(label, model, rho_text, groups, run["keys"], run["arrays"], run["result"],
                                   run["units"], unit_info, stability, profiles)
        group_summary += rows
        trajectories += traj
        named_groups[label] = {name: [run["keys"][i] for i in g] for name, g in groups.items()}
        grouped = int(table["Assignment_Status"].eq(mg.GROUPED).sum())
        reason = ""
        if not groups:
            reason = ("No set of comparable mutations passed the grouping criterion, so no group is reported "
                      "(this is a valid result, not an error).")
            if kind != "ordinary" and table["Comparable_Partners"].eq(0).all():
                reason = ("No pair of mutations had enough shared, deep-enough observations to detect a "
                          f"{mg.RESOLUTION:.2f} frequency difference, so the depth-aware method abstains.")
        summary["methods"].append({"method": label, "model": model, "ran": True, "groups": len(groups),
                                   "grouped": grouped, "analyzed": len(run["keys"]), "reason": reason,
                                   "limit_note": run["limit_note"]})
        if run["limit_note"]:
            summary["reasons"].append(run["limit_note"])
        if run["limit_note"]:
            excluded = eligibility[eligibility[column] & ~eligibility["Variant_Key"].isin(run["keys"])]
            assignments.append(pd.DataFrame({
                "Variant_Key": excluded["Variant_Key"], "Reference_ID": excluded["Reference_ID"],
                "Position": excluded["Position"], "Ref": excluded["Ref"], "Alt": excluded["Alt"], "Method": label,
                "Model": model, "Assumed_Rho": rho_text, "Group_ID": "", "Group_Size": 1,
                "Assignment_Status": mg.NOT_IN_SUBSET, "Comparable_Partners": np.nan, "Nearest_Partner": "",
                "Nearest_Distance": np.nan, "Known_Units": excluded["Known_Units"], "Group_Stability": np.nan}))

    # Comparison: ordinary versus depth-aware, and hypothetical rho scenarios (sensitivity only).
    if "depth" in runs:
        primary = runs["depth"]
        primary_labels = mg.label_vector(primary["result"]["groups"], len(primary["keys"]))
        if "ordinary" in runs:
            ordinary = runs["ordinary"]
            position = {key: i for i, key in enumerate(ordinary["keys"])}
            common = [key for key in primary["keys"] if key in position]
            if len(common) >= 2:
                ord_labels = mg.label_vector(ordinary["result"]["groups"], len(ordinary["keys"]))
                a = [primary_labels[primary["keys"].index(k)] for k in common]
                b = [ord_labels[position[k]] for k in common]
                comparison.append({"Comparison": "Ordinary vs depth-aware (primary)", "Mutations_Compared": len(common),
                                   "Adjusted_Rand_Index": mg.adjusted_rand_index(a, b),
                                   "Groups_First": len(ordinary["result"]["groups"]),
                                   "Groups_Second": len(primary["result"]["groups"]),
                                   "Note": "Agreement between two exploratory methods, not accuracy."})
        for rho in hypothetical:
            run = _run("depth", eligibility, profiles, units_table, {}, float(rho), minimum_shared, False)
            labels = mg.label_vector(run["result"]["groups"], len(run["keys"]))
            grouped = sum(len(g) for g in run["result"]["groups"])
            comparison.append({"Comparison": f"Hypothetical rho {rho:g} vs primary",
                               "Mutations_Compared": len(run["keys"]),
                               "Adjusted_Rand_Index": mg.adjusted_rand_index(primary_labels, labels),
                               "Groups_First": len(primary["result"]["groups"]),
                               "Groups_Second": len(run["result"]["groups"]),
                               "Note": (f"Hypothetical; not calibrated. {grouped} mutation(s) grouped. A larger assumed "
                                        "rho widens what counts as similar, so it can join unrelated mutations.")})
    if assignments:
        tables["sprint5_mutation_groups.csv"] = pd.concat(assignments, ignore_index=True)
    if experimental and runs:
        tables["sprint5_mutation_group_summary.csv"] = pd.DataFrame(group_summary, columns=SUMMARY_COLUMNS)
    if trajectories:
        tables["sprint5_mutation_group_trajectories.csv"] = pd.DataFrame(trajectories, columns=TRAJECTORY_COLUMNS)
    if comparison:
        tables["sprint5_grouping_comparison.csv"] = pd.DataFrame(comparison)
    tables["sprint5_grouping_parameters.csv"] = parameter_table(summary, minimum_shared, mode_label, hypothetical,
                                                                built, detection_only)
    depth_ran = "depth" in runs
    summary["primary_method"] = mg.DEPTH_AWARE_METHOD if depth_ran else ""
    summary["groups"] = len(runs["depth"]["result"]["groups"]) if depth_ran else 0
    summary["comparison"] = comparison
    figure_inputs = {"runs": runs, "profiles": profiles, "units": units_table, "named_groups": named_groups}
    return {"summary": summary, "tables": tables, "figures": figure_inputs, "named_groups": named_groups}


def parameter_table(summary, minimum_shared, mode_label, hypothetical, built, detection_only):
    rows = [
        ("Grouping mode", summary["mode"], "profiles = default (no groups); experimental = groups shown as not "
         "validated; off = Phase 3 skipped."),
        ("Validation status", "Did not meet the pre-registered criteria (P1, P2); see the Phase 3 report",
         "Groups are therefore never shown by default."),
        ("Minimum shared observation units", minimum_shared,
         "Two mutations are compared only over units where both are known, and only with at least this many."),
        ("Ordinary method", mg.ORDINARY_METHOD, "Baseline; ignores read depth."),
        ("Ordinary cut height (RMS frequency difference)", mg.ORDINARY_HEIGHT,
         "Every pair in an ordinary group differs by at most this much (complete linkage)."),
        ("Depth-aware method", mg.DEPTH_AWARE_METHOD, "MAGOS-inspired; tumor-specific assumptions removed."),
        ("Depth-aware measurement model", mode_label, "Binomial and hypothetical modes are not empirically calibrated."),
        ("Resolution (frequency difference to detect)", mg.RESOLUTION, "Power gate for comparable pairs."),
        ("Required power", mg.POWER, "Pairs below it are not comparable (abstain rather than group)."),
        ("Test level alpha", mg.ALPHA, "A clade is accepted when the homogeneity test gives P >= alpha."),
        ("Maximum mutations per method", mg.MAX_MUTATIONS, "Explicit limit for the pairwise step; disclosed if used."),
        ("Hypothetical rho values (sensitivity only)", ", ".join(f"{v:g}" for v in hypothetical) or "none",
         "Scenarios, never estimates."),
        ("Verified technical-replicate specimens", len(built["verified_pairs"]),
         "Specimen_ID + distinct Library_ID; their libraries are combined into one unit."),
        ("Units with a calibrated rho", len(built["rho_for_unit"]), "Upper 95% bound from Phase 2."),
        ("Supported non-detections available (A/C/G/T counts)", "no" if detection_only else "yes",
         "Without them an unreported mutation is unknown, not 0."),
        ("Mutation identity", "Reference_ID (segment) + position + Ref + Alt",
         "Different references or segments are never compared."),
    ]
    return pd.DataFrame(rows, columns=["Parameter", "Value", "Explanation"])


def _unit_order(units, units_table):
    info = units_table.set_index("Unit_ID")
    def key(unit):
        date = info.loc[unit, "Collection_Date"] if unit in info.index else pd.NaT
        return (pd.isna(date), date if pd.notna(date) else pd.Timestamp(0), str(info.loc[unit, "Location_Name"])
                if unit in info.index else "", unit)
    ordered = sorted(units, key=key)
    labels = []
    for unit in ordered:
        date = info.loc[unit, "Collection_Date"]
        labels.append((f"{pd.Timestamp(date).date()} " if pd.notna(date) else "") + f"{unit}")
    return ordered, labels


def add_grouping_outputs(result, evidence, matrix, calibration, samples, position_table, settings, add_figure, note):
    """Run the Phase 3 analysis and attach its tables, figures, and summary to `result`."""

    try:
        output = run_mutation_grouping(evidence, matrix, calibration, samples, position_table, settings)
    except (ValueError, KeyError, IndexError) as problem:  # never let an exploratory step stop the analysis
        note("depth_aware", f"Exploratory mutation grouping could not run: {problem}")
        return {"title": SECTION_TITLE, "explanation": EXPLANATION, "limitations": LIMITATIONS, "enabled": True,
                "reasons": [f"The grouping step could not run on this input ({problem}); every other analysis is "
                            "unaffected."], "methods": [], "groups": 0, "mutations": 0, "eligible": 0,
                "excluded": 0, "calibration": "", "primary_method": "", "comparison": [], "error": True}
    result["tables"].update(output["tables"])
    summary = output["summary"]
    runs = output["figures"].get("runs", {}) if output["figures"] else {}
    result["mutation_groups"] = output.get("named_groups", {})
    primary_kind = "depth" if "depth" in runs else ("ordinary" if "ordinary" in runs else None)
    if primary_kind is None:
        add_profile_heatmap(result, output, add_figure)
        return summary
    run = runs[primary_kind]
    units_table = output["figures"]["units"]
    groups = sorted(run["result"]["groups"], key=lambda g: (-len(g), min(g)))
    prefix = "G" if primary_kind == "depth" else "O"
    group_name = {index: f"{prefix}{number}" for number, group in enumerate(groups, start=1) for index in group}
    grouped_rows = [index for group in groups for index in group]
    known = np.isfinite(run["arrays"]["F"]).sum(axis=1)
    others = [index for index in np.argsort(-known, kind="mergesort") if index not in group_name]
    rows = (grouped_rows + others)[:HEATMAP_ROWS]
    if groups:
        ordered, labels = _unit_order(run["units"], units_table)
        columns = [run["units"].index(unit) for unit in ordered]
        F = run["arrays"]["F"][np.ix_(rows, columns)]
        shown_grouped = sum(1 for index in rows if index in group_name)
        method_name = "depth-aware" if primary_kind == "depth" else "ordinary"
        add_figure(result, "mutation_group_heatmap", "Exploratory mutation-frequency groups (heatmap)",
                   f"Rows: {len(rows)} of the {len(run['keys'])} mutations analyzed by the {method_name} method "
                   f"({shown_grouped} of {len(grouped_rows)} grouped mutations, then ungrouped ones with the most "
                   "known values); the full assignment table is sprint5_mutation_groups.csv. Columns: observation "
                   "units by date. Hatched gray cells are unknown, never 0. Groups are similar observed frequency "
                   "patterns, not lineages.",
                   lambda: figures.mutation_group_heatmap_figure(
                       F, [run["keys"][i] for i in rows], [group_name.get(i, "") for i in rows], labels))
        trajectories = result["tables"].get("sprint5_mutation_group_trajectories.csv")
        chosen = [f"{prefix}{number}" for number in range(1, min(len(groups), TRAJECTORY_GROUPS) + 1)]
        method_label = mg.DEPTH_AWARE_METHOD if primary_kind == "depth" else mg.ORDINARY_METHOD
        if trajectories is not None and trajectories["Collection_Date"].notna().any():
            subset = trajectories[trajectories["Method"].eq(method_label)]
            add_figure(result, "mutation_group_trajectories", "Exploratory groups over time",
                       f"The {len(chosen)} largest of {len(groups)} groups (all are in "
                       "sprint5_mutation_group_trajectories.csv). Each point is the mean frequency of the group's "
                       "mutations known in that sample; the bar spans its lowest to highest member. Unknown "
                       "members are left out, not counted as 0. Descriptive only.",
                       lambda: figures.group_trajectory_figure(subset, chosen))
    if "depth" in runs and "ordinary" in runs:
        depth_run, ordinary_run = runs["depth"], runs["ordinary"]
        depth_names = _names(depth_run, "G")
        ordinary_names = _names(ordinary_run, "O")
        common = [key for key in depth_run["keys"] if key in set(ordinary_run["keys"])]
        if common and (depth_run["result"]["groups"] or ordinary_run["result"]["groups"]):
            add_figure(result, "grouping_comparison", "Ordinary versus depth-aware grouping",
                       f"The {len(common)} mutations analyzed by both methods. Each cell counts mutations placed in "
                       "that depth-aware group (rows) and that ordinary group (columns). Agreement between two "
                       "exploratory methods is not evidence that either is correct; in simulation the ordinary "
                       "method grouped unrelated mutations far more often (see the Phase 3 report).",
                       lambda: figures.grouping_comparison_figure(
                           [depth_names.get(k, "—") for k in common], [ordinary_names.get(k, "—") for k in common]))
    return summary


def _names(run, prefix):
    groups = sorted(run["result"]["groups"], key=lambda g: (-len(g), min(g)))
    return {run["keys"][index]: f"{prefix}{number}" for number, group in enumerate(groups, start=1) for index in group}


def add_individual_overlap(result):
    """After the cross-source step: a descriptive Sprint 2 comparison of each group (if both sources exist)."""

    groups = result.get("mutation_groups")
    cells = result["tables"].get("cross_source_evidence_by_sample.csv")
    if not groups or cells is None or cells.empty:
        return None
    table = individual_overlap(groups, cells)
    if not table.empty:
        result["tables"]["sprint5_group_individual_overlap.csv"] = table
    return table


def add_profile_heatmap(result, output, add_figure):
    """Ungrouped heatmap of the eligible mutations with the most known values (profiles mode)."""

    figure_inputs = output.get("figures") or {}
    if "profiles" not in figure_inputs:
        return
    eligibility = output["tables"]["sprint5_mutation_eligibility.csv"]
    eligible = eligibility[eligibility["Eligible_Ordinary"]]
    if len(eligible) < 2:
        return
    chosen = eligible.sort_values(["Known_Units", "Reference_ID", "Position", "Alt"],
                                  ascending=[False, True, True, True], kind="mergesort").head(HEATMAP_ROWS)
    chosen = chosen.sort_values(["Reference_ID", "Position", "Alt"], kind="mergesort")
    profiles, units_table = figure_inputs["profiles"], figure_inputs["units"]
    units = sorted(profiles.loc[profiles["Variant_Key"].isin(chosen["Variant_Key"]), "Unit_ID"].unique())
    ordered, labels = _unit_order(units, units_table)
    arrays = mg.profile_arrays(profiles, list(chosen["Variant_Key"]), ordered, counts_only=False)
    F = arrays["F"]
    add_figure(result, "mutation_profile_heatmap", "Mutation-frequency profiles (not grouped)",
                  f"{len(chosen)} of the {len(eligible)} mutations known in at least "
                  f"{output['summary']['minimum_shared']} observation units (those with the most known values), in "
                  "genome order; every profile is in sprint5_mutation_profiles.csv. Columns are observation units by "
                  "date. Hatched gray cells are unknown, never 0. No groups are drawn: grouping did not meet its "
                  "validation criteria and is only available as an experimental view.",
                  lambda: figures.mutation_group_heatmap_figure(F, list(chosen["Variant_Key"]), [""] * len(chosen),
                                                                labels))
