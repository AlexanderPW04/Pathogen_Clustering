"""Sprint 5 Phase 3: exploratory grouping of mutations with similar frequency patterns.

EXPLORATORY. A group is "mutations whose observed frequencies are similar in
every shared observation unit". It is NOT a viral genome, haplotype, or
lineage, and it says nothing about linkage on the same molecule, ancestry, or
transmission. Inspired by MAGOS (Ahmadinejad et al. 2022, Mol Biol Evol
39:msac136), with the tumor-specific parts removed (see
docs/SPRINT5_PHASE3_IMPLEMENTATION_REPORT.md).

Inputs are the long profiles from `genomics/mutation_profiles.py`. Unknown
values stay unknown: two mutations are compared only over the units where
BOTH have a known value, and only when they share at least `minimum_shared`
units. A pair with too little overlap has NO distance (infinity), never
distance 0.

Linkage
-------
Complete linkage (the distance between two clusters is the largest pairwise
distance). It needs no Euclidean property, unlike Ward or centroid linkage,
so it is valid for distances computed over different shared units. Because no
pair can be joined across an infinite distance, every pair inside a group is
comparable. The implementation is a deterministic O(V^2)-per-merge NumPy loop
(SciPy is not a dependency); ties go to the lowest index.

1. Ordinary clustering (baseline; frequencies only, depth ignored)
   d_ij = sqrt( mean over shared units s of (f_is - f_js)^2 ), the RMS
   frequency difference. Groups are the clusters after joining every merge at
   height <= ORDINARY_HEIGHT (0.10): complete linkage guarantees every pair in
   a group differs by at most 0.10 RMS.

2. Depth-aware grouping (MAGOS-inspired; needs read counts)
   Measurement model per unit: VAF variance p (1 - p) c, with
   c = [n + rho (sum n_l^2 - n)] / n^2 (beta-binomial, Phase 2; c = 1 / n
   for one library with rho = 0, the binomial baseline).

   Pair statistic: in each shared unit s, under "same underlying frequency",
       z_s^2 = (f_is - f_js)^2 / [ p_s (1 - p_s) (c_is + c_js) ],
   with p_s the inverse-variance-weighted pooled frequency. Units where both
   values are exactly 0 (or 1) carry no variance and are dropped, which can
   only make a difference easier to detect. X2_ij = sum z_s^2 is a Pearson
   chi-square of homogeneity (2 x 2 per unit at rho = 0) with df = the number
   of units kept. The linkage distance is X2_ij / df: about 1 when the pair
   is consistent with one frequency and large when it is not.

   Power gate (added; not in MAGOS): not rejecting "same frequency" is weak
   evidence when depth is low. A pair is comparable only if the shared units
   would give at least POWER (80%) probability of detecting a difference of
   RESOLUTION (0.10) in every unit at level ALPHA. The noncentrality is
   lam_ij = sum_s RESOLUTION^2 / [ (c_is + c_js) p* (1 - p*) ], with p* the
   frequency within RESOLUTION / 2 of p_s that has the largest variance
   (conservative). A pair is comparable iff lam_ij >= lam_80(k_ij).

   Adaptive partitioning (MAGOS-style, re-derived): walk each tree of the
   forest from the root. Accept a node as a group if the design-effect-
   corrected Pearson chi-square test of equal frequencies across all its
   members (X2_G = sum_s sum_i (f_is - p_Gs)^2 / [p_Gs (1 - p_Gs) c_is],
   df = sum_s (m_s - 1) over units with two or more known members) gives
   P >= ALPHA. Otherwise split it into its two children. Leaves are single,
   ungrouped mutations. MAGOS's threshold was NOT transferred: false grouping
   is measured in validation/simulate_mutation_groups.py.

   Informative-pattern requirement (attempt 2, after attempt 1 failed the
   pre-registered false-grouping criterion): "not rejected" also accepts
   unrelated mutations that merely sit at the same, nearly constant level,
   because differences smaller than RESOLUTION cannot be resolved. Agreement
   is evidence of a shared PATTERN only if the pattern varies. A clade is
   therefore accepted only if, across the units where at least two members
   are known, its pooled frequency spans at least PATTERN_RANGE (0.20, two
   resolution steps). Groups whose frequencies are flat are left ungrouped
   (a recall cost accepted for fewer false groups).

Calibration modes
-----------------
- Empirically calibrated: every unit in the analysis has a rho estimated from
  its own verified technical replicates (Phase 2, upper 95% bound).
- Binomial baseline: rho = 0, independent reads; stated as an assumption.
- Hypothetical: a user-chosen rho for every unit, a sensitivity scenario only.
"""

import math

import numpy as np
import pandas as pd

from genomics.chi_square import chi2_sf, noncentrality_for_power
from genomics.mutation_profiles import FREQUENCY_ONLY, variance_factor

ORDINARY_METHOD = "Ordinary hierarchical clustering (RMS frequency difference, complete linkage)"
DEPTH_AWARE_METHOD = "Depth-aware grouping (beta-binomial chi-square, complete linkage, adaptive partitioning)"
CALIBRATED_MODE = "Empirically calibrated (rho from technical replicates)"
BINOMIAL_MODE = "Binomial baseline (independent reads assumed; rho = 0)"
HYPOTHETICAL_MODE = "Hypothetical sensitivity (assumed rho; not calibrated)"

DEFAULT_MINIMUM_SHARED = 2      # chosen by the pre-registered rule in the simulation study (Phase 3 report)
ORDINARY_HEIGHT = 0.10          # RMS frequency difference allowed within an ordinary group
RESOLUTION = 0.10               # frequency difference the depth-aware test must be able to detect ...
POWER = 0.80                    # ... with at least this power ...
ALPHA = 0.01                    # ... at this level (also the acceptance level, as in MAGOS)
PATTERN_RANGE = 0.20            # attempt 2: a group's shared profile must vary at least this much
MAX_MUTATIONS = 1000            # explicit size limit for the O(V^2) pairwise step
PAIR_BUDGET = 5.0e7             # V^2 x units above this is reduced to a disclosed subset (about 1-2 s)
STABILITY_BUDGET = 2.0e8        # V^2 x units^2 above this skips the leave-one-unit-out check

GROUPED = "Grouped"
NO_PARTNER = "Ungrouped: no comparable mutation"
NOT_SIMILAR = "Ungrouped: no comparable mutation with a consistent frequency pattern"
NOT_IN_SUBSET = "Not analyzed: outside the size-limited subset"

ASSIGNMENT_COLUMNS = [
    "Variant_Key", "Reference_ID", "Position", "Ref", "Alt", "Method", "Model", "Assumed_Rho",
    "Group_ID", "Group_Size", "Assignment_Status", "Comparable_Partners", "Nearest_Partner",
    "Nearest_Distance", "Known_Units", "Group_Stability",
]


# ----------------------------------------------------------------------------
# Arrays
# ----------------------------------------------------------------------------

def profile_arrays(profiles, keys, units, counts_only):
    """Dense V x U arrays for the selected mutations (NaN = unknown)."""

    row = {key: index for index, key in enumerate(keys)}
    column = {unit: index for index, unit in enumerate(units)}
    table = profiles[profiles["Variant_Key"].isin(row) & profiles["Unit_ID"].isin(column)]
    if counts_only:
        table = table[~table["Count_Basis"].eq(FREQUENCY_ONLY)]
    shape = (len(keys), len(units))
    arrays = {name: np.full(shape, np.nan) for name in ["F", "X", "N", "Q"]}
    i = table["Variant_Key"].map(row).to_numpy()
    j = table["Unit_ID"].map(column).to_numpy()
    arrays["F"][i, j] = table["VAF"].to_numpy(dtype=float)
    if counts_only:
        arrays["X"][i, j] = table["Alt_Read_Count"].astype("Float64").to_numpy(dtype=float, na_value=np.nan)
        arrays["N"][i, j] = table["Total_Depth"].astype("Float64").to_numpy(dtype=float, na_value=np.nan)
        arrays["Q"][i, j] = table["Depth_Square_Sum"].to_numpy(dtype=float)
    return arrays


def _same_reference(references):
    labels = np.asarray(references)
    return labels[:, None] == labels[None, :]


# ----------------------------------------------------------------------------
# Distances
# ----------------------------------------------------------------------------

def ordinary_distances(F, references, minimum_shared):
    """RMS difference over shared known units; inf with fewer than `minimum_shared` shared units."""

    known = np.isfinite(F).astype(float)
    values = np.where(known > 0, F, 0.0)
    shared = known @ known.T
    squares = (values ** 2) @ known.T
    total = squares + squares.T - 2 * (values @ values.T)
    with np.errstate(divide="ignore", invalid="ignore"):
        distance = np.sqrt(np.clip(total, 0, None) / shared)
    comparable = (shared >= minimum_shared) & _same_reference(references)
    distance = np.where(comparable, distance, np.inf)
    np.fill_diagonal(distance, 0.0)
    return {"distance": distance, "shared": shared.astype(int), "comparable": comparable}


def _p_star(pooled, resolution):
    """Frequency within resolution / 2 of the pooled value with the largest p (1 - p)."""

    return np.clip(np.clip(0.5, pooled - resolution / 2, pooled + resolution / 2), 0.0, 1.0)


def depth_aware_distances(arrays, references, unit_rho, minimum_shared, resolution=RESOLUTION,
                          alpha=ALPHA, power=POWER):
    """Pairwise chi-square distance X2/df with the minimum-overlap and power gates."""

    F, N, Q = arrays["F"], arrays["N"], arrays["Q"]
    count, units = F.shape
    chi2 = np.zeros((count, count))
    df = np.zeros((count, count))
    shared = np.zeros((count, count))
    lam = np.zeros((count, count))
    for s in range(units):
        known = np.isfinite(F[:, s]) & np.isfinite(N[:, s]) & (N[:, s] > 0)
        index = np.flatnonzero(known)
        if len(index) < 2:
            continue
        f = F[index, s]
        c = variance_factor(N[index, s], Q[index, s], unit_rho[s])
        weight = 1.0 / c
        pooled = (weight[:, None] * f[:, None] + weight[None, :] * f[None, :]) / (weight[:, None] + weight[None, :])
        spread = c[:, None] + c[None, :]
        variance = pooled * (1 - pooled) * spread
        positive = variance > 0
        with np.errstate(divide="ignore", invalid="ignore"):
            term = np.where(positive, (f[:, None] - f[None, :]) ** 2 / variance, 0.0)
        star = _p_star(pooled, resolution)
        block = np.ix_(index, index)
        chi2[block] += term
        df[block] += positive
        shared[block] += 1
        lam[block] += resolution ** 2 / (spread * star * (1 - star))
    needed = np.full(count and int(shared.max()) + 1, np.inf)
    for k in range(minimum_shared, len(needed)):
        needed[k] = noncentrality_for_power(k, alpha, power)
    k = shared.astype(int)
    overlap = shared >= minimum_shared
    resolvable = overlap & (lam >= needed[k])
    comparable = resolvable & _same_reference(references)
    with np.errstate(divide="ignore", invalid="ignore"):
        distance = np.where(comparable, chi2 / np.maximum(df, 1), np.inf)
    np.fill_diagonal(distance, 0.0)
    return {"distance": distance, "shared": k, "comparable": comparable, "overlap": overlap,
            "resolvable": resolvable, "chi2": chi2, "df": df.astype(int), "noncentrality": lam}


# ----------------------------------------------------------------------------
# Complete linkage
# ----------------------------------------------------------------------------

def complete_linkage(distance, stop_above=np.inf):
    """Agglomerate with complete linkage; return merges [(node_a, node_b, height, new_node)].

    Leaves are 0..V-1 and merge t creates node V + t. Merging stops at the first
    height that is infinite (or above `stop_above`), so the result is a forest.
    """

    matrix = np.array(distance, dtype=float, copy=True)
    count = len(matrix)
    np.fill_diagonal(matrix, np.inf)
    node = list(range(count))
    merges = []
    for step in range(max(count - 1, 0)):
        flat = int(np.argmin(matrix))
        i, j = divmod(flat, count)
        height = matrix[i, j]
        if not np.isfinite(height) or height > stop_above:
            break
        if j < i:
            i, j = j, i
        merged = np.maximum(matrix[i], matrix[j])
        matrix[i, :] = merged
        matrix[:, i] = merged
        matrix[i, i] = np.inf
        matrix[j, :] = np.inf
        matrix[:, j] = np.inf
        merges.append((node[i], node[j], float(height), count + step))
        node[i] = count + step
    return merges


def tree_members(merges, count):
    """Members, children, and height of every node of the forest; also its roots."""

    members = {leaf: [leaf] for leaf in range(count)}
    children, height = {}, {}
    used = set()
    for a, b, h, new in merges:
        members[new] = members[a] + members[b]
        children[new] = (a, b)
        height[new] = h
        used.update((a, b))
    roots = [node for node in members if node not in used]
    return members, children, height, roots


# ----------------------------------------------------------------------------
# Group test
# ----------------------------------------------------------------------------

def homogeneity_test(members, arrays, unit_rho):
    """Design-effect-corrected Pearson chi-square of equal frequencies across a set of mutations."""

    F, N, Q = arrays["F"][members], arrays["N"][members], arrays["Q"][members]
    known = np.isfinite(F) & np.isfinite(N) & (N > 0)
    c = variance_factor(np.where(known, N, 1.0), np.where(known, Q, 1.0), np.asarray(unit_rho)[None, :])
    weight = np.where(known, 1.0 / c, 0.0)
    total = weight.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        pooled = (weight * np.where(known, F, 0.0)).sum(axis=0) / total
    present = known.sum(axis=0)
    usable = (present >= 2) & (pooled > 0) & (pooled < 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(known, weight * (np.where(known, F, 0.0) - pooled[None, :]) ** 2, 0.0).sum(axis=0) / (
            pooled * (1 - pooled))
    statistic = float(np.where(usable, terms, 0.0).sum())
    df = int((present[usable] - 1).sum())
    p_value = chi2_sf(statistic, df) if df > 0 else float("nan")
    return statistic, df, p_value


# ----------------------------------------------------------------------------
# Group selection
# ----------------------------------------------------------------------------

def pattern_range(members, arrays, unit_rho=None):
    """Range of the group's pooled frequency over units where at least two members are known.

    With counts the pooled value is inverse-variance weighted (as in the test);
    for frequency-only arrays it is the mean member frequency.
    """

    F = arrays["F"][members]
    known = np.isfinite(F)
    if unit_rho is not None and "N" in arrays and np.isfinite(arrays["N"][members][known]).all():
        N, Q = arrays["N"][members], arrays["Q"][members]
        c = variance_factor(np.where(known, N, 1.0), np.where(known, Q, 1.0), np.asarray(unit_rho)[None, :])
        weight = np.where(known, 1.0 / c, 0.0)
    else:
        weight = known.astype(float)
    shared = known.sum(axis=0) >= 2
    if not shared.any():
        return 0.0
    total = weight.sum(axis=0)
    pooled = (weight * np.where(known, F, 0.0)).sum(axis=0)[shared] / total[shared]
    return float(pooled.max() - pooled.min())


def top_down(merges, count, accept):
    """Walk each tree from its root; keep the first node `accept(members, node)` approves, else split it."""

    members, children, height, roots = tree_members(merges, count)
    groups = []
    stack = sorted(roots, reverse=True)
    while stack:
        current = stack.pop()
        group = members[current]
        if len(group) < 2:
            continue
        if accept(sorted(group), height.get(current, 0.0)):
            groups.append(sorted(group))
        else:
            stack.extend(children[current])
    return groups


def ordinary_groups(distance, height=ORDINARY_HEIGHT, arrays=None, minimum_range=None):
    """Clusters after every complete-linkage merge at or below `height` (lists of indices).

    With `minimum_range` (simulation ablation only) a cluster must also pass the
    informative-pattern requirement, otherwise it is split.
    """

    merges = complete_linkage(distance, stop_above=height)
    if minimum_range is None:
        members, _, _, roots = tree_members(merges, len(distance))
        return [sorted(members[root]) for root in roots if len(members[root]) > 1], merges
    groups = top_down(merges, len(distance),
                      lambda group, _: pattern_range(group, arrays) >= minimum_range)
    return groups, merges


def adaptive_groups(distance, arrays, unit_rho, alpha=ALPHA, minimum_range=PATTERN_RANGE):
    """MAGOS-style top-down acceptance on the complete-linkage forest; returns (groups, tests, merges)."""

    merges = complete_linkage(distance)
    tests = {}

    def accept(group, _):
        statistic, df, p_value = homogeneity_test(group, arrays, unit_rho)
        if not (df > 0 and p_value >= alpha):
            return False
        if minimum_range is not None and pattern_range(group, arrays, unit_rho) < minimum_range:
            return False
        tests[tuple(group)] = (statistic, df, p_value)
        return True

    groups = top_down(merges, len(distance), accept)
    return groups, tests, merges


def label_vector(groups, count):
    """Integer labels: each group shares a label; every other mutation is its own singleton."""

    labels = np.arange(count) + 1_000_000
    for number, group in enumerate(groups):
        labels[list(group)] = number
    return labels


# ----------------------------------------------------------------------------
# Running one configuration
# ----------------------------------------------------------------------------

def run_method(method, arrays, references, unit_rho, minimum_shared, alpha=ALPHA, minimum_range=PATTERN_RANGE):
    """Return {"groups", "pairs", "tests", "merges"} for one method on prepared arrays.

    `minimum_range` is the informative-pattern requirement of the depth-aware
    method (None reproduces attempt 1). The ordinary baseline never uses it;
    "ordinary_pattern" is the simulation-only ablation that adds it.
    """

    if method in ("ordinary", "ordinary_pattern"):
        pairs = ordinary_distances(arrays["F"], references, minimum_shared)
        groups, merges = ordinary_groups(pairs["distance"], arrays=arrays,
                                         minimum_range=PATTERN_RANGE if method == "ordinary_pattern" else None)
        return {"groups": groups, "pairs": pairs, "tests": {}, "merges": merges}
    pairs = depth_aware_distances(arrays, references, unit_rho, minimum_shared, alpha=alpha)
    groups, tests, merges = adaptive_groups(pairs["distance"], arrays, unit_rho, alpha, minimum_range)
    return {"groups": groups, "pairs": pairs, "tests": tests, "merges": merges}


def leave_one_unit_out_stability(method, arrays, references, unit_rho, minimum_shared, groups,
                                 minimum_range=PATTERN_RANGE):
    """For each group, the share of leave-one-unit-out reruns that keep all its members together."""

    if not groups:
        return []
    units = arrays["F"].shape[1]
    kept = np.zeros(len(groups))
    for drop in range(units):
        keep = [index for index in range(units) if index != drop]
        reduced = {name: values[:, keep] for name, values in arrays.items()}
        result = run_method(method, reduced, references, np.asarray(unit_rho)[keep], minimum_shared,
                            minimum_range=minimum_range)
        labels = label_vector(result["groups"], len(references))
        for number, group in enumerate(groups):
            kept[number] += len(set(labels[group])) == 1
    return list(kept / units)


# ----------------------------------------------------------------------------
# Agreement metrics (also used by the simulation study)
# ----------------------------------------------------------------------------

def _pairs(values):
    values = np.asarray(values, dtype=float)
    return float((values * (values - 1) / 2).sum())


def adjusted_rand_index(truth, predicted):
    """Hubert-Arabie adjusted Rand index; 1.0 for two identical trivial partitions (as scikit-learn)."""

    truth, predicted = np.asarray(truth), np.asarray(predicted)
    count = len(truth)
    if count < 2:
        return 1.0
    _, t = np.unique(truth, return_inverse=True)
    _, p = np.unique(predicted, return_inverse=True)
    table = np.zeros((t.max() + 1, p.max() + 1))
    np.add.at(table, (t, p), 1)
    index = _pairs(table)
    rows, columns = _pairs(table.sum(axis=1)), _pairs(table.sum(axis=0))
    expected = rows * columns / _pairs([count])
    maximum = (rows + columns) / 2
    if maximum == expected:
        return 1.0
    return (index - expected) / (maximum - expected)


def normalized_mutual_information(truth, predicted):
    """Mutual information over the arithmetic mean of the entropies (scikit-learn's default)."""

    truth, predicted = np.asarray(truth), np.asarray(predicted)
    count = len(truth)
    _, t = np.unique(truth, return_inverse=True)
    _, p = np.unique(predicted, return_inverse=True)
    table = np.zeros((t.max() + 1, p.max() + 1))
    np.add.at(table, (t, p), 1)
    joint = table / count
    pt, pp = joint.sum(axis=1), joint.sum(axis=0)
    entropy_t = -float((pt[pt > 0] * np.log(pt[pt > 0])).sum())
    entropy_p = -float((pp[pp > 0] * np.log(pp[pp > 0])).sum())
    nonzero = joint > 0
    information = float((joint[nonzero] * np.log(joint[nonzero] / np.outer(pt, pp)[nonzero])).sum())
    if entropy_t == 0 and entropy_p == 0:
        return 1.0
    mean = (entropy_t + entropy_p) / 2
    return information / mean if mean > 0 else 0.0


def pairwise_precision_recall(truth, predicted):
    """(precision, recall) over mutation pairs; precision is NaN when no pair is placed together."""

    truth, predicted = np.asarray(truth), np.asarray(predicted)
    same_truth = truth[:, None] == truth[None, :]
    same_pred = predicted[:, None] == predicted[None, :]
    upper = np.triu(np.ones_like(same_truth, dtype=bool), 1)
    together = (same_pred & upper).sum()
    correct = (same_pred & same_truth & upper).sum()
    true_pairs = (same_truth & upper).sum()
    precision = correct / together if together else float("nan")
    recall = correct / true_pairs if true_pairs else float("nan")
    return float(precision), float(recall)


# ----------------------------------------------------------------------------
# Selection of mutations to analyze (explicit size limit)
# ----------------------------------------------------------------------------

def select_mutations(eligibility, column, units, max_mutations=MAX_MUTATIONS, pair_budget=PAIR_BUDGET):
    """Eligible mutations, limited (with disclosure) when V^2 x units would exceed the budget."""

    eligible = eligibility[eligibility[column]].copy()
    limit = int(min(max_mutations, math.floor(math.sqrt(pair_budget / max(units, 1)))))
    if len(eligible) <= limit:
        return eligible, ""
    known = "Known_Units_With_Counts" if column == "Eligible_Depth_Aware" else "Known_Units"
    chosen = eligible.sort_values([known, "Reference_ID", "Position", "Alt"],
                                  ascending=[False, True, True, True], kind="mergesort").head(limit)
    note = (f"{len(eligible)} mutations were eligible; to keep the pairwise step within the server's limits "
            f"only the {limit} with the most known observations were grouped. The other "
            f"{len(eligible) - limit} are listed as not analyzed, not as ungrouped.")
    return chosen, note


def assignment_table(keys, info, result, method_label, model, assumed_rho, prefix, known, stability=None):
    """One row per analyzed mutation with its group (or why it is ungrouped)."""

    groups = sorted(result["groups"], key=lambda g: (-len(g), min(g)))
    name_of, size_of, stable_of = {}, {}, {}
    original = [sorted(g) for g in result["groups"]]
    for number, group in enumerate(groups, start=1):
        for index in group:
            name_of[index] = f"{prefix}{number}"
            size_of[index] = len(group)
            if stability is not None:
                stable_of[index] = stability[original.index(sorted(group))]
    distance = result["pairs"]["distance"].copy()
    np.fill_diagonal(distance, np.inf)
    partners = np.isfinite(distance).sum(axis=1)
    rows = []
    for index, key in enumerate(keys):
        nearest = int(np.argmin(distance[index])) if partners[index] else None
        status = GROUPED if index in name_of else (NOT_SIMILAR if partners[index] else NO_PARTNER)
        rows.append({
            **{name: info[key][name] for name in ["Reference_ID", "Position", "Ref", "Alt"]},
            "Variant_Key": key, "Method": method_label, "Model": model,
            "Assumed_Rho": assumed_rho, "Group_ID": name_of.get(index, ""), "Group_Size": size_of.get(index, 1),
            "Assignment_Status": status, "Comparable_Partners": int(partners[index]),
            "Nearest_Partner": keys[nearest] if nearest is not None else "",
            "Nearest_Distance": float(distance[index, nearest]) if nearest is not None else np.nan,
            "Known_Units": int(known[key]), "Group_Stability": stable_of.get(index, np.nan),
        })
    return pd.DataFrame(rows, columns=ASSIGNMENT_COLUMNS), {name_of[g[0]]: g for g in groups}
