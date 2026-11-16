"""SYNTHETIC validation of exploratory mutation-frequency grouping (Sprint 5 Phase 3).

Everything here is simulated with known truth (seed 20261115). The results
show how the methods behave WHEN THE SIMULATION ASSUMPTIONS HOLD. They do not
show that real pooled pathogen samples follow these assumptions, and a
grouping that works here is not a lineage reconstruction.

Simulation model (fixed before any result was viewed; see
docs/SPRINT5_PHASE3_PROGRESS.md, Milestone 5)
------------------------------------------------------------------------
- Observation units: 2 sites x 6 dates = 12 specimens (other sizes in the
  "units" study).
- Sources: 3 latent sources per dataset. At site a and date t, source k has
  weight exp(b0 + b1 t) / (1 + sum_k exp(...)) (logistic shares, with
  b0 ~ N(-1, 1) and b1 ~ N(0, 0.6) per site and source, so trajectories
  differ between sites and sources).
- Planted structure: groups of 10, 6, and 3 mutations (unequal sizes)
  have their source's frequency. 2 "shared" mutations carry the sum of two
  sources (a partly overlapping pattern, truth: ungroupable), and 8
  "background" mutations each follow their own independent logistic
  trajectory (truth: ungroupable).
- Negative control: 29 mutations, each with its own independent trajectory.
- Look-alike control: 14 pairs of unrelated mutations with similar
  frequencies. The pair shares a level L ~ U(0.05, 0.5), and each mutation
  is L + N(0, 0.03) independently in every unit (truth: all ungroupable).
- Depth: nominal D in {20, 50, 100, 300} or "mixed" (each unit draws D from
  that set). Depth of a cell = D x unit factor LogNormal(0, 0.35) x cell
  factor LogNormal(0, 0.25), rounded, at least 1 (unequal coverage).
- Reads: beta-binomial with true rho in {0, 0.02} (rho = 0 is binomial).
- Missing coverage: each cell is unknown with probability pi in {0, 0.3, 0.6}.
- Detection-only input: no A/C/G/T counts, and a caller reports a mutation
  only with >= 2 alternative reads and VAF >= 0.02; every other cell is
  unknown (never 0).
- Technical replicates (calibration study): two libraries per unit, each
  with depth D and the same rho. rho is estimated per unit with the Phase 2
  moment estimator (upper 95% bound) from that unit's informative sites.

Attempt history: attempt 1 (no informative-pattern requirement) failed the
pre-registered criterion P1 in the minimum-overlap study; it is kept here as
"Depth-aware attempt 1" so that failure can be reproduced. Attempt 2 adds the
requirement (genomics/mutation_grouping.PATTERN_RANGE); the ablation
"Ordinary + pattern requirement" separates its effect from depth-awareness.

Truth labels: each planted group is one class. Every shared, background,
negative-control, and look-alike mutation is its own class. Ungrouped and
excluded mutations count as their own singletons in ARI and NMI.
"""

import argparse
import math
import sys
import time
import tracemalloc
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomics import mutation_grouping as mg  # noqa: E402
from genomics.overdispersion import (  # noqa: E402
    CALIBRATION_MIN_DEPTH,
    CALIBRATION_VAF_HIGH,
    CALIBRATION_VAF_LOW,
    USABLE_STATUSES,
    estimate_overdispersion,
)

SEED = 20261115
DEPTHS = [20, 50, 100, 300, "mixed"]
STRUCTURES = ["planted", "negative", "lookalike"]
ASSUMED_RHOS = [0.0, 0.002, 0.02, 0.05]
REPORT_MIN_READS, REPORT_MIN_VAF = 2, 0.02


# ----------------------------------------------------------------------------
# Data generation
# ----------------------------------------------------------------------------

def logistic_shares(rng, sources, sites, dates):
    """Source weights per unit: array (sources, units)."""

    t = np.arange(dates) - (dates - 1) / 2
    shares = []
    for _ in range(sites):
        b0 = rng.normal(-1, 1, (sources, 1))
        b1 = rng.normal(0, 0.6, (sources, 1))
        e = np.exp(b0 + b1 * t[None, :])
        shares.append(e / (1 + e.sum(axis=0, keepdims=True)))
    return np.concatenate(shares, axis=1)


def own_trajectories(rng, count, sites, dates):
    """Independent logistic frequency trajectories, one per mutation: (count, units)."""

    t = np.arange(dates) - (dates - 1) / 2
    rows = []
    for _ in range(sites):
        b0 = rng.normal(-1.5, 1, (count, 1))
        b1 = rng.normal(0, 0.6, (count, 1))
        rows.append(1 / (1 + np.exp(-(b0 + b1 * t[None, :]))))
    return np.concatenate(rows, axis=1)


def true_frequencies(rng, structure, sites=2, dates=6, sizes=(10, 6, 3), shared=2, background=8):
    """Return (frequency array (V, U), truth labels)."""

    units = sites * dates
    if structure == "planted":
        shares = logistic_shares(rng, len(sizes), sites, dates)
        rows, labels = [], []
        for k, size in enumerate(sizes):
            rows += [shares[k]] * size
            labels += [k] * size
        for s in range(shared):
            a, b = s % len(sizes), (s + 1) % len(sizes)
            rows.append(np.clip(shares[a] + shares[b], 0, 0.999))
            labels.append(100 + s)
        rows += list(own_trajectories(rng, background, sites, dates))
        labels += list(range(200, 200 + background))
        return np.array(rows), np.array(labels)
    total = sum(sizes) + shared + background
    if structure == "negative":
        return own_trajectories(rng, total, sites, dates), np.arange(total) + 300
    if structure == "lookalike":
        pairs = total // 2
        level = rng.uniform(0.05, 0.5, pairs)
        rows = np.clip(np.repeat(level, 2)[:, None] + rng.normal(0, 0.03, (2 * pairs, units)), 0.005, 0.995)
        return rows, np.arange(2 * pairs) + 400
    raise ValueError(structure)


def draw_depths(rng, depth, shape):
    count, units = shape
    if depth == "mixed":
        nominal = rng.choice([20, 50, 100, 300], units).astype(float)
    else:
        nominal = np.full(units, float(depth))
    factor = nominal * rng.lognormal(0, 0.35, units)
    return np.maximum(1, np.round(factor[None, :] * rng.lognormal(0, 0.25, shape))).astype(float)


def beta_binomial(rng, n, p, rho):
    if rho <= 0:
        return rng.binomial(n.astype(int), p).astype(float)
    a = p * (1 - rho) / rho
    b = (1 - p) * (1 - rho) / rho
    return rng.binomial(n.astype(int), rng.beta(a, b)).astype(float)


def simulate(rng, structure, depth, rho, missing, detection_only=False, replicates=False, sites=2, dates=6,
             sizes=(10, 6, 3), shared=2, background=8):
    """One synthetic dataset: arrays for grouping, truth labels, and per-library counts when replicated."""

    p, truth = true_frequencies(rng, structure, sites, dates, sizes, shared, background)
    libraries = 2 if replicates else 1
    xs, ns = [], []
    for _ in range(libraries):
        n = draw_depths(rng, depth, p.shape)
        xs.append(beta_binomial(rng, n, p, rho))
        ns.append(n)
    x, n = sum(xs), sum(ns)
    q = sum(value ** 2 for value in ns)
    known = rng.random(p.shape) >= missing
    if detection_only:
        known &= (x >= REPORT_MIN_READS) & (x / n >= REPORT_MIN_VAF)
    arrays = {"F": np.where(known, x / n, np.nan), "X": np.where(known, x, np.nan),
              "N": np.where(known, n, np.nan), "Q": np.where(known, q, np.nan)}
    return {"arrays": arrays, "truth": truth, "true_p": p,
            "libraries": [(np.where(known, a, np.nan), np.where(known, b, np.nan)) for a, b in zip(xs, ns)]}


def calibrated_rhos(data, bootstrap=1000):
    """Per-unit rho (upper 95% bound) from the two simulated libraries; NaN when not estimable."""

    (x1, n1), (x2, n2) = data["libraries"]
    rhos = []
    for s in range(x1.shape[1]):
        a, b, c, d = x1[:, s], n1[:, s], x2[:, s], n2[:, s]
        mean = (a / b + c / d) / 2
        use = (np.isfinite(a) & (b > CALIBRATION_MIN_DEPTH) & (d > CALIBRATION_MIN_DEPTH)
               & (mean > CALIBRATION_VAF_LOW) & (mean < CALIBRATION_VAF_HIGH))
        sites = pd.DataFrame({"Alt_Reads_1": a[use], "Depth_1": b[use], "Alt_Reads_2": c[use], "Depth_2": d[use]})
        result = estimate_overdispersion(sites, replicates=bootstrap, seed=SEED + s, include_plugin=False)
        rhos.append(result["Rho_Used_For_Intervals"] if result["Status"] in USABLE_STATUSES else np.nan)
    return np.array(rhos, dtype=float)


# ----------------------------------------------------------------------------
# Evaluation
# ----------------------------------------------------------------------------

def score(truth, groups, eligible_count, runtime):
    labels = mg.label_vector(groups, len(truth))
    precision, recall = mg.pairwise_precision_recall(truth, labels)
    singleton_truth = pd.Series(truth).map(pd.Series(truth).value_counts()).eq(1).to_numpy()
    grouped = np.zeros(len(truth), bool)
    for group in groups:
        grouped[list(group)] = True
    same_truth = truth[:, None] == truth[None, :]
    together = labels[:, None] == labels[None, :]
    upper = np.triu(np.ones_like(same_truth), 1)
    negatives = (~same_truth & upper).sum()
    return {
        "ARI": mg.adjusted_rand_index(truth, labels),
        "NMI": mg.normalized_mutual_information(truth, labels),
        "Precision": precision, "Recall": recall,
        "False_Grouping_Share": (float(grouped[singleton_truth].mean()) if singleton_truth.any() else np.nan),
        "Pair_False_Positive_Rate": float((together & ~same_truth & upper).sum() / negatives) if negatives else np.nan,
        "Groups": len(groups), "Grouped_Mutations": int(grouped.sum()), "Eligible": eligible_count,
        "Excluded": len(truth) - eligible_count, "Abstained": len(groups) == 0, "Runtime_s": runtime,
    }


def run_methods(data, minimum_shared, methods, unit_rhos=None, stability=False):
    """Run each configured method; return one metric dict per method."""

    arrays, truth = data["arrays"], data["truth"]
    count, units = arrays["F"].shape
    references = ["SIM"] * count
    known_units = np.isfinite(arrays["F"]).sum(axis=1)
    eligible = known_units >= minimum_shared
    index = np.flatnonzero(eligible)
    reduced = {name: values[index] for name, values in arrays.items()}
    rows = []
    for name, kind, rho in methods:
        if kind == "none":
            rows.append({"Method": name, **score(truth, [], len(index), 0.0)})
            continue
        if kind == "calibrated":
            if unit_rhos is None or not np.isfinite(unit_rhos).all():
                rows.append({"Method": name, "Calibration_Available": False,
                             **score(truth, [], len(index), 0.0)})
                continue
            rho_vector = unit_rhos
        else:
            rho_vector = np.full(units, rho)
        start = time.perf_counter()
        method_kind, extra = (("depth", {"minimum_range": None}) if kind == "depth_attempt1"
                              else (kind, {}))
        result = (mg.run_method(method_kind, reduced, references[:len(index)], rho_vector, minimum_shared, **extra)
                  if len(index) >= 2 else {"groups": []})
        runtime = time.perf_counter() - start
        groups = [list(index[g]) for g in result["groups"]]
        row = {"Method": name, **score(truth, groups, len(index), runtime)}
        if kind == "calibrated":
            row["Calibration_Available"] = True
        if stability and result["groups"]:
            values = mg.leave_one_unit_out_stability(method_kind, reduced, references[:len(index)], rho_vector,
                                                     minimum_shared, result["groups"], **extra)
            row["Stability"] = float(np.mean(values))
        rows.append(row)
    return rows


def ordinary_with_height(height):
    """An ordinary-clustering variant with another cut height (baseline sensitivity)."""

    def run(arrays, references, minimum_shared):
        pairs = mg.ordinary_distances(arrays["F"], references, minimum_shared)
        return mg.ordinary_groups(pairs["distance"], height)[0]
    return run


ABLATION = ("Ordinary + pattern requirement (ablation)", "ordinary_pattern", 0.0)
ATTEMPT1 = ("Depth-aware attempt 1 (no pattern requirement)", "depth_attempt1", 0.0)


def main_methods():
    methods = [("Ordinary (h=0.10)", "ordinary", 0.0), ABLATION]
    methods += [(f"Depth-aware (rho={rho:g})", "depth", rho) for rho in ASSUMED_RHOS]
    methods.append(("No grouping", "none", 0.0))
    return methods


def summarise(table, by):
    """Mean and 95% normal interval across seeds for every metric column."""

    metrics = ["ARI", "NMI", "Precision", "Recall", "False_Grouping_Share", "Pair_False_Positive_Rate",
               "Groups", "Grouped_Mutations", "Excluded", "Abstained", "Runtime_s", "Stability"]
    metrics = [name for name in metrics if name in table]
    rows = []
    for keys, group in table.groupby(by, sort=False, dropna=False):
        row = dict(zip(by, keys if isinstance(keys, tuple) else (keys,)))
        row["Datasets"] = len(group)
        for name in metrics:
            values = group[name].astype(float).dropna()
            row[name] = float(values.mean()) if len(values) else np.nan
            row[name + "_CI95"] = (1.96 * float(values.std(ddof=1)) / math.sqrt(len(values))
                                   if len(values) > 1 else np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# Studies
# ----------------------------------------------------------------------------

def minimum_shared_study(seeds=30, candidates=(2, 3, 4, 6)):
    """Step 1: pick the minimum number of shared units (pre-registered rule, see the progress file)."""

    rows = []
    for structure in STRUCTURES:
        for depth in [50, 300]:
            for missing in [0.3, 0.6]:
                for seed in range(seeds):
                    rng = np.random.default_rng([SEED, 1, seed, STRUCTURES.index(structure), depth, int(missing * 10)])
                    data = simulate(rng, structure, depth, 0.0, missing)
                    for m in candidates:
                        for row in run_methods(data, m, [("Ordinary (h=0.10)", "ordinary", 0.0), ABLATION,
                                                         ("Depth-aware (rho=0)", "depth", 0.0), ATTEMPT1]):
                            rows.append({"Structure": structure, "Depth": depth, "Missing": missing,
                                         "Minimum_Shared": m, "Seed": seed, **row})
    return pd.DataFrame(rows)


def choose_minimum_shared(summary, candidates=(2, 3, 4, 6), limit=0.10):
    """Smallest candidate whose depth-aware method meets P1 in every non-abstaining control cell."""

    depth_aware = summary[summary["Method"].eq("Depth-aware (rho=0)")
                          & summary["Structure"].isin(["negative", "lookalike"])]
    for m in candidates:
        cells = depth_aware[depth_aware["Minimum_Shared"].eq(m)]
        if (cells["False_Grouping_Share"].fillna(0) <= limit).all():
            return m
    return None  # no candidate meets P1: reported as a failure, never silently replaced


def main_study(minimum_shared, seeds=30, detection_only=False):
    """Main grid: structure x depth x true rho x missingness (or the detection-only grid)."""

    rows = []
    missing_levels = [0.0] if detection_only else [0.0, 0.3, 0.6]
    for structure in STRUCTURES:
        for depth in DEPTHS:
            for rho in [0.0, 0.02]:
                for missing in missing_levels:
                    for seed in range(seeds):
                        rng = np.random.default_rng([SEED, 2 + detection_only, seed, STRUCTURES.index(structure),
                                                     DEPTHS.index(depth), int(rho * 1000), int(missing * 10)])
                        data = simulate(rng, structure, depth, rho, missing, detection_only=detection_only)
                        stability = (structure == "planted" and depth in (100, 300) and missing <= 0.3
                                     and seed < 10)
                        for row in run_methods(data, minimum_shared, main_methods(), stability=stability):
                            rows.append({"Structure": structure, "Depth": str(depth), "True_Rho": rho,
                                         "Missing": missing, "Detection_Only": detection_only, "Seed": seed,
                                         **row})
                        for height in [0.05, 0.15]:
                            start = time.perf_counter()
                            arrays = data["arrays"]
                            eligible = np.flatnonzero(np.isfinite(arrays["F"]).sum(axis=1) >= minimum_shared)
                            reduced = {k: v[eligible] for k, v in arrays.items()}
                            groups = ordinary_with_height(height)(reduced, ["SIM"] * len(eligible), minimum_shared)
                            groups = [list(eligible[g]) for g in groups]
                            rows.append({"Structure": structure, "Depth": str(depth), "True_Rho": rho,
                                         "Missing": missing, "Detection_Only": detection_only, "Seed": seed,
                                         "Method": f"Ordinary (h={height:.2f})",
                                         **score(data["truth"], groups, len(eligible), time.perf_counter() - start)})
    return pd.DataFrame(rows)


def calibration_study(minimum_shared, seeds=30, bootstrap=1000):
    """Two libraries per unit: binomial vs true rho vs rho calibrated from the replicates."""

    rows = []
    for structure in ["planted", "negative", "lookalike"]:
        for depth in [100, 300, 1000]:
            for rho in [0.005, 0.02]:
                for missing in [0.0, 0.3]:
                    for seed in range(seeds):
                        rng = np.random.default_rng([SEED, 4, seed, STRUCTURES.index(structure), depth,
                                                     int(rho * 1000), int(missing * 10)])
                        data = simulate(rng, structure, depth, rho, missing, replicates=True)
                        rhos = calibrated_rhos(data, bootstrap)
                        methods = [("Ordinary (h=0.10)", "ordinary", 0.0), ("Depth-aware (rho=0)", "depth", 0.0),
                                   ("Depth-aware (true rho)", "depth", rho),
                                   ("Depth-aware (calibrated rho)", "calibrated", 0.0),
                                   ("No grouping", "none", 0.0)]
                        for row in run_methods(data, minimum_shared, methods, unit_rhos=rhos):
                            rows.append({"Structure": structure, "Depth_Per_Library": depth, "True_Rho": rho,
                                         "Missing": missing, "Seed": seed,
                                         "Calibrated_Units": int(np.isfinite(rhos).sum()),
                                         "Median_Calibrated_Rho": float(np.nanmedian(rhos)) if np.isfinite(rhos).any()
                                         else np.nan, **row})
    return pd.DataFrame(rows)


def units_study(minimum_shared, seeds=30):
    """Effect of the number of observation units and of more mutations (depth 100, rho 0, 30% missing)."""

    rows = []
    for structure in ["planted", "negative"]:
        for sites, dates in [(1, 3), (1, 6), (2, 6), (4, 6)]:
            for scale in [1, 3]:
                for seed in range(seeds):
                    rng = np.random.default_rng([SEED, 5, seed, STRUCTURES.index(structure), sites, dates, scale])
                    data = simulate(rng, structure, 100, 0.0, 0.3, sites=sites, dates=dates,
                                    sizes=tuple(s * scale for s in (10, 6, 3)), shared=2 * scale,
                                    background=8 * scale)
                    for row in run_methods(data, minimum_shared, [("Ordinary (h=0.10)", "ordinary", 0.0),
                                                                  ("Depth-aware (rho=0)", "depth", 0.0)]):
                        rows.append({"Structure": structure, "Units": sites * dates, "Mutations": len(data["truth"]),
                                     "Seed": seed, **row})
    return pd.DataFrame(rows)


def null_calibration_study(seeds=2000):
    """Type I error of the pairwise test for two mutations with the SAME true frequency (rho known)."""

    rows = []
    for depth in [20, 50, 100, 300]:
        for p in [0.02, 0.05, 0.2, 0.5]:
            for rho in [0.0, 0.02]:
                rng = np.random.default_rng([SEED, 6, depth, int(p * 100), int(rho * 1000)])
                rejected = tested = 0
                for _ in range(seeds):
                    n = draw_depths(rng, depth, (2, 6))
                    x = beta_binomial(rng, n, np.full((2, 6), p), rho)
                    arrays = {"F": x / n, "X": x, "N": n, "Q": n ** 2}
                    statistic, df, p_value = mg.homogeneity_test([0, 1], arrays, np.full(6, rho))
                    if df > 0:
                        tested += 1
                        rejected += p_value < mg.ALPHA
                rows.append({"Depth": depth, "True_VAF": p, "Rho": rho, "Tests": tested,
                             "Type_I_Error": rejected / tested if tested else np.nan, "Nominal": mg.ALPHA})
    return pd.DataFrame(rows)


def runtime_study():
    """Runtime and peak traced memory of the pairwise step and grouping at larger sizes."""

    rows = []
    for mutations, units in [(100, 12), (300, 52), (1000, 52), (1000, 192)]:
        rng = np.random.default_rng([SEED, 7, mutations, units])
        p = rng.uniform(0.02, 0.9, (mutations, units))
        n = np.round(rng.lognormal(math.log(200), 0.5, (mutations, units)))
        x = rng.binomial(n.astype(int), p).astype(float)
        known = rng.random((mutations, units)) > 0.3
        arrays = {"F": np.where(known, x / n, np.nan), "X": np.where(known, x, np.nan),
                  "N": np.where(known, n, np.nan), "Q": np.where(known, n ** 2, np.nan)}
        for kind in ["ordinary", "depth"]:
            tracemalloc.start()
            start = time.perf_counter()
            mg.run_method(kind, arrays, ["SIM"] * mutations, np.zeros(units), 4)
            seconds = time.perf_counter() - start
            peak = tracemalloc.get_traced_memory()[1]
            tracemalloc.stop()
            rows.append({"Mutations": mutations, "Units": units, "Method": kind, "Seconds": seconds,
                         "Peak_Traced_MB": peak / 2 ** 20})
    return pd.DataFrame(rows)


def main(arguments=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output", nargs="?", default="validation_output", type=Path)
    parser.add_argument("--seeds", type=int, default=30)
    args = parser.parse_args(arguments)
    args.output.mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()
    m_table = minimum_shared_study(args.seeds)
    m_summary = summarise(m_table, ["Structure", "Depth", "Missing", "Minimum_Shared", "Method"])
    m_summary.to_csv(args.output / "grouping_minimum_shared.csv", index=False)
    chosen = choose_minimum_shared(m_summary)
    if chosen is None:
        print("No candidate minimum meets P1 for the depth-aware method; the main grid uses the largest "
              "candidate and the failure is reported.")
        chosen = 6
    print(f"Minimum shared units used: {chosen}")

    outputs = {}
    main_table = main_study(chosen, args.seeds)
    detection_table = main_study(chosen, args.seeds, detection_only=True)
    outputs["grouping_main"] = summarise(main_table, ["Structure", "Depth", "True_Rho", "Missing", "Method"])
    outputs["grouping_detection_only"] = summarise(detection_table, ["Structure", "Depth", "True_Rho", "Method"])
    outputs["grouping_calibration"] = summarise(calibration_study(chosen, args.seeds),
                                                ["Structure", "Depth_Per_Library", "True_Rho", "Missing", "Method"])
    outputs["grouping_units"] = summarise(units_study(chosen, args.seeds), ["Structure", "Units", "Mutations", "Method"])
    outputs["grouping_null_calibration"] = null_calibration_study()
    outputs["grouping_runtime"] = runtime_study()
    combined = pd.concat([main_table, detection_table], ignore_index=True)
    outputs["grouping_paired_ari"] = pd.concat(
        [paired_differences(combined), paired_differences(combined, ABLATION[0])], ignore_index=True)
    for name, table in outputs.items():
        table.to_csv(args.output / f"{name}.csv", index=False)
    main_table.to_csv(args.output / "grouping_main_per_dataset.csv", index=False)
    print(f"Wrote {len(outputs) + 2} tables to {args.output} in {time.perf_counter() - started:.0f} s")
    return outputs


def paired_differences(table, baseline="Ordinary (h=0.10)"):
    """Per condition: depth-aware ARI minus ordinary ARI on the same datasets (P3)."""

    keys = ["Structure", "Depth", "True_Rho", "Missing", "Detection_Only", "Seed"]
    base = table[table["Method"].eq(baseline)].set_index(keys)["ARI"]
    rows = []
    for method in [m for m in table["Method"].unique() if m.startswith("Depth-aware") and m != baseline]:
        other = table[table["Method"].eq(method)].set_index(keys)["ARI"]
        difference = (other - base).dropna().reset_index()
        for condition, group in difference.groupby(keys[:-1], sort=False):
            values = group["ARI"].to_numpy()
            half = 1.96 * values.std(ddof=1) / math.sqrt(len(values)) if len(values) > 1 else np.nan
            mean = values.mean()
            verdict = "Improves" if mean - half > 0 else ("Worse" if mean + half < 0 else "Similar")
            rows.append({**dict(zip(keys[:-1], condition)), "Method": method, "Baseline": baseline,
                         "Mean_ARI_Difference": mean,
                         "CI95_Low": mean - half, "CI95_High": mean + half, "Verdict": verdict})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
