"""Synthetic validation of the Sprint 5 Phase 2 overdispersion-aware intervals.

Everything here is SYNTHETIC: read counts are drawn from beta-binomial models
with known allele frequencies, known overdispersion rho, and fixed seeds. It
tests whether the methods behave as intended WHEN THE MODEL'S ASSUMPTIONS
HOLD. It does not show that real wastewater libraries follow a beta-binomial
model, that one rho fits every site, or that any real sample is calibrated.

The performance criteria were fixed before this study was run (see
docs/SPRINT5_PHASE2_PROGRESS.md):

C1  Known rho, exact beta-binomial sums: for each rho, mean coverage 0.95 +/- 0.01
    and minimum >= 0.90 over depth {20, 50, 100, 300} x p {0.01 ... 0.99}.
C2  The Phase 1 scenario (rho = 0.02, depths 20-300, p 0.05/0.2/0.5): Wilson's
    under-coverage is reproduced, and the new interval's simulated coverage
    is within 4 Monte Carlo SE of its exact coverage.
C3  Estimated rho from m replicate sites: mean coverage >= 0.93 and every
    cell >= 0.90; designs that fail are not supported.
C4  The published plug-in likelihood estimator is evaluated on the same data.

Tables (written to <folder>/overdispersion_*.csv when a folder is given):
known_rho, phase1_scenario, estimated_rho, estimator_bias, sensitivity,
heterogeneous_rho, runtime.

Usage (from the project folder; about 4-6 minutes):
    python validation/simulate_overdispersion.py [output folder]
"""

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genomics import beta_binomial as bb  # noqa: E402
from genomics import overdispersion as od  # noqa: E402
from genomics.frequency_uncertainty import wilson_interval  # noqa: E402
from validation.simulate_vaf_intervals import SEED as PHASE1_SEED  # noqa: E402

SEED = 20261115
DEPTHS = [20, 50, 100, 300]
TRUE_VAFS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95, 0.98, 0.99]
RHOS = [0.0, 0.001, 0.005, 0.02, 0.05]
TEST_VAFS = [0.05, 0.2, 0.5]
CALIBRATION_SITES = [5, 10, 20, 50]
CALIBRATION_RHOS = [0.0, 0.002, 0.01, 0.02, 0.05]
CALIBRATION_DATASETS = 400
CALIBRATION_DEPTHS = (200, 3000)        # replicate depths, all above the 100x site filter
CALIBRATION_VAF_RANGE = (0.02, 0.80)    # site frequencies inside the 1%-90% filter


def draw_counts(generator, depth, p, rho, size=None):
    """Beta-binomial (or binomial when rho = 0) alternative read counts."""

    if rho == 0:
        return generator.binomial(depth, p, size=size)
    alpha, beta = p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho
    return generator.binomial(depth, generator.beta(alpha, beta, size=size))


def intervals_for_all_counts(depth, rho, method):
    k = np.arange(depth + 1)
    if method == "Wilson":
        return wilson_interval(k, np.full(depth + 1, depth))
    if method == "Beta-binomial score":
        return bb.score_interval(k, np.full(depth + 1, depth), rho)
    bounds = np.array([bb.likelihood_ratio_interval(int(value), depth, rho) for value in k])
    return bounds[:, 0], bounds[:, 1]


def exact_coverage_and_width(p, depth, true_rho, lower, upper):
    probabilities = bb.pmf_all(depth, p, true_rho)
    covered = (lower <= p) & (p <= upper)
    return float(probabilities[covered].sum()), float(np.sum(probabilities * (upper - lower)))


def known_rho_table():
    """C1: exact coverage with the true rho supplied to the interval (no Monte Carlo error)."""

    rows = []
    for rho in RHOS:
        for depth in DEPTHS:
            bounds = {method: intervals_for_all_counts(depth, rho, method)
                      for method in ["Wilson", "Beta-binomial score", "Beta-binomial likelihood ratio"]}
            for p in TRUE_VAFS:
                row = {"True_Rho": rho, "Depth": depth, "True_VAF": p}
                for method, (lower, upper) in bounds.items():
                    coverage, width = exact_coverage_and_width(p, depth, rho, lower, upper)
                    row[f"{method} coverage"] = coverage
                    row[f"{method} mean width"] = width
                rows.append(row)
    return pd.DataFrame(rows)


def phase1_scenario_table(replicates=20000):
    """C2: the Phase 1 overdispersion scenario, same seed and draws, with both intervals."""

    generator = np.random.default_rng(PHASE1_SEED + 1)  # the seed used by Phase 1's overdispersion_table
    rho = 0.02
    rows = []
    for depth in DEPTHS:
        for p in TEST_VAFS:
            alpha, beta = p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho
            reads = generator.binomial(depth, generator.beta(alpha, beta, size=replicates))
            depths = np.full(replicates, depth)
            w_lower, w_upper = wilson_interval(reads, depths)
            s_lower, s_upper = bb.score_interval(reads, depths, rho)
            exact_lower, exact_upper = intervals_for_all_counts(depth, rho, "Beta-binomial score")
            exact, _ = exact_coverage_and_width(p, depth, rho, exact_lower, exact_upper)
            simulated = float(np.mean((s_lower <= p) & (p <= s_upper)))
            rows.append({"Depth": depth, "True_VAF": p, "True_Rho": rho, "Replicates": replicates,
                         "Wilson_Coverage": float(np.mean((w_lower <= p) & (p <= w_upper))),
                         "BetaBinomial_Simulated_Coverage": simulated,
                         "BetaBinomial_Exact_Coverage": exact,
                         "Monte_Carlo_SE": float(np.sqrt(exact * (1 - exact) / replicates)),
                         "Wilson_Mean_Width": float(np.mean(w_upper - w_lower)),
                         "BetaBinomial_Mean_Width": float(np.mean(s_upper - s_lower))})
    return pd.DataFrame(rows)


def calibration_dataset(generator, sites, rho, heterogeneous=False):
    p = generator.uniform(*CALIBRATION_VAF_RANGE, sites)
    n1 = generator.integers(*CALIBRATION_DEPTHS, sites)
    n2 = generator.integers(*CALIBRATION_DEPTHS, sites)
    site_rho = generator.gamma(1.0, rho, sites) if heterogeneous and rho > 0 else np.full(sites, rho)
    x1 = np.array([draw_counts(generator, n, q, r) for n, q, r in zip(n1, p, site_rho)])
    x2 = np.array([draw_counts(generator, n, q, r) for n, q, r in zip(n2, p, site_rho)])
    return x1, n1, x2, n2


def _test_cells():
    """Exact pmfs for the test observations, computed once per (depth, p, rho)."""

    return {(depth, p, rho): bb.pmf_all(depth, p, rho)
            for depth in DEPTHS for p in TEST_VAFS for rho in set(CALIBRATION_RHOS)}


def _coverage_given_rho(cells, true_rho, used_rho):
    """Exact coverage and expected width in each test cell for an interval built with `used_rho`."""

    results = {}
    for depth in DEPTHS:
        k = np.arange(depth + 1)
        lower, upper = bb.score_interval(k, np.full(depth + 1, depth), used_rho)
        for p in TEST_VAFS:
            probabilities = cells[(depth, p, true_rho)]
            covered = (lower <= p) & (p <= upper)
            results[(depth, p)] = (float(probabilities[covered].sum()),
                                   float(np.sum(probabilities * (upper - lower))))
    return results


def estimated_rho_table(datasets=CALIBRATION_DATASETS, seed=SEED, sites_list=None, rhos=None,
                        bootstrap=od.BOOTSTRAP_REPLICATES):
    """C3: coverage when rho is estimated from simulated technical replicates."""

    generator = np.random.default_rng(seed)
    cells = _test_cells()
    rows, bias_rows = [], []
    for sites in sites_list or CALIBRATION_SITES:
        for rho in rhos or CALIBRATION_RHOS:
            sums = {rule: {} for rule in ["Point estimate", "Upper 95% bound"]}
            estimates, uppers, unstable = [], [], 0
            for dataset in range(datasets):
                x1, n1, x2, n2 = calibration_dataset(generator, sites, rho)
                estimate = od.moment_estimate(x1, n1, x2, n2)
                boot = od.bootstrap_estimates(x1, n1, x2, n2, bootstrap, seed + dataset)
                upper = float(np.quantile(boot, 0.975))
                if estimate >= od.RHO_SEARCH_MAX or np.mean(boot >= od.RHO_SEARCH_MAX) > od.BOUNDARY_SHARE_LIMIT:
                    unstable += 1
                estimates.append(estimate)
                uppers.append(upper)
                for rule, used in [("Point estimate", estimate), ("Upper 95% bound", upper)]:
                    for cell, (coverage, width) in _coverage_given_rho(cells, rho, used).items():
                        total = sums[rule].setdefault(cell, [0.0, 0.0])
                        total[0] += coverage
                        total[1] += width
            for rule, totals in sums.items():
                coverage = np.array([value[0] / datasets for value in totals.values()])
                width = np.array([value[1] / datasets for value in totals.values()])
                rows.append({"Sites": sites, "True_Rho": rho, "Rule": rule, "Datasets": datasets,
                             "Mean_Coverage": float(coverage.mean()), "Min_Cell_Coverage": float(coverage.min()),
                             "Mean_Width": float(width.mean()),
                             "Supported": bool(coverage.mean() >= 0.93 and coverage.min() >= 0.90)})
            estimates, uppers = np.array(estimates), np.array(uppers)
            bias_rows.append({"Sites": sites, "True_Rho": rho, "Estimator": "Moment (this application)",
                              "Mean_Estimate": float(estimates.mean()), "Bias": float(estimates.mean() - rho),
                              "Relative_Bias": float((estimates.mean() - rho) / rho) if rho else np.nan,
                              "RMSE": float(np.sqrt(np.mean((estimates - rho) ** 2))),
                              "Share_At_Zero": float(np.mean(estimates == 0)),
                              "Upper_Bound_Covers_True_Rho": float(np.mean(uppers >= rho)),
                              "Unstable_Datasets": unstable, "Datasets": datasets})
    return pd.DataFrame(rows), pd.DataFrame(bias_rows)


def plugin_bias_table(datasets=60, seed=SEED + 7, sites=20):
    """C4: the published plug-in likelihood versus the moment estimator on identical data."""

    generator = np.random.default_rng(seed)
    cells = _test_cells()
    rows = []
    for rho in CALIBRATION_RHOS:
        moment, plugin, coverage = [], [], []
        for _ in range(datasets):
            x1, n1, x2, n2 = calibration_dataset(generator, sites, rho)
            moment.append(od.moment_estimate(x1, n1, x2, n2))
            estimate, _ = od.plugin_likelihood_estimate(x1, n1, x2, n2)
            plugin.append(estimate)
            coverage.append(np.mean([value[0] for value in _coverage_given_rho(cells, rho, estimate).values()]))
        rows.append({"Sites": sites, "True_Rho": rho, "Datasets": datasets,
                     "Moment_Mean": float(np.mean(moment)), "Plugin_Likelihood_Mean": float(np.mean(plugin)),
                     "Plugin_Relative_Bias": float((np.mean(plugin) - rho) / rho) if rho else np.nan,
                     "Moment_Relative_Bias": float((np.mean(moment) - rho) / rho) if rho else np.nan,
                     "Coverage_Using_Plugin_Point_Estimate": float(np.mean(coverage))})
    return pd.DataFrame(rows)


def sensitivity_table():
    """Assumed (hypothetical) rho versus the true rho: what a wrong assumption does to coverage."""

    rows = []
    values = [0.0, 0.002, 0.02, 0.05]
    for depth in [50, 300]:
        for true_rho in values:
            probabilities = bb.pmf_all(depth, 0.2, true_rho)
            for assumed in values:
                lower, upper = intervals_for_all_counts(depth, assumed, "Beta-binomial score")
                covered = (lower <= 0.2) & (0.2 <= upper)
                rows.append({"Depth": depth, "True_VAF": 0.2, "True_Rho": true_rho, "Assumed_Rho": assumed,
                             "Exact_Coverage": float(probabilities[covered].sum())})
    return pd.DataFrame(rows)


def heterogeneous_rho_table(datasets=200, seed=SEED + 11, sites=50):
    """Misspecification: rho differs between sites (Gamma, coefficient of variation 1) around its mean."""

    generator = np.random.default_rng(seed)
    rows = []
    for rho in [0.01, 0.02, 0.05]:
        covered = {depth: [] for depth in DEPTHS}
        for dataset in range(datasets):
            x1, n1, x2, n2 = calibration_dataset(generator, sites, rho, heterogeneous=True)
            boot = od.bootstrap_estimates(x1, n1, x2, n2, 400, seed + dataset)
            used = float(np.quantile(boot, 0.975))
            for depth in DEPTHS:
                test_rho = generator.gamma(1.0, rho, 200)
                p = generator.choice(TEST_VAFS, 200)
                reads = np.array([draw_counts(generator, depth, q, r) for q, r in zip(p, test_rho)])
                lower, upper = bb.score_interval(reads, np.full(200, depth), used)
                covered[depth].append(np.mean((lower <= p) & (p <= upper)))
        for depth, values in covered.items():
            rows.append({"Mean_Rho": rho, "Site_Rho_CV": 1.0, "Sites": sites, "Depth": depth,
                         "Coverage": float(np.mean(values)), "Datasets": datasets})
    return pd.DataFrame(rows)


def runtime_table():
    generator = np.random.default_rng(SEED + 13)
    reads = generator.integers(0, 300, 8000)
    depth = np.full(8000, 300)
    start = time.perf_counter()
    bb.score_interval(reads, depth, 0.02)
    score_seconds = time.perf_counter() - start
    x1, n1, x2, n2 = calibration_dataset(generator, 50, 0.02)
    start = time.perf_counter()
    od.moment_estimate(x1, n1, x2, n2)
    od.bootstrap_estimates(x1, n1, x2, n2)
    calibration_seconds = time.perf_counter() - start
    start = time.perf_counter()
    od.plugin_likelihood_estimate(x1, n1, x2, n2)
    plugin_seconds = time.perf_counter() - start
    start = time.perf_counter()
    for k in range(0, 301, 10):
        bb.likelihood_ratio_interval(k, 300, 0.02)
    lr_seconds = (time.perf_counter() - start) / 31
    return pd.DataFrame([
        {"Task": "Score interval, 8,000 observations", "Seconds": score_seconds},
        {"Task": "Moment estimate + 1,000 bootstrap resamples, 50 sites", "Seconds": calibration_seconds},
        {"Task": "Plug-in likelihood (grid + refinement), 50 sites", "Seconds": plugin_seconds},
        {"Task": "Likelihood-ratio interval, one observation at depth 300", "Seconds": lr_seconds},
    ])


def criteria_summary(known, phase1, estimated):
    lines = []
    for rho, group in known.groupby("True_Rho"):
        for method in ["Wilson", "Beta-binomial score", "Beta-binomial likelihood ratio"]:
            values = group[f"{method} coverage"]
            passed = abs(values.mean() - 0.95) <= 0.01 and values.min() >= 0.90
            lines.append(f"C1 rho={rho:<6} {method:<31} mean {values.mean():.4f} min {values.min():.4f} "
                         f"{'PASS' if passed else 'FAIL'}")
    gap = (phase1["BetaBinomial_Simulated_Coverage"] - phase1["BetaBinomial_Exact_Coverage"]).abs()
    lines.append(f"C2 Wilson coverage at 300x: {phase1[phase1['Depth'] == 300]['Wilson_Coverage'].round(3).tolist()}; "
                 f"new interval within 4 SE of exact: {bool((gap <= 4 * phase1['Monte_Carlo_SE']).all())}")
    for (sites, rule), group in estimated.groupby(["Sites", "Rule"]):
        lines.append(f"C3 sites={sites:<3} {rule:<16} supported for every rho: {bool(group['Supported'].all())} "
                     f"(min mean coverage {group['Mean_Coverage'].min():.4f}, "
                     f"min cell {group['Min_Cell_Coverage'].min():.4f})")
    return lines


def main(folder=None):
    started = time.perf_counter()
    known = known_rho_table()
    phase1 = phase1_scenario_table()
    estimated, bias = estimated_rho_table()
    tables = {"overdispersion_known_rho.csv": known,
              "overdispersion_phase1_scenario.csv": phase1,
              "overdispersion_estimated_rho.csv": estimated,
              "overdispersion_estimator_bias.csv": bias,
              "overdispersion_plugin_comparison.csv": plugin_bias_table(),
              "overdispersion_sensitivity.csv": sensitivity_table(),
              "overdispersion_heterogeneous_rho.csv": heterogeneous_rho_table(),
              "overdispersion_runtime.csv": runtime_table()}
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 30)
    for name, table in tables.items():
        if name != "overdispersion_known_rho.csv":
            print(f"\n== {name} (SYNTHETIC; seed {SEED}) ==")
            print(table.round(4).to_string(index=False))
        if folder:
            Path(folder).mkdir(parents=True, exist_ok=True)
            table.to_csv(Path(folder) / name, index=False)
    print("\n== Criteria ==")
    print("\n".join(criteria_summary(known, phase1, estimated)))
    print(f"\nTotal {time.perf_counter() - started:.0f} s")
    return tables


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
