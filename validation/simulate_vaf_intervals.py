"""Synthetic validation of the Sprint 5 depth-aware VAF intervals (Phase 1).

Everything here is SYNTHETIC: allele frequencies and read counts are drawn
from known models with fixed random seeds. It checks the statistical
behaviour of the 95% Wilson interval under the model it assumes (binomial
read sampling); it says nothing about real wastewater samples and does not
validate lineage detection.

Four checks, written to <output folder>/simulation_*.csv and summarized on screen:

1. Bias, coverage, and width at depths 20, 50, 100, 300 for true VAFs from
   0.01 to 0.99, by simulation (seed 20261114) and by exact binomial sums.
2. Boundaries: coverage for VAF near 0 and 1.
3. Overdispersion (beta-binomial reads, a common feature of amplicon data):
   how far coverage falls when the binomial assumption is wrong.
4. Missing coverage in a time series: a mutation whose true VAF rises
   linearly; low-VAF samples are not reported. The slope from detections only
   (the old Sprint 3 behaviour) is compared with the slope using supported
   non-detections from base counts.

Usage (from the project folder):
    python validation/simulate_vaf_intervals.py [output folder]
"""

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genomics.frequency_uncertainty import wilson_interval  # noqa: E402


SEED = 20261114
DEPTHS = [20, 50, 100, 300]
TRUE_VAFS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95, 0.98, 0.99]
REPLICATES = 20000


def exact_coverage(p, n):
    """P(lower <= p <= upper) summed over every possible count x = 0..n (no simulation)."""

    x = np.arange(n + 1)
    lower, upper = wilson_interval(x, np.full(n + 1, n))
    probability = np.array([math.comb(n, k) * p ** k * (1 - p) ** (n - k) for k in x])
    covered = (lower <= p) & (p <= upper)
    return float(np.sum(probability[covered]))


def coverage_table(replicates=REPLICATES, seed=SEED):
    generator = np.random.default_rng(seed)
    rows = []
    for depth in DEPTHS:
        for p in TRUE_VAFS:
            reads = generator.binomial(depth, p, size=replicates)
            lower, upper = wilson_interval(reads, np.full(replicates, depth))
            observed = reads / depth
            rows.append({"Depth": depth, "True_VAF": p,
                         "Mean_Observed_VAF": observed.mean(), "Bias": observed.mean() - p,
                         "Simulated_Coverage": float(np.mean((lower <= p) & (p <= upper))),
                         "Exact_Coverage": exact_coverage(p, depth),
                         "Mean_Width": float(np.mean(upper - lower)),
                         "Share_With_Zero_Reads": float(np.mean(reads == 0))})
    return pd.DataFrame(rows)


def overdispersion_table(replicates=REPLICATES, seed=SEED + 1, rho=0.02):
    """Beta-binomial reads with intra-sample correlation rho; the Wilson interval assumes rho = 0."""

    generator = np.random.default_rng(seed)
    rows = []
    for depth in DEPTHS:
        for p in [0.05, 0.2, 0.5]:
            alpha, beta = p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho
            reads = generator.binomial(depth, generator.beta(alpha, beta, size=replicates))
            lower, upper = wilson_interval(reads, np.full(replicates, depth))
            rows.append({"Depth": depth, "True_VAF": p, "Overdispersion_rho": rho,
                         "Coverage": float(np.mean((lower <= p) & (p <= upper)))})
    return pd.DataFrame(rows)


def missing_coverage_table(replicates=2000, seed=SEED + 2, report_limit=0.03):
    """Slope bias when undetected low-VAF samples are dropped instead of counted as observed values."""

    generator = np.random.default_rng(seed)
    days = np.arange(0, 70, 7, dtype=float)
    rows = []
    for label, truth in [("Rising 0.00 -> 0.30", np.linspace(0.0, 0.3, len(days))),
                         ("Falling 0.20 -> 0.00", np.linspace(0.2, 0.0, len(days)))]:
        true_slope = np.polyfit(days, truth, 1)[0] * 7
        for depth in [50, 300]:
            detections_only, with_non_detections = [], []
            for _ in range(replicates):
                reads = generator.binomial(depth, truth)
                vaf = reads / depth
                reported = vaf >= report_limit  # a caller that reports only frequencies >= 3%
                if reported.sum() >= 3:
                    detections_only.append(np.polyfit(days[reported], vaf[reported], 1)[0] * 7)
                # With A/C/G/T counts every covered sample has a known value (0 when no reads).
                with_non_detections.append(np.polyfit(days, vaf, 1)[0] * 7)
            rows.append({"Scenario": label, "Depth": depth, "True_Slope_Per_Week": true_slope,
                         "Mean_Slope_Detections_Only": float(np.mean(detections_only)),
                         "Mean_Slope_With_Supported_Non_Detections": float(np.mean(with_non_detections)),
                         "Replicates_With_3_Detections": len(detections_only), "Replicates": replicates})
    return pd.DataFrame(rows)


def main(folder=None):
    tables = {"simulation_coverage.csv": coverage_table(),
              "simulation_overdispersion.csv": overdispersion_table(),
              "simulation_missing_coverage.csv": missing_coverage_table()}
    pd.set_option("display.width", 200)
    for name, table in tables.items():
        print(f"\n== {name} (SYNTHETIC; seed {SEED}) ==")
        print(table.round(4).to_string(index=False))
        if folder:
            Path(folder).mkdir(parents=True, exist_ok=True)
            table.to_csv(Path(folder) / name, index=False)
    return tables


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
