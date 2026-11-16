"""Read-sampling uncertainty for allele frequencies measured from read counts.

The simple model used in Sprint 5 Phase 1: at one position, each of the n
usable reads independently carries the alternative allele with probability p,
so the alternative read count x is Binomial(n, p) and the observed variant
allele frequency (VAF) is x / n.

The 95% Wilson score interval (Wilson 1927; Newcombe 1998, method 3) is

    centre    = (p_hat + z^2 / (2n)) / (1 + z^2 / n)
    half      = z / (1 + z^2 / n) * sqrt(p_hat (1 - p_hat) / n + z^2 / (4 n^2))
    interval  = [centre - half, centre + half]

with p_hat = x / n and z = 1.959964 (the 97.5th percentile of the standard
normal distribution). Unlike the simple p_hat +/- z * SE interval it stays
inside [0, 1] and does not collapse to width 0 when x = 0 or x = n.

What the interval does and does not include:

- It describes read-sampling variation only, under the binomial model.
- It does NOT include sequencing or PCR errors, amplicon dropout or primer
  bias, mapping errors, RNA degradation, or the biological variation between
  wastewater samples. Real pooled data are usually more variable
  (overdispersed), so these intervals are a lower bound on uncertainty, not a
  calibrated measure of all wastewater uncertainty.
- It is not the MAGOS error model and is not used to group mutations.
"""

import numpy as np
import pandas as pd


Z_95 = 1.959963984540054  # standard normal 97.5th percentile
DEMONSTRATION_VAF = 0.20
DEMONSTRATION_DEPTHS = [20, 50, 100, 300, 1000]


def wilson_interval(alt_reads, depth, z=Z_95):
    """Wilson score interval for alt_reads / depth; arrays allowed.

    Returns (lower, upper) as float arrays. Entries with a missing value,
    depth < 1, or a count outside 0..depth are NaN: no interval is made from
    data that cannot be a binomial count.
    """

    x = np.asarray(alt_reads, dtype=float)
    n = np.asarray(depth, dtype=float)
    valid = np.isfinite(x) & np.isfinite(n) & (n >= 1) & (x >= 0) & (x <= n) & (x == np.round(x))
    safe_n = np.where(valid, n, 1.0)
    p_hat = np.where(valid, x, 0.0) / safe_n
    z2 = z * z
    denominator = 1.0 + z2 / safe_n
    centre = (p_hat + z2 / (2.0 * safe_n)) / denominator
    half = (z / denominator) * np.sqrt(p_hat * (1.0 - p_hat) / safe_n + z2 / (4.0 * safe_n ** 2))
    lower = np.clip(centre - half, 0.0, 1.0)
    upper = np.clip(centre + half, 0.0, 1.0)
    # The exact boundary values (0 when x = 0, 1 when x = n), free of rounding.
    lower = np.where(p_hat == 0.0, 0.0, lower)
    upper = np.where(p_hat == 1.0, 1.0, upper)
    return np.where(valid, lower, np.nan), np.where(valid, upper, np.nan)


def expected_interval_width(vaf, depth, z=Z_95):
    """Width of the Wilson interval for a true proportion `vaf` at `depth` reads (for theory curves).

    The same formula with p_hat = vaf; the count vaf x depth need not be whole.
    """

    n = np.asarray(depth, dtype=float)
    z2 = z * z
    return 2 * (z / (1 + z2 / n)) * np.sqrt(vaf * (1 - vaf) / n + z2 / (4 * n ** 2))


def depth_dependence_table(vaf=DEMONSTRATION_VAF, depths=None):
    """The same VAF observed at several depths: the interval narrows as depth grows.

    Only depths where vaf * depth is a whole number of reads are used, so
    every row is a real possible observation (for 0.20: 4/20, 10/50, ...).
    """

    rows = []
    for depth in depths or DEMONSTRATION_DEPTHS:
        alt = vaf * depth
        if abs(alt - round(alt)) > 1e-9:
            continue
        lower, upper = wilson_interval([round(alt)], [depth])
        rows.append({"Alt_Reads": int(round(alt)), "Depth": depth, "VAF": round(alt) / depth,
                     "Wilson_Lower": float(lower[0]), "Wilson_Upper": float(upper[0]),
                     "Interval_Width": float(upper[0] - lower[0])})
    return pd.DataFrame(rows, columns=["Alt_Reads", "Depth", "VAF", "Wilson_Lower", "Wilson_Upper",
                                       "Interval_Width"])
