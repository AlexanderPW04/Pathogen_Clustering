"""Beta-binomial model for allele read counts (Sprint 5, Phase 2).

Parameterization (written out because several are in use)
---------------------------------------------------------
At one position in one sequencing library:

    n   = total usable depth (reads covering the position),
    k   = reads carrying the alternative allele,
    p   = underlying mean alternative allele frequency,
    rho = overdispersion (intra-library correlation), 0 <= rho < 1.

Binomial model (rho = 0): every read independently carries the allele with
probability p, so k ~ Binomial(n, p) and Var(k) = n p (1 - p).

Beta-binomial model (rho > 0): the allele probability of the library itself
varies around p, q ~ Beta(alpha, beta), and k | q ~ Binomial(n, q), with

    alpha = p (1 - rho) / rho,   beta = (1 - p) (1 - rho) / rho,
    so    alpha + beta = (1 - rho) / rho  and  rho = 1 / (alpha + beta + 1).

Then E[k] = n p (unchanged) and

    Var(k)   = n p (1 - p) [1 + (n - 1) rho]
    Var(k/n) = p (1 - p) / n x [1 + (n - 1) rho].

The factor 1 + (n - 1) rho is Kish's design effect for clustered binomial
data. As n grows, Var(k/n) approaches p (1 - p) rho instead of 0, so deeper
sequencing cannot remove variation that arises before the reads are counted
(for example uneven amplification of the input RNA molecules). This is the
same rho as in Tonkin-Hill et al. (2021, eLife 10:e66857) and the Phase 1
simulation (validation/simulate_vaf_intervals.py).

Numerically stable probabilities
--------------------------------
The probability mass function

    P(k) = C(n, k) B(k + alpha, n - k + beta) / B(alpha, beta)

is evaluated as a product of rising factorials in which 1/rho cancels:

    log P(k) = log C(n, k)
               + sum_{i<k}     log(p (1 - rho) + i rho)
               + sum_{i<n-k}   log((1 - p)(1 - rho) + i rho)
               - sum_{i<n}     log((1 - rho) + i rho).

Every term is a logarithm of a number between 0 and n, so there is no
overflow at large depth, no cancellation of huge gamma functions when rho is
tiny, and rho = 0 gives exactly the binomial k log p + (n - k) log(1 - p).
p = 0 or 1 gives the correct point masses. The cost is proportional to n,
which is fine for the depths in processed pathogen tables.

What the model does not do
--------------------------
The beta-binomial describes extra variation of a frequency measured from one
library around its specimen's mean. It does not correct mapping artifacts,
primer bias that is the same in every library, reference mismatch,
contamination, or real biological differences between samples or dates.
"""

import math

import numpy as np


Z_95 = 1.959963984540054  # standard normal 97.5th percentile (same as Phase 1)
CHI2_1_95 = Z_95 ** 2      # 95th percentile of chi-square with 1 degree of freedom

BINOMIAL = "Binomial"
BETA_BINOMIAL = "Beta-binomial"
NO_MODEL = "No calibrated uncertainty model"

MODEL_TEXT = {
    BINOMIAL: "Independent reads with one common allele probability (read sampling only).",
    BETA_BINOMIAL: ("Allows extra variation beyond independent read sampling, set by the overdispersion "
                    "rho; it does not correct mapping, primer, reference, or contamination errors."),
    NO_MODEL: "The data do not justify a reliable overdispersion-aware uncertainty calculation.",
}


def check_rho(rho):
    """Return rho as a float, or raise ValueError when it is not a valid overdispersion value."""

    try:
        value = float(rho)
    except (TypeError, ValueError):
        raise ValueError(f"Overdispersion rho must be a number; got {rho!r}.") from None
    if not math.isfinite(value) or value < 0 or value >= 1:
        raise ValueError(f"Overdispersion rho must be at least 0 and below 1; got {rho!r}.")
    return value


def shape_parameters(p, rho):
    """(alpha, beta) of the Beta distribution for mean p and overdispersion rho > 0."""

    rho = check_rho(rho)
    if rho == 0:
        raise ValueError("rho = 0 is the binomial limit; alpha and beta are infinite.")
    scale = (1 - rho) / rho
    return p * scale, (1 - p) * scale


def rho_from_shapes(alpha, beta):
    return 1.0 / (alpha + beta + 1.0)


def design_effect(depth, rho):
    """Variance inflation 1 + (n - 1) rho of the beta-binomial relative to the binomial."""

    return 1.0 + (np.asarray(depth, dtype=float) - 1.0) * check_rho(rho)


def count_mean(depth, p):
    return np.asarray(depth, dtype=float) * np.asarray(p, dtype=float)


def count_variance(depth, p, rho):
    n, p = np.asarray(depth, dtype=float), np.asarray(p, dtype=float)
    return n * p * (1 - p) * design_effect(n, rho)


def frequency_variance(depth, p, rho):
    """Var(k / n) = p (1 - p) / n x [1 + (n - 1) rho]."""

    n, p = np.asarray(depth, dtype=float), np.asarray(p, dtype=float)
    return p * (1 - p) / n * design_effect(n, rho)


def _check_count(k, n):
    if not (isinstance(n, (int, np.integer)) and isinstance(k, (int, np.integer))):
        if not (float(n).is_integer() and float(k).is_integer()):
            raise ValueError("Read counts must be whole numbers.")
    k, n = int(k), int(n)
    if n < 1 or k < 0 or k > n:
        raise ValueError(f"Need 0 <= alternative reads <= depth and depth >= 1; got {k} of {n}.")
    return k, n


def _rising_log_sums(count, start, step):
    """sum_{i<count} log(start + i step) for i = 0..count-1 (0 when count = 0)."""

    if count == 0:
        return 0.0
    terms = start + step * np.arange(count, dtype=float)
    with np.errstate(divide="ignore"):
        return float(np.sum(np.log(terms)))


def _log_choose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def log_pmf(k, n, p, rho):
    """log P(K = k) for one observation; -inf when the probability is 0."""

    k, n = _check_count(k, n)
    rho = check_rho(rho)
    p = float(p)
    if not 0 <= p <= 1:
        raise ValueError(f"Allele frequency p must be between 0 and 1; got {p}.")
    log_value = (_log_choose(n, k)
                 + _rising_log_sums(k, p * (1 - rho), rho)
                 + _rising_log_sums(n - k, (1 - p) * (1 - rho), rho)
                 - _rising_log_sums(n, 1 - rho, rho))
    return log_value


def log_pmf_all(n, p, rho):
    """log P(K = k) for every k = 0..n at once (an array of length n + 1)."""

    _, n = _check_count(0, n)
    rho = check_rho(rho)
    i = np.arange(n, dtype=float)
    with np.errstate(divide="ignore"):
        up = np.concatenate([[0.0], np.cumsum(np.log(p * (1 - rho) + i * rho))])
        down = np.concatenate([[0.0], np.cumsum(np.log((1 - p) * (1 - rho) + i * rho))])
    total = float(np.sum(np.log((1 - rho) + i * rho)))
    k = np.arange(n + 1)
    choose = np.array([_log_choose(n, value) for value in k])
    return choose + up[k] + down[n - k] - total


def pmf_all(n, p, rho):
    return np.exp(log_pmf_all(n, p, rho))


def count_prediction_interval(p, n, rho, level=0.95):
    """Central prediction range (lowest, highest) for a FUTURE alternative count at depth n.

    This answers "where would the count of one new library fall if its mean
    frequency is p?". It is not a confidence interval for p.
    """

    probabilities = pmf_all(n, p, rho)
    cumulative = np.cumsum(probabilities)
    tail = (1 - level) / 2
    lowest = int(np.searchsorted(cumulative, tail - 1e-12))
    highest = int(np.searchsorted(cumulative, 1 - tail - 1e-12))
    return lowest, min(highest, n)


def score_interval(alt_reads, depth, rho, z=Z_95):
    """Overdispersion-aware score interval for the mean frequency p; arrays allowed.

    It solves (x/n - p)^2 = z^2 Var(x/n) with the beta-binomial variance
    p (1 - p) / n x [1 + (n - 1) rho], which is the Wilson score interval with
    the depth replaced by the effective depth n / [1 + (n - 1) rho]. With
    rho = 0 it is exactly the Phase 1 Wilson interval. It is an approximation
    (a quasi-likelihood score interval), not an exact beta-binomial interval;
    its coverage is checked by exact beta-binomial sums in
    validation/simulate_overdispersion.py.

    Returns (lower, upper); NaN where the count, depth, or rho is invalid.
    """

    x = np.asarray(alt_reads, dtype=float)
    n = np.asarray(depth, dtype=float)
    r = np.asarray(rho, dtype=float)
    valid = (np.isfinite(x) & np.isfinite(n) & np.isfinite(r) & (n >= 1) & (x >= 0) & (x <= n)
             & (x == np.round(x)) & (r >= 0) & (r < 1))
    safe_n = np.where(valid, n, 1.0)
    p_hat = np.where(valid, x, 0.0) / safe_n
    effective = safe_n / (1.0 + (safe_n - 1.0) * np.where(valid, r, 0.0))
    z2 = z * z
    denominator = 1.0 + z2 / effective
    centre = (p_hat + z2 / (2.0 * effective)) / denominator
    half = (z / denominator) * np.sqrt(p_hat * (1.0 - p_hat) / effective + z2 / (4.0 * effective ** 2))
    lower = np.clip(centre - half, 0.0, 1.0)
    upper = np.clip(centre + half, 0.0, 1.0)
    lower = np.where(p_hat == 0.0, 0.0, lower)
    upper = np.where(p_hat == 1.0, 1.0, upper)
    return np.where(valid, lower, np.nan), np.where(valid, upper, np.nan)


def expected_score_width(vaf, depth, rho, z=Z_95):
    """Width of the score interval for a true frequency `vaf` at `depth` (theory curves)."""

    n = np.asarray(depth, dtype=float)
    effective = n / design_effect(n, rho)
    z2 = z * z
    return 2 * (z / (1 + z2 / effective)) * np.sqrt(vaf * (1 - vaf) / effective + z2 / (4 * effective ** 2))


def _log_likelihood_curve(k, n, rho):
    """A function p -> log P(k | n, p, rho), with the p-free parts computed once."""

    constant = _log_choose(n, k) - _rising_log_sums(n, 1 - rho, rho)
    i_up = np.arange(k, dtype=float)
    i_down = np.arange(n - k, dtype=float)

    def curve(p):
        with np.errstate(divide="ignore"):
            return (constant + float(np.sum(np.log(p * (1 - rho) + i_up * rho)))
                    + float(np.sum(np.log((1 - p) * (1 - rho) + i_down * rho))))
    return curve


def _maximize(function, low, high, iterations=80):
    """Golden-section search for the maximum of a unimodal function on [low, high]."""

    ratio = (math.sqrt(5) - 1) / 2
    a, b = low, high
    c, d = b - ratio * (b - a), a + ratio * (b - a)
    fc, fd = function(c), function(d)
    for _ in range(iterations):
        if fc < fd:
            a, c, fc = c, d, fd
            d = a + ratio * (b - a)
            fd = function(d)
        else:
            b, d, fd = d, c, fc
            c = b - ratio * (b - a)
            fc = function(c)
    return (a + b) / 2


def _crossing(function, low, high, target, iterations=60):
    """Bisection for function(p) = target where function(low) >= target > function(high) or reverse."""

    f_low = function(low) - target
    for _ in range(iterations):
        middle = (low + high) / 2
        f_middle = function(middle) - target
        if (f_middle >= 0) == (f_low >= 0):
            low, f_low = middle, f_middle
        else:
            high = middle
    return (low + high) / 2


def likelihood_ratio_interval(k, n, rho, critical=CHI2_1_95):
    """Likelihood-ratio interval for p from one count with KNOWN rho (exact beta-binomial likelihood).

    The set of p with 2 [max log L - log L(p)] <= 3.84. Used in the
    validation study to compare with the score interval; slower (iterative).
    """

    k, n = _check_count(k, n)
    rho = check_rho(rho)
    curve = _log_likelihood_curve(k, n, rho)
    if k == 0:
        best = 0.0
    elif k == n:
        best = 1.0
    else:
        best = _maximize(curve, 1e-12, 1 - 1e-12)
    target = curve(best) - critical / 2
    lower = 0.0 if k == 0 or curve(0.0) >= target else _crossing(curve, 0.0, best, target)
    upper = 1.0 if k == n or curve(1.0) >= target else _crossing(curve, best, 1.0, target)
    return lower, upper
