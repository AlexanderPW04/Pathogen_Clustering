"""Chi-square tail probabilities without SciPy (Sprint 5 Phase 3).

SciPy is not a project dependency, so the few distribution functions the
mutation grouping needs are implemented here from the regularized incomplete
gamma function (series and continued fraction, as in Press et al., Numerical
Recipes, section 6.2). They are tested against SciPy when it is installed and
against closed forms otherwise (tests/test_mutation_grouping.py).

- `chi2_sf(x, k)`: P(X > x) for X ~ chi-square with k degrees of freedom.
- `chi2_isf(alpha, k)`: the upper-alpha critical value.
- `ncx2_sf(x, k, lam)`: tail of the noncentral chi-square, as the Poisson
  mixture sum_j Pois(j; lam / 2) * chi2_sf(x, k + 2 j).
- `noncentrality_for_power(k, alpha, power)`: the noncentrality at which a
  level-alpha chi-square test on k degrees of freedom has the given power.
"""

import math
from functools import lru_cache

_EPS = 1e-15
_TINY = 1e-300


def _lower_series(a, x):
    term = total = 1.0 / a
    ap = a
    for _ in range(10000):
        ap += 1
        term *= x / ap
        total += term
        if abs(term) < abs(total) * _EPS:
            break
    return total * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _upper_fraction(a, x):
    b = x + 1 - a
    c = 1 / _TINY
    d = 1 / b
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = _TINY if abs(d) < _TINY else d
        c = b + an / c
        c = _TINY if abs(c) < _TINY else c
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < _EPS:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def regularized_upper_gamma(a, x):
    """Q(a, x) = Gamma(a, x) / Gamma(a) for a > 0, x >= 0."""

    if x <= 0:
        return 1.0
    if x < a + 1:
        return max(0.0, 1.0 - _lower_series(a, x))
    return _upper_fraction(a, x)


def chi2_sf(x, k):
    """Upper tail P(chi2_k > x); k may be any positive number."""

    if k <= 0:
        raise ValueError("degrees of freedom must be positive")
    if not math.isfinite(x):
        return 0.0 if x > 0 else 1.0
    return regularized_upper_gamma(k / 2.0, max(x, 0.0) / 2.0)


@lru_cache(maxsize=4096)
def chi2_isf(alpha, k):
    """Critical value c with P(chi2_k > c) = alpha (bisection; accurate to about 1e-10)."""

    low, high = 0.0, max(10.0, 2.0 * k)
    while chi2_sf(high, k) > alpha:
        high *= 2
    for _ in range(200):
        middle = (low + high) / 2
        if chi2_sf(middle, k) > alpha:
            low = middle
        else:
            high = middle
        if high - low < 1e-12 * max(1.0, high):
            break
    return (low + high) / 2


def ncx2_sf(x, k, lam):
    """Upper tail of the noncentral chi-square (Poisson mixture of central tails)."""

    if lam <= 0:
        return chi2_sf(x, k)
    half = lam / 2.0
    centre = int(half)
    total = 0.0
    # Sum outward from the Poisson mode until the remaining weight is negligible.
    for direction in (range(centre, centre + 100000), range(centre - 1, -1, -1)):
        for j in direction:
            weight = math.exp(-half + j * math.log(half) - math.lgamma(j + 1)) if half > 0 else float(j == 0)
            total += weight * chi2_sf(x, k + 2 * j)
            if weight < 1e-17 and j != centre:
                break
    return min(1.0, total)


@lru_cache(maxsize=4096)
def noncentrality_for_power(k, alpha, power):
    """Smallest noncentrality lam with ncx2_sf(chi2_isf(alpha, k), k, lam) >= power (bisection)."""

    critical = chi2_isf(alpha, k)
    low, high = 0.0, 10.0
    while ncx2_sf(critical, k, high) < power:
        high *= 2
    for _ in range(100):
        middle = (low + high) / 2
        if ncx2_sf(critical, k, middle) < power:
            low = middle
        else:
            high = middle
        if high - low < 1e-9:
            break
    return high
