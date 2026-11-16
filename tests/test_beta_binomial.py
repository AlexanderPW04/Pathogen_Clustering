"""Independent numerical checks of genomics/beta_binomial.py (Sprint 5 Phase 2).

The expected values do not reuse the module's rising-factorial formula. They
come from the beta-binomial's definition as a Beta mixture of binomials
(numerical integration), the gamma-function form of the pmf, closed forms
(Beta(1, 1) gives a discrete uniform; one read gives a Bernoulli), published
Wilson values (Newcombe 1998), an analytic likelihood-ratio bound, and SciPy
when it happens to be installed (it is not a project dependency).
"""

import math
import unittest

import numpy as np

from genomics import beta_binomial as bb
from genomics.frequency_uncertainty import wilson_interval

try:  # optional independent reference; skipped where SciPy is not installed
    from scipy import stats as scipy_stats
except ImportError:  # pragma: no cover - depends on the environment
    scipy_stats = None


def mixture_pmf(k, n, p, rho, points=200001):
    """P(K = k) = integral of Binomial(k; n, q) x Beta(q; alpha, beta) dq (trapezoid rule)."""

    alpha, beta = p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho
    q = np.linspace(0, 1, points)[1:-1]
    log_beta_density = ((alpha - 1) * np.log(q) + (beta - 1) * np.log1p(-q)
                        - (math.lgamma(alpha) + math.lgamma(beta) - math.lgamma(alpha + beta)))
    binomial = math.comb(n, k) * q ** k * (1 - q) ** (n - k)
    return float(np.trapezoid(binomial * np.exp(log_beta_density), q))


def gamma_form_pmf(k, n, p, rho):
    alpha, beta = p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho
    log_value = (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
                 + math.lgamma(k + alpha) + math.lgamma(n - k + beta) - math.lgamma(n + alpha + beta)
                 - math.lgamma(alpha) - math.lgamma(beta) + math.lgamma(alpha + beta))
    return math.exp(log_value)


class ParameterizationTest(unittest.TestCase):
    def test_shapes_and_rho_round_trip(self):
        alpha, beta = bb.shape_parameters(0.3, 0.1)
        self.assertAlmostEqual(alpha, 2.7)
        self.assertAlmostEqual(beta, 6.3)
        self.assertAlmostEqual(bb.rho_from_shapes(alpha, beta), 0.1)
        self.assertAlmostEqual(alpha / (alpha + beta), 0.3)

    def test_invalid_rho_is_rejected(self):
        for value in [-0.01, 1, 1.5, float("nan"), float("inf"), "x", None]:
            with self.assertRaises(ValueError):
                bb.check_rho(value)
        self.assertEqual(bb.check_rho("0.02"), 0.02)

    def test_invalid_counts_are_rejected(self):
        for k, n in [(-1, 10), (11, 10), (1, 0), (2.5, 10)]:
            with self.assertRaises(ValueError):
                bb.log_pmf(k, n, 0.2, 0.01)
        with self.assertRaises(ValueError):
            bb.log_pmf(1, 10, 1.2, 0.01)


class ProbabilityTest(unittest.TestCase):
    def test_matches_the_beta_mixture_definition(self):
        for k, n, p, rho in [(3, 20, 0.3, 0.1), (0, 15, 0.2, 0.05), (12, 30, 0.5, 0.2)]:
            self.assertAlmostEqual(math.exp(bb.log_pmf(k, n, p, rho)), mixture_pmf(k, n, p, rho), places=7)

    def test_matches_the_gamma_function_form(self):
        for k, n, p, rho in [(40, 300, 0.1, 0.02), (150, 300, 0.5, 0.005), (900, 5000, 0.2, 0.01)]:
            self.assertAlmostEqual(math.exp(bb.log_pmf(k, n, p, rho)) / gamma_form_pmf(k, n, p, rho), 1, places=8)

    def test_beta_one_one_is_discrete_uniform(self):
        # p = 0.5, rho = 1/3 gives alpha = beta = 1: every count 0..n has probability 1/(n + 1).
        np.testing.assert_allclose(bb.pmf_all(9, 0.5, 1 / 3), np.full(10, 0.1), rtol=1e-12)

    def test_one_read_is_bernoulli_for_any_rho(self):
        for rho in [0, 0.01, 0.5]:
            np.testing.assert_allclose(bb.pmf_all(1, 0.3, rho), [0.7, 0.3], rtol=1e-12)

    def test_rho_zero_is_exactly_binomial(self):
        n, p = 40, 0.15
        binomial = [math.comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(n + 1)]
        np.testing.assert_allclose(bb.pmf_all(n, p, 0.0), binomial, rtol=1e-11)
        self.assertAlmostEqual(bb.log_pmf(6, n, p, 0.0), math.log(binomial[6]), places=10)

    def test_tiny_rho_is_continuous_with_the_binomial(self):
        n, p = 300, 0.2
        # The true relative gap is about rho k^2 / (2p) (2e-7 in the far tail at k = 300), not rounding.
        np.testing.assert_allclose(bb.pmf_all(n, p, 1e-12), bb.pmf_all(n, p, 0.0), rtol=1e-6, atol=1e-300)

    def test_boundary_frequencies_are_point_masses(self):
        self.assertEqual(bb.log_pmf(0, 50, 0.0, 0.02), 0.0)
        self.assertEqual(bb.log_pmf(3, 50, 0.0, 0.02), -math.inf)
        self.assertEqual(bb.log_pmf(50, 50, 1.0, 0.02), 0.0)
        self.assertEqual(bb.log_pmf(49, 50, 1.0, 0.0), -math.inf)

    def test_large_depth_is_finite_and_sums_to_one(self):
        probabilities = bb.pmf_all(50000, 0.003, 0.02)
        self.assertTrue(np.isfinite(probabilities).all())
        # 50,000 summed log terms carry about 1e-9 of rounding error.
        self.assertAlmostEqual(float(probabilities.sum()), 1.0, places=7)
        self.assertTrue(math.isfinite(bb.log_pmf(150, 50000, 0.003, 0.02)))

    def test_mean_and_variance_from_the_distribution(self):
        # Moments summed from the pmf (not from the closed-form functions).
        for n, p, rho in [(20, 0.2, 0.0), (100, 0.05, 0.02), (300, 0.5, 0.1)]:
            k = np.arange(n + 1)
            probabilities = bb.pmf_all(n, p, rho)
            mean = float(np.sum(k * probabilities))
            variance = float(np.sum((k - mean) ** 2 * probabilities))
            self.assertAlmostEqual(mean, n * p, places=8)
            self.assertAlmostEqual(variance, n * p * (1 - p) * (1 + (n - 1) * rho), places=6)
            self.assertAlmostEqual(float(bb.count_variance(n, p, rho)), variance, places=6)
            self.assertAlmostEqual(float(bb.frequency_variance(n, p, rho)), variance / n ** 2, places=10)

    def test_variance_of_the_frequency_approaches_p_q_rho(self):
        # Deep sequencing cannot remove library-level variation.
        self.assertAlmostEqual(float(bb.frequency_variance(10 ** 9, 0.2, 0.02)), 0.2 * 0.8 * 0.02, places=9)
        self.assertEqual(float(bb.design_effect(1, 0.3)), 1.0)
        self.assertAlmostEqual(float(bb.design_effect(300, 0.02)), 6.98)

    @unittest.skipIf(scipy_stats is None, "SciPy is not installed (optional reference only)")
    def test_matches_scipy_betabinom(self):
        for n, p, rho in [(20, 0.2, 0.02), (300, 0.05, 0.001), (5000, 0.3, 0.02)]:
            alpha, beta = bb.shape_parameters(p, rho)
            reference = scipy_stats.betabinom.pmf(np.arange(n + 1), n, alpha, beta)
            np.testing.assert_allclose(bb.pmf_all(n, p, rho), reference, atol=1e-12)


class PredictionIntervalTest(unittest.TestCase):
    def test_binomial_quantiles(self):
        # Independent: cumulative binomial sums with math.comb.
        n, p = 100, 0.2
        cumulative = np.cumsum([math.comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(n + 1)])
        expected = (int(np.argmax(cumulative >= 0.025)), int(np.argmax(cumulative >= 0.975)))
        self.assertEqual(bb.count_prediction_interval(p, n, 0.0), expected)

    def test_overdispersion_widens_the_prediction_range(self):
        narrow = bb.count_prediction_interval(0.2, 300, 0.0)
        wide = bb.count_prediction_interval(0.2, 300, 0.02)
        self.assertLess(wide[0], narrow[0])
        self.assertGreater(wide[1], narrow[1])
        probabilities = bb.pmf_all(300, 0.2, 0.02)
        self.assertGreaterEqual(float(probabilities[wide[0]:wide[1] + 1].sum()), 0.95)


class ScoreIntervalTest(unittest.TestCase):
    def test_rho_zero_equals_the_phase1_wilson_interval(self):
        x = np.array([0, 1, 4, 15, 81, 20, 263])
        n = np.array([20, 29, 20, 148, 263, 20, 263])
        lower, upper = bb.score_interval(x, n, 0.0)
        wilson_lower, wilson_upper = wilson_interval(x, n)
        np.testing.assert_allclose(lower, wilson_lower, atol=1e-15)
        np.testing.assert_allclose(upper, wilson_upper, atol=1e-15)
        # Newcombe (1998) table: 81/263 -> 0.2553 to 0.3662.
        self.assertAlmostEqual(float(lower[4]), 0.2553, places=4)
        self.assertAlmostEqual(float(upper[4]), 0.3662, places=4)

    def test_equals_wilson_at_the_effective_depth(self):
        # 40 / 200 with rho = 0.02: effective depth 200 / 4.98. Independent: Wilson formula by hand.
        n_eff, p_hat, z = 200 / 4.98, 0.2, bb.Z_95
        centre = (p_hat + z * z / (2 * n_eff)) / (1 + z * z / n_eff)
        half = z / (1 + z * z / n_eff) * math.sqrt(p_hat * (1 - p_hat) / n_eff + z * z / (4 * n_eff ** 2))
        lower, upper = bb.score_interval([40], [200], 0.02)
        self.assertAlmostEqual(float(lower[0]), centre - half, places=12)
        self.assertAlmostEqual(float(upper[0]), centre + half, places=12)

    def test_widens_with_rho_and_stays_in_range(self):
        widths = []
        for rho in [0, 0.001, 0.01, 0.05, 0.3, 0.99]:
            lower, upper = bb.score_interval([0, 3, 150, 300], [300] * 4, rho)
            self.assertTrue(((lower >= 0) & (upper <= 1) & (lower <= upper)).all())
            self.assertEqual(float(lower[0]), 0.0)
            self.assertEqual(float(upper[3]), 1.0)
            widths.append(float(upper[1] - lower[1]))
        self.assertTrue(all(a < b for a, b in zip(widths, widths[1:])))

    def test_invalid_inputs_give_missing_intervals(self):
        lower, upper = bb.score_interval([np.nan, 5, 3, 2.5, 2, 2], [20, 4, 0, 20, 20, 20],
                                         [0.01, 0.01, 0.01, 0.01, -0.1, 1.0])
        self.assertTrue(np.isnan(lower).all() and np.isnan(upper).all())

    def test_expected_width_matches_interval_at_whole_counts(self):
        lower, upper = bb.score_interval([60], [300], 0.02)
        self.assertAlmostEqual(float(bb.expected_score_width(0.2, 300, 0.02)), float(upper[0] - lower[0]), places=12)


class LikelihoodRatioIntervalTest(unittest.TestCase):
    def test_zero_count_binomial_bound_is_analytic(self):
        # L(p) = (1 - p)^n, so the upper bound solves -2 n log(1 - p) = 3.8415.
        lower, upper = bb.likelihood_ratio_interval(0, 20, 0.0)
        self.assertEqual(lower, 0.0)
        self.assertAlmostEqual(upper, 1 - math.exp(-bb.CHI2_1_95 / 40), places=8)

    def test_contains_the_estimate_and_widens_with_rho(self):
        binomial = bb.likelihood_ratio_interval(60, 300, 0.0)
        overdispersed = bb.likelihood_ratio_interval(60, 300, 0.02)
        self.assertLess(binomial[0], 0.2)
        self.assertGreater(binomial[1], 0.2)
        self.assertLess(overdispersed[0], binomial[0])
        self.assertGreater(overdispersed[1], binomial[1])
        self.assertEqual(bb.likelihood_ratio_interval(300, 300, 0.02)[1], 1.0)


if __name__ == "__main__":
    unittest.main()
