"""Estimating beta-binomial overdispersion from technical sequencing replicates (Sprint 5 Phase 2).

Why replicates are required
---------------------------
Frequencies from different samples, dates, or genomic sites differ for
biological reasons, so their spread is not sequencing noise. Overdispersion
of the measurement process can only be estimated from repeated measurements
of the SAME specimen: technical replicates with separate library preparation
(Tonkin-Hill et al. 2021, eLife 10:e66857, "Beta-binomial modelling of
replicate samples"). Without such replicates no rho is estimated here.

Input design
------------
Two pooled samples (separate `Sample_ID`s) that share a `Specimen_ID` are
treated as two technical replicates of one specimen, provided they have the
same collection date and location AND (Phase 3) each has its own documented
`Library_ID` and the two differ. Matching metadata alone does not show that
the libraries were prepared separately: undocumented pairs, one library
listed twice, and exact duplicate records are reported as unverified and
receive no rho. One rho is estimated for each replicate
pair from all of its informative sites (as in Tonkin-Hill et al.). Pooling the
sites assumes the extra variation comes from the library process (input
molecules, reverse transcription, amplification), which is shared by every
site of that library; sites in the same amplicon may be correlated, which the
site bootstrap below does not capture.

Informative sites (Tonkin-Hill et al. criteria)
-----------------------------------------------
A single-alternative-allele site with depth above 100 in both replicates and
a mean frequency of the two replicates above 1% and below 90%. Both counts
must be measured: supplied alternate counts or A/C/G/T base counts. A
replicate that did not report the allele counts as 0 reads only when its
A/C/G/T counts say so; otherwise the site is skipped (dropping it can bias rho
low, so the number skipped is reported). Reconstructed counts are not used.

Estimators
----------
For site j with counts x1, x2, depths n1, n2, frequencies f = x/n, and a
common mean p_j, the beta-binomial model gives

    E[(f1 - f2)^2] = p_j (1 - p_j) (c1 + c2),   c_i = [1 + (n_i - 1) rho] / n_i
    E[X (N - X)]   = p_j (1 - p_j) (N^2 - S),
    with X = x1 + x2, N = n1 + n2, S = N + rho [n1 (n1 - 1) + n2 (n2 - 1)].

1. Moment estimator (used for calibration). rho solves

       sum_j (f1j - f2j)^2 = sum_j X_j (N_j - X_j) (c1j + c2j) / (N_j^2 - S_j).

   Both sides have the same expectation at the true rho for any p_j, so no
   site frequency has to be estimated first. The right side increases with
   rho, so the root is found by bisection. rho = 0 is reported when the
   replicates agree at least as well as independent read sampling predicts.
   A 95% range comes from resampling sites (percentile bootstrap). In the
   simulation study it contained the true rho less often than 95% (about
   85% with 10 sites, 92% with 50), so it is reported as a rough range.

Rho used for intervals: the 97.5th bootstrap percentile ("upper 95% bound"),
not the point estimate, so the intervals allow for the uncertainty of the
estimate. With rho estimated this way from 5-50 simulated sites, the
intervals' coverage was 0.955-0.972 (point estimate: 0.924-0.952).

2. Plug-in likelihood (Tonkin-Hill et al., reproduced for comparison). It
   maximizes prod_i prod_j BB(x_ij; n_ij, p = X_j / N_j, rho) over a grid of
   rho, then refines between grid points. Because p_j is fitted to the same
   two counts, the replicates look about half as discordant as they are, and
   rho is underestimated by roughly a factor of 2 (shown in
   validation/simulate_overdispersion.py). The authors used rho to rank
   samples against a threshold, where this scale bias matters less; it is
   not used here for intervals.
"""

import math

import numpy as np
import pandas as pd

from genomics import beta_binomial as bb
from genomics.pooled_observations import COUNT_COLUMNS
from genomics.unified_schema import NUCLEOTIDES


CALIBRATION_MIN_DEPTH = 100     # depth must be ABOVE this in both replicates (Tonkin-Hill et al.)
CALIBRATION_VAF_LOW = 0.01      # mean replicate frequency must be above 1% ...
CALIBRATION_VAF_HIGH = 0.90     # ... and below 90% (Tonkin-Hill et al.)
# Fewest informative sites for a usable estimate. In the simulation study (C3) the
# interval rule below met the coverage criterion from 5 sites, but with 5 sites the
# bootstrap upper bound of rho was above the true rho only 73-76% of the time; from
# 10 sites both the point-estimate and upper-bound rules passed, so 10 is required.
MINIMUM_SITES = 10
BOOTSTRAP_REPLICATES = 1000
BOOTSTRAP_SEED = 20261115
RHO_SEARCH_MAX = 0.5            # estimates this large are reported as unusable, not as calibrated
BOUNDARY_SHARE_LIMIT = 0.05     # unstable if more than 5% of bootstrap estimates hit RHO_SEARCH_MAX

ESTIMATED = "Estimated from technical replicates"
ESTIMATED_ZERO = "Estimated from technical replicates (no extra-binomial variation detected)"
TOO_FEW_SITES = "Not estimated: too few informative sites"
UNSTABLE = "Not estimated: unstable estimate"
NOT_REPLICATE_PAIR = "Not estimated: not a usable replicate pair"
UNVERIFIED_PAIR = "Not estimated: replicate libraries not documented (unverified)"
USABLE_STATUSES = [ESTIMATED, ESTIMATED_ZERO]

MEASURED_SOURCES = {"Supplied", "Base counts"}

SITE_COLUMNS = ["Specimen_ID", "Replicate_1", "Replicate_2", "Reference_ID", "Position", "Ref", "Alt",
                "Variant_Key", "Alt_Reads_1", "Depth_1", "Count_Source_1", "Alt_Reads_2", "Depth_2",
                "Count_Source_2", "VAF_1", "VAF_2", "Mean_VAF", "Used", "Exclusion_Reason"]

CALIBRATION_COLUMNS = [
    "Specimen_ID", "Replicate_1", "Replicate_2", "Collection_Date", "Location_Name", "Status",
    "Status_Detail", "Sites_Considered", "Sites_Used", "Sites_Without_Measured_Counts",
    "Depth_Min", "Depth_Median", "Depth_Max", "Mean_VAF_Min", "Mean_VAF_Max",
    "Replicate_VAF_Correlation", "Mean_Absolute_VAF_Difference", "Rho_Estimate",
    "Rho_Lower_95", "Rho_Upper_95", "Rho_Used_For_Intervals", "Bootstrap_At_Search_Limit",
    "Plugin_Likelihood_Rho", "Plugin_Likelihood_Ratio_vs_Binomial", "Share_In_Prediction_Range",
    "Estimator",
]


# ----------------------------------------------------------------------------
# Estimators on arrays of paired counts
# ----------------------------------------------------------------------------

def _moment_terms(x1, n1, x2, n2):
    x1, n1, x2, n2 = (np.asarray(value, dtype=float) for value in (x1, n1, x2, n2))
    discordance = (x1 / n1 - x2 / n2) ** 2
    total, depth = x1 + x2, n1 + n2
    spread = total * (depth - total)
    a = 1 / n1 + 1 / n2                              # c1 + c2 at rho = 0
    b = (n1 - 1) / n1 + (n2 - 1) / n2                # slope of c1 + c2 in rho
    e = depth ** 2 - depth                           # N^2 - S at rho = 0
    g = n1 * (n1 - 1) + n2 * (n2 - 1)                # slope of S in rho
    return discordance, spread, a, b, e, g


def moment_estimate(x1, n1, x2, n2, weights=None, upper=RHO_SEARCH_MAX, iterations=60):
    """Moment estimate of rho from replicate pairs; `weights` (B x m) gives B estimates at once.

    Returns a float (or an array of length B with weights). Values equal to
    `upper` mean the root lies at or beyond the search limit.
    """

    discordance, spread, a, b, e, g = _moment_terms(x1, n1, x2, n2)
    w = np.ones((1, len(discordance))) if weights is None else np.asarray(weights, dtype=float)
    target = w @ discordance

    def expected(rho):
        terms = spread * (a + b * rho[:, None]) / (e - g * rho[:, None])
        return np.sum(w * terms, axis=1)

    low = np.zeros(len(w))
    high = np.full(len(w), float(upper))
    at_zero = target <= expected(low)
    at_limit = target >= expected(high)
    for _ in range(iterations):
        middle = (low + high) / 2
        below = expected(middle) < target
        low = np.where(below, middle, low)
        high = np.where(below, high, middle)
    estimate = np.where(at_zero, 0.0, np.where(at_limit, float(upper), (low + high) / 2))
    return float(estimate[0]) if weights is None else estimate


def bootstrap_estimates(x1, n1, x2, n2, replicates=BOOTSTRAP_REPLICATES, seed=BOOTSTRAP_SEED):
    """Moment estimates of rho for `replicates` resamples of the sites (with replacement)."""

    m = len(np.asarray(x1))
    generator = np.random.default_rng(seed)
    weights = generator.multinomial(m, np.full(m, 1.0 / m), size=replicates)
    return moment_estimate(x1, n1, x2, n2, weights=weights)


def plugin_log_likelihood(x1, n1, x2, n2, rho):
    """Tonkin-Hill et al. approximate log-likelihood at one rho (p_j = pooled frequency)."""

    pooled = (np.asarray(x1) + np.asarray(x2)) / (np.asarray(n1) + np.asarray(n2))
    total = 0.0
    for counts, depths in ((x1, n1), (x2, n2)):
        for k, n, p in zip(counts, depths, pooled):
            total += bb.log_pmf(int(k), int(n), float(p), rho)
    return total


def plugin_likelihood_estimate(x1, n1, x2, n2, grid_points=41):
    """(rho, likelihood-ratio statistic against rho = 0) by grid search plus golden-section refinement."""

    grid = np.concatenate([[0.0], 10 ** np.linspace(-6, math.log10(RHO_SEARCH_MAX), grid_points)])
    values = np.array([plugin_log_likelihood(x1, n1, x2, n2, rho) for rho in grid])
    best = int(np.argmax(values))
    if 0 < best < len(grid) - 1:
        low, high = math.log10(max(grid[best - 1], 1e-7)), math.log10(grid[best + 1])
        ratio = (math.sqrt(5) - 1) / 2
        c, d = high - ratio * (high - low), low + ratio * (high - low)
        fc = plugin_log_likelihood(x1, n1, x2, n2, 10 ** c)
        fd = plugin_log_likelihood(x1, n1, x2, n2, 10 ** d)
        for _ in range(30):
            if fc < fd:
                low, c, fc = c, d, fd
                d = low + ratio * (high - low)
                fd = plugin_log_likelihood(x1, n1, x2, n2, 10 ** d)
            else:
                high, d, fd = d, c, fc
                c = high - ratio * (high - low)
                fc = plugin_log_likelihood(x1, n1, x2, n2, 10 ** c)
        rho = 10 ** ((low + high) / 2)
        best_value = max(plugin_log_likelihood(x1, n1, x2, n2, rho), values[best])
        if best_value == values[best]:
            rho = grid[best]
    else:
        rho, best_value = grid[best], values[best]
    return float(rho), float(2 * (best_value - values[0]))


def share_in_prediction_range(x1, n1, x2, n2, rho, level=0.95):
    """Share of replicate counts inside the 95% beta-binomial prediction range at the pooled frequency.

    A goodness-of-fit diagnostic only: the pooled frequency is estimated from
    the same counts, so the share is somewhat optimistic.
    """

    pooled = (np.asarray(x1) + np.asarray(x2)) / (np.asarray(n1) + np.asarray(n2))
    inside = []
    for counts, depths in ((x1, n1), (x2, n2)):
        for k, n, p in zip(counts, depths, pooled):
            lowest, highest = bb.count_prediction_interval(float(p), int(n), rho, level)
            inside.append(lowest <= k <= highest)
    return float(np.mean(inside)) if inside else float("nan")


def estimate_overdispersion(sites, minimum_sites=MINIMUM_SITES, replicates=BOOTSTRAP_REPLICATES,
                            seed=BOOTSTRAP_SEED, include_plugin=True):
    """Estimate rho for one replicate pair from its USED sites; return a dict of results and diagnostics."""

    used = sites[sites["Used"]] if "Used" in sites else sites
    x1, n1 = used["Alt_Reads_1"].to_numpy(float), used["Depth_1"].to_numpy(float)
    x2, n2 = used["Alt_Reads_2"].to_numpy(float), used["Depth_2"].to_numpy(float)
    result = {"Sites_Used": len(used), "Rho_Estimate": np.nan, "Rho_Lower_95": np.nan,
              "Rho_Upper_95": np.nan, "Rho_Used_For_Intervals": np.nan, "Bootstrap_At_Search_Limit": np.nan,
              "Plugin_Likelihood_Rho": np.nan, "Plugin_Likelihood_Ratio_vs_Binomial": np.nan,
              "Share_In_Prediction_Range": np.nan, "Estimator": "Moment estimator with site bootstrap"}
    if len(used) == 0:
        return {**result, "Status": TOO_FEW_SITES, "Status_Detail": "No informative sites."}
    depths = np.concatenate([n1, n2])
    f1, f2 = x1 / n1, x2 / n2
    result.update({
        "Depth_Min": int(depths.min()), "Depth_Median": float(np.median(depths)), "Depth_Max": int(depths.max()),
        "Mean_VAF_Min": float(((f1 + f2) / 2).min()), "Mean_VAF_Max": float(((f1 + f2) / 2).max()),
        "Replicate_VAF_Correlation": float(np.corrcoef(f1, f2)[0, 1]) if len(used) > 2 else np.nan,
        "Mean_Absolute_VAF_Difference": float(np.mean(np.abs(f1 - f2))),
    })
    if len(used) < minimum_sites:
        return {**result, "Status": TOO_FEW_SITES,
                "Status_Detail": (f"{len(used)} informative site(s); at least {minimum_sites} are needed for an "
                                  "estimate whose intervals passed the simulation study.")}
    estimate = moment_estimate(x1, n1, x2, n2)
    boot = bootstrap_estimates(x1, n1, x2, n2, replicates, seed)
    at_limit = float(np.mean(boot >= RHO_SEARCH_MAX))
    result.update({"Rho_Estimate": estimate, "Rho_Lower_95": float(np.quantile(boot, 0.025)),
                   "Rho_Upper_95": float(np.quantile(boot, 0.975)), "Bootstrap_At_Search_Limit": at_limit})
    if include_plugin:
        plugin_rho, statistic = plugin_likelihood_estimate(x1, n1, x2, n2)
        result.update({"Plugin_Likelihood_Rho": plugin_rho, "Plugin_Likelihood_Ratio_vs_Binomial": statistic})
    if estimate >= RHO_SEARCH_MAX or at_limit > BOUNDARY_SHARE_LIMIT:
        return {**result, "Status": UNSTABLE,
                "Status_Detail": (f"The estimate reaches the search limit rho = {RHO_SEARCH_MAX} in "
                                  f"{at_limit:.0%} of site resamples; the replicates disagree too much for a "
                                  "stable beta-binomial fit (check for sample swaps or failed libraries).")}
    result["Rho_Used_For_Intervals"] = result["Rho_Upper_95"]
    result["Share_In_Prediction_Range"] = share_in_prediction_range(x1, n1, x2, n2, estimate)
    if estimate == 0:
        return {**result, "Status": ESTIMATED_ZERO,
                "Status_Detail": ("The replicates agree as closely as independent read sampling predicts "
                                  "(rho = 0); the intervals use the upper 95% bound of rho.")}
    return {**result, "Status": ESTIMATED,
            "Status_Detail": "Intervals use the upper 95% bound of rho, which allows for estimation uncertainty."}


# ----------------------------------------------------------------------------
# Building replicate pairs from the Sprint 5 evidence table
# ----------------------------------------------------------------------------

def _base_counts(position_table):
    """{(Sample_ID, Reference_ID, Position): (depth, {base: count})} when A/C/G/T counts exist."""

    if position_table is None or not position_table.attrs.get("has_base_counts", False):
        return {}
    keys = zip(position_table["Sample_ID"], position_table["Reference_ID"], position_table["Position"])
    counts = position_table[COUNT_COLUMNS].to_numpy()
    depth = position_table["Depth"].to_numpy()
    return {key: (int(d), dict(zip(NUCLEOTIDES, (int(v) for v in row))))
            for key, d, row in zip(keys, depth, counts)}


# Phase 3 audit: a shared Specimen_ID, date, and location do not by themselves show
# that two records are separately prepared libraries (they could be one library
# listed twice, or a re-run of one library). Calibration therefore also needs a
# documented, distinct Library_ID for each sample.
MISSING_LIBRARY = ("No Library_ID documents that the two samples are separately prepared libraries, so they are "
                   "treated as unverified replicates and not used for calibration.")
SAME_LIBRARY = ("Both samples have the same Library_ID: they are repeated measurements of one library, not "
                "independent technical replicates.")
DUPLICATE_RECORDS = ("The two samples report identical read counts and depths at every allele, which indicates "
                     "duplicate records of one measurement rather than separately prepared libraries.")
UNVERIFIED_PROBLEMS = (MISSING_LIBRARY, SAME_LIBRARY, DUPLICATE_RECORDS)


def _identical_measurements(group, samples):
    """True when both samples report the same alleles with the same counts and depths (duplicate records)."""

    columns = ["Variant_Key", "Alt_Read_Count", "Total_Depth"]
    if not all(column in group for column in columns):
        return False
    first, second = (group[group["Sample_ID"].eq(sample)][columns].astype(str)
                     .drop_duplicates().sort_values(columns).reset_index(drop=True) for sample in samples)
    return len(first) > 0 and first.equals(second)


def replicate_pairs(evidence):
    """Specimens with replicate samples: list of (specimen, samples, problem or "").

    A pair is usable only if exactly two samples share the Specimen_ID, they have
    the same collection date, location, and reference, each has its own
    documented Library_ID, the two Library_IDs differ, and the records are not
    exact duplicates. Problems in UNVERIFIED_PROBLEMS mean the replicate claim
    is undocumented rather than contradicted.
    """

    if "Specimen_ID" not in evidence:
        return []
    table = evidence[evidence["Specimen_ID"].fillna("").astype(str).str.strip().ne("")]
    libraries = (table["Library_ID"].fillna("").astype(str).str.strip() if "Library_ID" in table
                 else pd.Series("", index=table.index))
    pairs = []
    for specimen, group in table.groupby("Specimen_ID", sort=True):
        samples = sorted(group["Sample_ID"].unique())
        problem = ""
        if len(samples) != 2:
            problem = (f"{len(samples)} sample(s) share this Specimen_ID; exactly two technical replicates "
                       "are supported.")
        else:
            info = group.drop_duplicates("Sample_ID")
            dates = info["Collection_Date"].astype("string").fillna("<missing>")
            library = libraries[group.index].groupby(group["Sample_ID"]).first()
            if dates.nunique() > 1:
                problem = "The two samples have different collection dates, so they are not technical replicates."
            elif info["Location_Name"].astype(str).str.casefold().nunique() > 1:
                problem = "The two samples have different locations, so they are not technical replicates."
            elif group.groupby("Sample_ID")["Reference_ID"].agg(lambda v: tuple(sorted(set(v)))).nunique() > 1:
                problem = "The two samples use different reference sequences."
            elif library.eq("").any():
                problem = MISSING_LIBRARY
            elif library.nunique() == 1:
                problem = SAME_LIBRARY
            elif _identical_measurements(group, samples):
                problem = DUPLICATE_RECORDS
        pairs.append((specimen, samples, problem))
    return pairs


def calibration_sites(evidence, position_table, specimen, samples):
    """One row per alternative allele observed in either replicate, with the counts and why it is (not) used."""

    first, second = samples
    lookup = _base_counts(position_table)
    usable = evidence[evidence["Sample_ID"].isin(samples) & evidence["Counted_Once"]]
    by_key = {(row["Sample_ID"], row["Variant_Key"]): row for row in usable.to_dict("records")}
    alleles = (usable.groupby(["Reference_ID", "Position"])["Alt"].nunique())
    rows = []
    for (reference_id, position, ref, alt, key), _ in usable.groupby(
            ["Reference_ID", "Position", "Ref", "Alt", "Variant_Key"], sort=True):
        row = {"Specimen_ID": specimen, "Replicate_1": first, "Replicate_2": second,
               "Reference_ID": reference_id, "Position": int(position), "Ref": ref, "Alt": alt,
               "Variant_Key": key, "Used": False, "Exclusion_Reason": ""}
        values, reasons = [], []
        for sample in samples:
            record = by_key.get((sample, key))
            if record is not None:
                source = record["Alt_Count_Source"]
                if source not in MEASURED_SOURCES:
                    reasons.append(f"{sample}: counts are {str(source).lower()}, not measured")
                elif not str(record["Interval_Status"]).startswith("Calculated"):
                    reasons.append(f"{sample}: {record['Interval_Status']}")
                elif not record["Analysis_Included"]:
                    reasons.append(f"{sample}: excluded by quality control")
                values.append((record["Alt_Read_Count"], record["Total_Depth"], source))
            else:
                cell = lookup.get((sample, reference_id, int(position)))
                if cell is None:
                    reasons.append(f"{sample}: allele not reported and no A/C/G/T counts (count unknown)")
                    values.append((pd.NA, pd.NA, "Unknown"))
                else:
                    values.append((cell[1][alt], cell[0], "Base counts (allele not reported)"))
        (x1, d1, s1), (x2, d2, s2) = values
        row.update({"Alt_Reads_1": x1, "Depth_1": d1, "Count_Source_1": s1,
                    "Alt_Reads_2": x2, "Depth_2": d2, "Count_Source_2": s2})
        if alleles.get((reference_id, position), 1) > 1:
            reasons.append("more than one alternative allele at this position")
        if not reasons:
            if pd.isna(d1) or pd.isna(d2) or d1 <= CALIBRATION_MIN_DEPTH or d2 <= CALIBRATION_MIN_DEPTH:
                reasons.append(f"depth not above {CALIBRATION_MIN_DEPTH} in both replicates")
            else:
                f1, f2 = x1 / d1, x2 / d2
                mean = (f1 + f2) / 2
                row.update({"VAF_1": f1, "VAF_2": f2, "Mean_VAF": mean})
                if not CALIBRATION_VAF_LOW < mean < CALIBRATION_VAF_HIGH:
                    reasons.append(f"mean frequency {mean:.3f} outside {CALIBRATION_VAF_LOW:.0%} to "
                                   f"{CALIBRATION_VAF_HIGH:.0%}")
        row["Used"] = not reasons
        row["Exclusion_Reason"] = "; ".join(reasons)
        rows.append(row)
    table = pd.DataFrame(rows, columns=SITE_COLUMNS)
    for column in ["Alt_Reads_1", "Depth_1", "Alt_Reads_2", "Depth_2"]:
        table[column] = pd.array(pd.to_numeric(table[column], errors="coerce"), dtype="Float64").round().astype("Int64")
    table["Used"] = table["Used"].astype(bool)
    return table


def calibrate_replicates(evidence, position_table, minimum_sites=MINIMUM_SITES):
    """Return (calibration table, site table) for every specimen with replicate samples."""

    rows, site_tables = [], []
    for specimen, samples, problem in replicate_pairs(evidence):
        info = evidence[evidence["Sample_ID"].isin(samples)].iloc[0]
        base = {"Specimen_ID": specimen, "Replicate_1": samples[0] if samples else "",
                "Replicate_2": samples[1] if len(samples) > 1 else "",
                "Collection_Date": info["Collection_Date"], "Location_Name": info["Location_Name"]}
        if problem:
            status = UNVERIFIED_PAIR if problem in UNVERIFIED_PROBLEMS else NOT_REPLICATE_PAIR
            rows.append({**base, "Status": status, "Status_Detail": problem, "Sites_Considered": 0,
                         "Sites_Used": 0, "Sites_Without_Measured_Counts": 0})
            continue
        sites = calibration_sites(evidence, position_table, specimen, samples)
        site_tables.append(sites)
        unmeasured = int(sites["Exclusion_Reason"].str.contains("count unknown|not measured").sum())
        estimate = estimate_overdispersion(sites, minimum_sites)
        if unmeasured and estimate["Status"] in USABLE_STATUSES:
            estimate["Status_Detail"] += (f" {unmeasured} site(s) were skipped because a count was not measured "
                                          "in one replicate; if those sites disagree more, rho is underestimated.")
        rows.append({**base, "Sites_Considered": len(sites), "Sites_Without_Measured_Counts": unmeasured, **estimate})
    calibration = pd.DataFrame(rows, columns=CALIBRATION_COLUMNS)
    sites = pd.concat(site_tables, ignore_index=True) if site_tables else pd.DataFrame(columns=SITE_COLUMNS)
    return calibration, sites
