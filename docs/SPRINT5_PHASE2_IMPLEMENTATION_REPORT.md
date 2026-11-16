# Sprint 5 Phase 2 — Implementation Report

**Overdispersion-aware pathogen variant analysis.** Branch `sprint4-integration`, HEAD `b641533`.
All Phase 1 and Phase 2 work is local and uncommitted. Nothing was committed, pushed, merged, or
deployed, and no paid API calls were made (`AI_ENRICHMENT_ENABLED=false` in every run).
Labels: **verified** means it was run or read in this session. **Synthetic** means simulated data
with known truth.

## 1. Summary of implementation

Phase 1 left one documented gap. Under an overdispersed simulation, the 95% Wilson (binomial)
interval covered the true allele frequency only about 52% of the time at 300×. Phase 2 adds:

1. **A numerically stable beta-binomial model** (`genomics/beta_binomial.py`). It covers the
   pmf, moments, the design effect, prediction ranges for future counts, a score interval for the
   mean frequency, and a likelihood-ratio interval used only for evaluation. It is exactly binomial
   at ρ = 0 and needs no new dependency.
2. **Overdispersion estimated only from technical replicates** (`genomics/overdispersion.py`).
   - Two samples that share a `Specimen_ID`, collection date, and location are one replicate pair.
   - One ρ is estimated per pair, from the sites that pass the criteria of Tonkin-Hill et al.
     (2021).
   - The estimator is a method of moments with a site bootstrap. The published plug-in
     likelihood is reproduced for comparison but **not** used, because it underestimates ρ by
     about half (§3, §5).
3. **Three kinds of result, kept strictly apart** (`genomics/overdispersion_analysis.py`):
   - the Phase 1 Wilson baseline, always kept and never overwritten;
   - calibrated beta-binomial intervals, only for the replicated samples, using the upper 95%
     bound of ρ;
   - hypothetical sensitivity scenarios, using user-specified ρ values, always labelled
     "Hypothetical; not calibrated".
   Datasets without replicates, including the published `mmc2.xlsx`, never receive a ρ estimate.
4. **A pre-registered simulation study** (`validation/simulate_overdispersion.py`). It
   reproduces the Phase 1 undercoverage with identical random draws and shows the new interval
   corrects it under the model's assumptions. It also identifies which estimators and interval
   rules are supported.
5. **Website integration.** Section 9 now opens with the panel "Sequencing Uncertainty and
   Variant Frequency". It also adds:
   - two figures;
   - a calibrated overlay on the time series;
   - a hypothetical-ρ setting;
   - status rows;
   - provenance entries;
   - four new downloads;
   - 11 appended columns in `sprint5_vaf_observations.csv`.

Out of scope and not implemented: MAGOS clustering, lineage identification, annotation,
syn/nonsyn classification, FASTQ alignment, new callers, effect prediction, outbreak prediction,
and new Claude features.

## 2. Exact mathematical models

At one position in one library, n is the usable depth, k the alternative read count, p the
underlying mean alternative allele frequency, and ρ the overdispersion (intra-library
correlation).

| | Binomial (Phase 1) | Beta-binomial (Phase 2) |
|---|---|---|
| Model | k ~ Binomial(n, p) | q ~ Beta(α, β), k \| q ~ Binomial(n, q) |
| Parameters | p | p, ρ, with α = p(1−ρ)/ρ, β = (1−p)(1−ρ)/ρ, ρ = 1/(α+β+1) |
| E[k] | np | np (unchanged) |
| Var(k) | np(1−p) | np(1−p)·[1 + (n−1)ρ] |
| Var(k/n) | p(1−p)/n | p(1−p)/n·[1 + (n−1)ρ] → p(1−p)ρ as n → ∞ |

The factor 1 + (n−1)ρ is Kish's design effect for clustered binomial data. Because Var(k/n)
approaches p(1−p)ρ instead of 0, deeper sequencing cannot remove variation that arises before
the reads are counted.

**Numerically stable pmf.** P(k) = C(n,k)·B(k+α, n−k+β)/B(α, β) is evaluated as

log P(k) = log C(n,k) + Σ_{i<k} log(p(1−ρ)+iρ) + Σ_{i<n−k} log((1−p)(1−ρ)+iρ) − Σ_{i<n} log((1−ρ)+iρ).

The 1/ρ factors cancel. As a result:

- ρ = 0 gives the binomial exactly;
- tiny ρ does not cancel huge gamma functions;
- p = 0 or 1 gives point masses;
- 50,000× sums to 1 within 2×10⁻⁹.

**Interval for the mean frequency, "beta-binomial score interval".** It solves
(k/n − p)² = z²·p(1−p)/n·[1+(n−1)ρ]. This is the Wilson interval with n replaced by the
effective depth n/[1+(n−1)ρ], so ρ = 0 gives the Phase 1 interval exactly (tested). It is a
quasi-likelihood score interval. It is labelled **approximate**, and its coverage was checked by
exact beta-binomial sums (§5); it is not called an exact beta-binomial interval.

**What is estimated versus predicted.** The score interval is a confidence interval for the
library's mean allele frequency p. `count_prediction_interval(p, n, ρ)` answers a different
question: where would the count of one *future* library fall? It is used only as a goodness-of-fit
diagnostic for replicate pairs (`Share_In_Prediction_Range`) and is never shown as a confidence
interval.

**Multiallelic sites.** An allele-versus-rest count is a marginal of a Dirichlet-multinomial,
which is beta-binomial with the same ρ. Intervals are therefore made per allele, as in Phase 1.
Multiallelic positions are excluded from calibration.

**Limits of the model.** It describes extra variation between libraries around the specimen's
mean. It does **not** correct:

- mapping artifacts;
- primer bias that is the same in every library;
- reference mismatch;
- contamination;
- biological differences between samples, sites, or dates.

## 3. Parameter-estimation methodology

**Source methods (read in full from the Europe PMC full-text XML, not summaries).**

- **Tonkin-Hill et al. 2021 (eLife 10:e66857).**
  - 1,181 SARS-CoV-2 samples were sequenced in duplicate, with separate reverse transcription,
    PCR, and library preparation.
  - For each replicate pair, the variable sites with coverage above 100× in both replicates and a
    mean VAF above 1% and below 90% were fitted with one beta-binomial ρ.
  - The fit is maximum likelihood by grid search on the "approximate likelihood"
    Π_i Π_j BB(x_ij; n_ij, p = (x1j+x2j)/(n1j+n2j), ρ).
  - ρ ≤ 0.02 was used as an acceptable-discordance QC threshold, and 58% of samples met it.
  - ρ correlated strongly with Ct (higher Ct, less input RNA, more discordance).
- **MAGOS (Ahmadinejad et al. 2022, MBE).**
  - Each subclone cluster's VAFs are modelled as Beta(α, β), with α+β = mean depth and
    α/(α+β) = mean VAF.
  - This gives var(v) = v̄(1−v̄)/(ē+1), which is binomial-like with **no overdispersion
    parameter**.
  - Depth weights the clustering distance. It is a tumor-subclone method and is context for later
    clustering work only.

**Estimator used: method of moments on replicate pairs.** For site j, with X = x1+x2,
N = n1+n2, and c_i = [1+(n_i−1)ρ]/n_i, the model gives two identities:

- E[(f1−f2)²] = p_j(1−p_j)(c1+c2);
- E[X(N−X)] = p_j(1−p_j)(N² − S), with S = N + ρ[n1(n1−1)+n2(n2−1)].

Both were verified in the tests by exact enumeration of the joint pmf. ρ̂ solves

Σ_j (f1j−f2j)² = Σ_j X_j(N_j−X_j)(c1j+c2j)/(N_j²−S_j).

Both sides have the same expectation at the true ρ whatever the site frequencies are, so p_j
never has to be fitted first. The right side increases with ρ, so the root is found by bisection
on [0, 0.5]. ρ̂ = 0 is reported when the replicates agree at least as well as independent read
sampling predicts.

**Uncertainty of ρ̂.** A percentile bootstrap over sites gives 1,000 resamples (seed 20261115),
reported as `Rho_Lower_95` and `Rho_Upper_95`.

**ρ used for intervals: the upper bound `Rho_Upper_95`, not ρ̂.** This rule was chosen by the
simulation study (§5, C3) because it carries the estimation uncertainty into the interval.

**Why the published plug-in likelihood is not used.** Fitting p_j to the same two counts absorbs
one degree of freedom per site (an incidental-parameter problem), so the replicates appear about
half as discordant as they are. On identical simulated data its relative bias was −46% to −69%,
and intervals built from it covered only 0.90–0.93. It is still computed and reported
(`Plugin_Likelihood_Rho`, plus a likelihood-ratio statistic against ρ = 0). For Tonkin-Hill et
al.'s purpose, ranking samples against a threshold, this scale bias matters less.

**Pooling assumption.** One ρ per library pair assumes the extra variation comes from the
library process (input molecules, reverse transcription, amplification), which every site of that
library shares. Sites in the same amplicon may be correlated; the site bootstrap does not capture
this. No ρ is pooled across specimens, because Tonkin-Hill et al. found ρ varies with Ct. A
"borrowed" value (the largest calibrated upper bound) is offered only as a labelled sensitivity
scenario for unreplicated samples.

**Diagnostics reported per specimen** (`sprint5_overdispersion_calibration.csv`):

- sites considered, sites used, and sites lacking measured counts;
- depth and VAF ranges;
- the replicate VAF correlation and the mean absolute difference;
- ρ̂ and its bootstrap range;
- the share of resamples at the search limit;
- the plug-in ρ and its likelihood-ratio statistic;
- the share of replicate counts inside the 95% prediction range;
- the status, and the reason when there is no estimate.

## 4. Calibration requirements

A ρ is estimated, and calibrated intervals are produced, only when **all** of the following hold:

1. **Two samples share a `Specimen_ID`.** It is an optional column in the SNV table or the pooled
   metadata, and must be one value per sample. Each replicate keeps its own `Sample_ID`.
2. **Both samples have the same collection date and location, and the same reference sequence.**
   A pair collected a week apart is rejected: "not technical replicates". If both dates are
   missing, the declared `Specimen_ID` is accepted.
3. **At least 10 informative sites.** Each site must meet all of these:
   - a single alternative allele;
   - depth above 100 in both replicates;
   - a mean VAF of the two replicates between 1% and 90%;
   - **measured** counts in both replicates. These are supplied alternate counts, or the
     A/C/G/T count of an allele that one replicate did not report. Reconstructed counts and
     unknown counts are not used. An unreported allele without A/C/G/T counts is "count unknown",
     never 0.
4. **The estimate is stable.** ρ̂ must be below 0.5, and at most 5% of bootstrap resamples may
   reach 0.5.

Every unmet condition is reported with its reason, and those samples keep the binomial baseline.
Exactly two replicates per specimen are supported; three or more are reported as unsupported
rather than pooled ad hoc.

## 5. Validation and simulation results

All results here are **synthetic** (seed 20261115). They show that the methods behave as intended
*when the beta-binomial assumptions hold*. They do not show that real wastewater libraries follow
the model, and calibration on synthetic data does not prove calibration on real samples.

**Criteria fixed before the run** (recorded in `docs/SPRINT5_PHASE2_PROGRESS.md` first):

- **C1, known ρ.** Over depth {20, 50, 100, 300} × p {0.01 … 0.99}, by exact sums: mean coverage
  0.95 ± 0.01 and minimum ≥ 0.90.
- **C2, Phase 1 scenario.** Wilson's undercoverage is reproduced, and the new interval's simulated
  coverage is within 4 Monte Carlo SE of its exact coverage.
- **C3, estimated ρ.** Mean coverage ≥ 0.93 and every cell ≥ 0.90.
- **C4.** The published estimator is evaluated on the same data.

**C1: ρ known, exact coverage (no Monte Carlo error)**

| True ρ | Wilson mean / min | **Score interval** mean / min | Likelihood-ratio mean / min |
|---|---|---|---|
| 0 | 0.947 / 0.911 (pass) | **0.947 / 0.911 (pass; identical to Wilson)** | 0.949 / 0.852 (fail) |
| 0.001 | 0.935 / 0.906 (fail) | **0.950 / 0.907 (pass)** | 0.947 / 0.834 (fail) |
| 0.005 | 0.893 / 0.771 (fail) | **0.953 / 0.916 (pass)** | 0.955 / 0.850 (fail) |
| 0.02 | 0.800 / 0.524 (fail) | **0.951 / 0.919 (pass)** | 0.951 / 0.829 (fail) |
| 0.05 | 0.716 / 0.352 (fail) | **0.952 / 0.940 (pass)** | 0.953 / 0.896 (fail) |

The likelihood-ratio interval inverts the exact beta-binomial likelihood. It undercovers at small
counts near the boundaries, so it was **not promoted**.

**C2: the Phase 1 scenario, reproduced exactly** (ρ = 0.02, 20,000 runs, same seed and same
random draws; the Wilson coverages match Phase 1's table to 0.0)

| Depth | p | Wilson coverage | New interval simulated / exact | Wilson width → new width |
|---|---|---|---|---|
| 20 | 0.05 / 0.2 / 0.5 | 0.897 / 0.915 / 0.918 | 0.964 / 0.973 / 0.969 (exact 0.963 / 0.974 / 0.967) | 0.32 → 0.37 (p 0.2) |
| 100 | 0.05 / 0.2 / 0.5 | 0.812 / 0.719 / 0.721 | 0.954 / 0.952 / 0.946 | 0.15 → 0.26 (p 0.2) |
| 300 | 0.05 / 0.2 / 0.5 | **0.553 / 0.533 / 0.522** | **0.952 / 0.951 / 0.948** | 0.089 → 0.231 (p 0.2) |

Every simulated value lies within 4 SE of its exact coverage (SE ≈ 0.0015). The interval is not
calibrated merely by being wider: its exact coverage is close to 0.95 throughout, and the widening
it needs grows with depth, as the design effect predicts.

**C3: ρ estimated from simulated replicate pairs** (400 calibration datasets per cell; coverage of
test observations at depth 20–300 and p 0.05/0.2/0.5, computed exactly given each estimate)

| Sites | Rule | Mean coverage over true ρ 0–0.05 | Worst cell | Supported |
|---|---|---|---|---|
| 5 | Point estimate ρ̂ | 0.924–0.952 | 0.893 | **No** (at ρ = 0.05) |
| 5 | Upper 95% bound | 0.958–0.965 | 0.926 | Yes |
| 10 | Point estimate | 0.939–0.951 | 0.922 | Yes |
| 10 | **Upper 95% bound (used)** | 0.958–0.972 | 0.925 | **Yes** |
| 20 | Upper 95% bound | 0.956–0.972 | 0.925 | Yes |
| 50 | Upper 95% bound | 0.955–0.969 | 0.924 | Yes |

Moment-estimator properties:

- Relative bias is within ±4.1% for every design with ρ ≥ 0.01 (target ≤ 10%).
- RMSE at ρ = 0.02 is 0.013 with 5 sites, 0.0097 with 10, 0.0068 with 20, and 0.0042 with 50.
- The bootstrap upper bound was above the true ρ in 73–76% of datasets with 5 sites, about 84%
  with 10, 89% with 20, and 92–94% with 50. It is a rough range, not an exact 97.5% bound.
- Because 5 sites gave such an unreliable bound, **at least 10 informative sites are required.**
- There were no unstable datasets.

**C4: published plug-in likelihood versus moments** (20 sites, 60 datasets; identical data)

| True ρ | Moment mean | Plug-in mean (relative bias) | Coverage using the plug-in estimate |
|---|---|---|---|
| 0.002 | 0.0020 | 0.0006 (−69%) | 0.934 |
| 0.01 | 0.0103 | 0.0051 (−49%) | 0.916 |
| 0.02 | 0.0201 | 0.0107 (−46%) | 0.910 |
| 0.05 | 0.0467 | 0.0266 (−47%) | 0.897 |

**Sensitivity: a wrong assumption** (exact, p = 0.2, 300×)

| True ρ \ assumed ρ | 0 | 0.002 | 0.02 | 0.05 |
|---|---|---|---|---|
| 0 | 0.949 | 0.989 | 1.000 | 1.000 |
| 0.02 | 0.535 | 0.659 | **0.952** | 0.996 |
| 0.05 | 0.362 | 0.461 | 0.810 | **0.959** |

A hypothetical ρ gives nominal coverage only when it happens to equal the truth. This is why
sensitivity results are never called calibrated.

**Misspecification.** When ρ varies between sites (Gamma, coefficient of variation 1, mean 0.01,
0.02, or 0.05), the upper-bound rule with 50 sites still covered 0.961–0.968 at every depth.

**Boundary behaviour.** Intervals stay in [0, 1], with exact 0 at k = 0 and 1 at k = n. There are
no NaNs for valid input, and NaN is returned for invalid counts or ρ (tested). The p grid
includes 0.01 and 0.99, and the minimum C1 coverage there was 0.907.

**Runtime** (one core, this laptop):

| Task | Time |
|---|---|
| 8,000 score intervals | 0.2 ms |
| One calibration with 50 sites and 1,000 bootstrap resamples | 13 ms |
| Plug-in likelihood | 0.2 s |
| One LR interval at 300× | 1.5 ms |
| The full study | 90 s |

**Synthetic replicate demo** (`examples/integrated/make_replicate_example.py`, true values known):

| Specimen | True ρ | Result |
|---|---|---|
| SPEC_A | 0.01 | ρ̂ = 0.0105 (range 0.0064–0.0152), 37 sites; plug-in 0.0055 |
| SPEC_B | 0.002 | ρ̂ = 0.0021 (range 0.0006–0.0039), 28 sites; plug-in 0.0009 |
| SPEC_C | — | "too few informative sites" (5 sites) |
| SPEC_D | — | "different collection dates, so they are not technical replicates" |

## 6. Comparison with Phase 1 Wilson intervals

- **The Wilson interval is unchanged.** Every Phase 1 column of `sprint5_vaf_observations.csv` is
  identical, the 37 columns in the same order. The web test recomputes Wilson bounds from the
  downloaded counts.
- **The two models agree exactly at ρ = 0.** The score interval reduces to Wilson there (C1 row 1
  and a unit test).
- **Under overdispersion, Wilson is too narrow, and more so with depth:**
  - ρ = 0.02: coverage 0.92 at 20×, 0.72 at 100×, 0.53 at 300×;
  - the beta-binomial interval stays near 0.95;
  - the required widening is √(design effect): about 1.2× at 20× and 2.6× at 300× for ρ = 0.02.
- **Interpretation on the page:**
  - The Wilson interval is labelled the baseline and the narrowest plausible interval.
  - A calibrated interval is shown only where technical replicates support it.
  - Otherwise the page shows the binomial baseline plus clearly hypothetical scenarios.

## 7. Real-data compatibility

| Input | Phase 2 result (verified through `/review` → `/run`) |
|---|---|
| `mmc2.xlsx` (52 samples, 7,973 observations, reconstructed counts) | Runs; **no ρ estimate** (no replicates). The panel shows "Technical overdispersion could not be estimated from this dataset. The binomial interval is shown as a baseline, and any overdispersion sensitivity analysis is hypothetical." Every `Rho` is blank. 7,970 Wilson intervals (unchanged). The sensitivity table has 23,910 rows (7,970 × 3 hypothetical values). Median width ratios versus Wilson: 1.55× (ρ 0.002), 4.61× (0.02), 8.24× (0.05). |
| Synthetic pooled example (no `Specimen_ID`) | Not calibrated, with the required statement; sensitivity only |
| Synthetic technical replicates | Calibrated for SPEC_A and SPEC_B (127 intervals); explicit reasons for SPEC_C, SPEC_D, and WW_SINGLE |
| Replicates without dates | Still calibrated (`Specimen_ID` declared); temporal trends say dates are unavailable |
| Missing depth | No intervals; reasons shown; no errors |
| Sprint 2 + 3 combined | Panel present; cross-source and Sprint 2 outputs unchanged |
| Invalid ρ setting (e.g. "2") | A message is shown; no crash |

**Original data.** The `mmc2.xlsx` copies were read only (see §11 for checksums). The project
never writes to an input file.

## 8. Files and functions modified

**New files**

| File | Contents |
|---|---|
| `genomics/beta_binomial.py` | `check_rho`, `shape_parameters`, `design_effect`, `count_mean`, `count_variance`, `frequency_variance`, `log_pmf`, `log_pmf_all`, `pmf_all`, `count_prediction_interval`, `score_interval`, `expected_score_width`, `likelihood_ratio_interval` |
| `genomics/overdispersion.py` | `moment_estimate` (vectorized over bootstrap weights), `bootstrap_estimates`, `plugin_log_likelihood`, `plugin_likelihood_estimate`, `share_in_prediction_range`, `estimate_overdispersion`, `replicate_pairs`, `calibration_sites`, `calibrate_replicates`; constants for the site criteria, `MINIMUM_SITES = 10`, and statuses |
| `genomics/overdispersion_analysis.py` | `parse_rho_values`, `add_overdispersion_columns`, `sensitivity_scenarios`, `sensitivity_table`, `model_comparison`, `example_observations`, `run_overdispersion_analysis`; labels and the required statements |
| `validation/simulate_overdispersion.py` | The full study: `known_rho_table`, `phase1_scenario_table`, `estimated_rho_table`, `plugin_bias_table`, `sensitivity_table`, `heterogeneous_rho_table`, `runtime_table`, `criteria_summary` |
| `examples/integrated/make_replicate_example.py` | The synthetic replicate generator, which writes `synthetic_replicate_snvs.csv`, `synthetic_replicate_base_counts.csv`, and `synthetic_replicate_reference.fasta` |
| `tests/test_beta_binomial.py` | 23 tests |
| `tests/test_overdispersion.py` | 19 tests |
| `tests/test_overdispersion_simulation.py` | 8 tests |
| `tests/test_sprint5_phase2_web.py` | 10 tests |
| `docs/SPRINT5_PHASE2_PROGRESS.md` | The progress record |
| `docs/SPRINT5_PHASE2_IMPLEMENTATION_REPORT.md` | This report |

**Modified files** (Phase 1 files that were already uncommitted are edited in place)

| File | Change |
|---|---|
| `genomics/depth_aware_analysis.py` | Calls `run_overdispersion_analysis` and attaches `Specimen_ID`. Adds the four tables and the two figures. Adds calibrated bars to time-series points. Adds `result["depth_aware"]["uncertainty"]`. |
| `genomics/figures.py` | New `uncertainty_models_figure` and `replicate_concordance_figure`. `vaf_time_series_figure` draws the calibrated band and its legend entry only when present. |
| `genomics/integrated_analysis.py` | Default `overdispersion_rhos`. Two status rows. Two provenance parameters. |
| `wastewater/processed_snv.py` | Optional `Specimen_ID`, which metadata may fill. It must be one value per sample and is kept as an extra column. |
| `genomics/pooled_observations.py` | `Specimen_ID` added to `PROVENANCE_COLUMNS`. |
| `genomics/column_detection.py` | `SPECIMEN_ALIASES`, with a `Specimen_ID` field for the SNV and pooled-metadata roles. |
| `app.py` | `overdispersion_rho_values` in `default_values`, validated in `check_settings`, passed as `settings["overdispersion_rhos"]`. |
| `integrated_view.py` | `OVERDISPERSION_LIMITATIONS`, `uncertainty_view`, a new stat, a finding line, the calibration table, the figure order, and the downloads. |
| `templates/index.html` | The `uncertainty_panel` macro, rendered first in section 9. The advanced-settings fieldset "Sequencing uncertainty (Sprint 5)". Panel style. |
| `templates/_toolbar.html` | Re-runs keep `overdispersion_rho_values`. |
| `dashboard_guide.py` | Guides for `uncertainty_models` and `replicate_concordance`. |

**Intentional schema changes**

- `sprint5_vaf_observations.csv` gains 11 columns **after** the 37 Phase 1 columns:
  - `Specimen_ID`
  - `Baseline_Model`
  - `Overdispersion_Model`
  - `Rho`
  - `Rho_Source`
  - `Calibration_Status`
  - `BetaBinomial_Lower`
  - `BetaBinomial_Upper`
  - `BetaBinomial_Width`
  - `BetaBinomial_Interval_Type`
  - `Uncertainty_Limitations`
- Per-row labels are short so the `mmc2` table stays under the 8 MB download limit (§10). The full
  wording is on the page and in this report.
- New downloads:
  - `sprint5_uncertainty_models.csv`;
  - `sprint5_overdispersion_sensitivity.csv` (one row per observation and hypothetical ρ);
  - `sprint5_overdispersion_calibration.csv` (only with replicates);
  - `sprint5_calibration_sites.csv` (only with replicates).
- No other download changed.

## 9. User-interface improvements

Section 9 now opens with **"Sequencing Uncertainty and Variant Frequency"**, before any technical
table. The panel shows:

1. **A plain explanation**: "Sequencing depth affects how precisely a mutation's frequency can be
   estimated. Additional sequencing variation may make the uncertainty larger than expected from
   read counts alone."
2. **A badge and headline**, either:
   - **Empirically calibrated**: the number of specimens, the ρ range, and how many observations
     have calibrated intervals; or
   - **Not calibrated**: "Technical overdispersion could not be estimated from this dataset. The
     binomial interval is shown as a baseline, and any overdispersion sensitivity analysis is
     hypothetical."
3. **A compact model comparison**: one row per model (binomial baseline, calibrated
   beta-binomial, each hypothetical ρ, borrowed ρ), with ρ, calibration status, number of
   observations, median width, and width relative to the basic interval.
4. **Up to 6 example observations**, spread over depth: frequency, depth, basic interval,
   overdispersion-aware interval (calibrated, or marked hypothetical), and whether it is
   calibrated.
5. **"Why some uncertainty results are unavailable"**: a collapsible list of reasons per
   specimen.

Figures:

- The time series stays first. It shows calibrated bands only for replicated samples.
- "Binomial versus beta-binomial interval width":
  - left, theory at VAF 0.20 against depth;
  - right, each observation's Wilson width against its overdispersion-aware width.
- "Technical replicate agreement" appears only with replicates.
- The Phase 1 figures are kept.

The new setting, "Hypothetical overdispersion values (rho)", is in Advanced settings. Its help
text says these values are scenarios, not estimates, and explains how to declare technical
replicates.

## 10. Performance measurements

**Measured locally** (Windows working set, `OPENBLAS_NUM_THREADS=1`). Every run used the guided web
path `/review` → `/run` on `mmc2.xlsx` and the same scripts as Phase 1. Render was **not** measured.

| Run (`mmc2.xlsx`) | Phase 1 (reported) | Phase 2 (measured, project `.venv`) |
|---|---|---|
| One run, peak memory | 287 MB | **296 MB** (two runs: 296 and 296) |
| Same interpreter, Phase 2 step disabled (A/B) | — | 291 MB and 309 MB |
| Four runs in one process | 307 MB | **306 MB** (296 → 303 → 305 → 306) |
| 512 MB cap and 10% of one CPU core | 282 MB, 107.5 s | **285 MB, 110.4 s** (job committed 288 MB) |
| One run, time (unconstrained) | 10.7 s | 10.4 s |
| Results page | 2.50 MB | 2.65 MB (panel, two figures) |
| `sprint5_vaf_observations.csv` | about 5.0 MB | 6.27 MB (11 columns added) |
| `sprint5_overdispersion_sensitivity.csv` | — | 3.86 MB (23,910 rows) |

**Interpretation:**

- **Memory.** With Phase 2 disabled in the same interpreter, the peak was 291–309 MB. With it
  enabled, it was 296 MB both times. Run-to-run variation (about ±10 MB) is therefore larger than
  any Phase 2 effect, and no memory regression could be detected. The highest point moved from the
  Sprint 3 `observation_status_figure` (still 287 MB) to the Phase 1 `vaf_depth` figure (295 MB),
  drawn while the Phase 2 tables are in memory. System Python (numpy 2.2.6, matplotlib 3.9.2)
  measured 304 MB for the same run; the interpreter alone accounts for about 8 MB.
- **Time.** The constrained run was 2.9 s (2.7%) slower, which is within normal variation.
  110.4 s is well inside gunicorn's 180 s timeout.
- **Regression found and fixed during integration.** The first version added long per-row text
  and made `mmc2`'s `sprint5_vaf_observations.csv` 8.6 MB. That is over the 8 MB per-file
  download limit, so the file silently disappeared from the downloads. The labels were shortened,
  and the Phase 1 web test for `mmc2` (which downloads this file) now passes again. No result or
  plot was dropped to meet a target.
- **Not measured:** Linux and Render memory, and Render's proxy timeout.

## 11. Full test results

**Command** (from the project folder; `FONTENELE_MMC2` points to `mmc2.xlsx`; MAFFT is on PATH at
`C:\Tools\mafft-win`, as in Phase 1):

```bash
python -m unittest discover -s tests -v
```

| Run | Result |
|---|---|
| Baseline before Phase 2 (verified at start; MAFFT not on PATH) | Ran 475 tests in 330.5 s — OK (skipped=2) |
| Intermediate full run after integration (MAFFT not on PATH) | Ran 535 tests in 301.9 s — OK (skipped=2) |
| **Final full run** (MAFFT on PATH; after the last edit, a caption change) | **Ran 535 tests in 297.1 s — OK (skipped=1)**: the only skip is the optional NCBI files test, as in Phase 1; all three real-`mmc2.xlsx` tests and the real-MAFFT test ran. An earlier run on the same code before the caption edit gave the same result (301.97 s). |
| Phase 2 tests only, project `.venv` (no SciPy) | Ran 60 tests in 32.8 s — OK (skipped=2: the optional SciPy check, and mmc2 without a path) |

**New tests (60):**

- `tests/test_beta_binomial.py` (23): parameterization; pmf against the mixture integral and the
  gamma form; closed forms; the binomial limit; boundaries; 50,000×; moments; SciPy; prediction
  ranges; score interval equal to Wilson at ρ = 0 and at the effective depth; monotonicity; invalid
  input; the LR interval's analytic bound.
- `tests/test_overdispersion.py` (19): the moment identities by exact enumeration; recovery of ρ
  and zero; bootstrap reproducibility; plug-in optimum against a brute-force grid; plug-in bias;
  statuses (too few, zero, unstable, upper-bound rule); replicate-pair rules; site rules (measured
  counts, base counts, unknown ≠ 0, reconstructed, multiallelic, depth, VAF); setting parser;
  integration (no invented ρ, calibration only for replicated samples, no counts → no interval).
- `tests/test_overdispersion_simulation.py` (8): fast versions of C1–C4 and the sensitivity table.
- `tests/test_sprint5_phase2_web.py` (10): through `/review` → `/run`. Covers replicates,
  no replicates, user and empty settings, an invalid setting, missing depth, replicates without
  dates, Sprint 2+3 combined, `mmc2`, and the `Specimen_ID` loader (table and metadata).

**Existing tests:** none was changed in Phase 2. All 475 pre-existing tests pass, including the
three that use the real `mmc2.xlsx` and the Phase 1 web test that downloads
`sprint5_vaf_observations.csv` for `mmc2`. That test caught the 8 MB download regression described
in §10.

**Output regression check (verified).** `run_integrated_analysis.py --example` was rerun with the
`.venv` interpreter and MAFFT, and compared with the Phase 1 output:

- every Phase 1 CSV is identical, apart from two added status rows and two added provenance
  parameters;
- the first 37 columns of `sprint5_vaf_observations.csv` are identical;
- all 23 Phase 1 figures are byte-identical.

This covers MAFFT alignment, Sprint 2 distances and screening, Sprint 3 QC, Yue–Clayton, PCoA,
cross-source matching, Phase 1 VAF provenance, and missing-coverage logic. Sprint 4 literature
code was not modified.

## 12. Limitations and assumptions

1. **Synthetic validation only.** The coverage results assume beta-binomial libraries. Real
   amplicon data can also vary in other ways:
   - primer dropout;
   - amplicon-specific bias;
   - position-specific errors;
   - correlation between nearby sites.
   No real replicated wastewater data were available in this repository, so **no real sample has
   been shown to be calibrated**.
2. **One ρ per library pair.** It is assumed common to every informative site. When ρ varied
   between sites, the upper-bound rule still covered 0.96–0.97 in simulation. Correlation between
   sites in one amplicon is not modelled, so the site bootstrap may understate the uncertainty
   of ρ.
3. **The ρ range from the bootstrap is approximate.** It contained the true ρ less often than
   nominal (about 84% with 10 sites, 92–94% with 50). The interval rule still met the coverage
   criterion.
4. **Calibration excludes some sites.** These are sites unreported in both replicates, below 1%,
   or above 90%. It uses only measured counts, so a dataset with reported frequencies alone (such
   as `mmc2`) cannot be calibrated even if replicates were declared. Sites skipped for lack of
   counts are counted, because dropping them can bias ρ low.
5. **Calibrated intervals apply only to the replicated samples.** ρ depends on input quantity
   (Ct) in Tonkin-Hill et al., so it is not transferred to other samples. The "borrowed" value is
   a labelled sensitivity scenario only.
6. **Exactly two replicates per specimen are supported.** Larger designs are reported as
   unsupported.
7. **The score interval is an approximation.** Its exact coverage was verified on the stated
   grid, but not for every depth and frequency. Wilson-style oscillation is present: the minimum
   was 0.907 at p = 0.01.
8. **What is and is not included.** The intervals describe read sampling plus library-to-library
   variation for one specimen. They exclude primer bias shared by every library, mapping and
   reference errors, contamination, and biological variation between samples, sites, and dates.
   They do not identify lineages, outbreaks, or transmission. No mutation clustering is done.
9. **Phase 1 limitations still apply.** Reconstructed `mmc2` counts are derived, not measured, and
   time series without A/C/G/T counts use detections only.

## 13. How to test the new features locally

From the project folder:

1. Regenerate the synthetic replicate demo (deterministic; the files are already in
   `examples/integrated/`):

   ```bash
   python examples/integrated/make_replicate_example.py
   ```

2. Run the simulation study (about 90 s). The output folder is optional:

   ```bash
   python validation/simulate_overdispersion.py validation_output
   ```

3. Run the Phase 2 tests:

   ```bash
   python -m unittest tests.test_beta_binomial tests.test_overdispersion tests.test_overdispersion_simulation tests.test_sprint5_phase2_web -v
   ```

4. In the website:
   - Choose **Pooled sequencing**.
   - Upload `synthetic_replicate_snvs.csv`, `synthetic_replicate_base_counts.csv`, and
     `synthetic_replicate_reference.fasta`.
   - Run, and open section 9. The panel should say **Empirically calibrated**, and the replicate
     figure and calibration table should appear.
   - Then upload `synthetic_pooled_snvs.csv` or `mmc2.xlsx`. The panel should say
     **Not calibrated**, with the required statement.
   - To change or remove the hypothetical values, use **Advanced settings → Sequencing
     uncertainty**.
5. **To declare your own technical replicates:**
   - Give each library its own `Sample_ID`.
   - Add a `Specimen_ID` column, with the same value for the two libraries of one specimen,
     either in the SNV table or in the pooled metadata.
   - Keep the same date and location.
   - Supply alternate read counts, and ideally an A/C/G/T table, so an allele missed in one
     replicate is counted rather than unknown.

## 14. Requirements for proceeding to mutation clustering

Before MAGOS-style clustering of mutation frequencies:

1. **Real replicate data.** Obtain technical-replicate libraries for at least some wastewater
   specimens, from the same protocol as the surveillance samples, with measured counts. Estimate
   ρ and check whether it depends on depth, input quantity, or amplicon. Without this, any
   clustering distance that uses VAF variance rests on an unverified noise model.
2. **A variance model per observation.** MAGOS sets var(v) = v̄(1−v̄)/(ē+1), with no
   overdispersion. For wastewater, the natural replacement is the beta-binomial variance
   p(1−p)/n·[1+(n−1)ρ], using a calibrated ρ. It should be decided and validated whether ρ is
   per library, per batch, or per amplicon.
3. **Biological versus technical variation.** Clustering groups mutations whose frequencies move
   together across samples and dates, which is biological signal. It must use technical ρ only as
   the noise term, never as a fitted between-sample variance.
4. **Coverage-aware absence.** Use supported non-detections from A/C/G/T counts. Unknown values
   must never be 0.
5. **Validation design.** Use simulated mixtures with known co-varying groups plus replicate
   noise. Pre-registered criteria are needed for cluster recovery, such as the adjusted Rand
   index. MAGOS's own acceptance threshold, J ≥ 0.80, comes from tumor data and must not be
   assumed for wastewater.
6. **No lineage claims.** Clusters are groups of co-varying mutations, not lineage
   identifications.

## 15. Render deployment considerations

- **No new dependency.** `requirements.txt` and the Dockerfile are unchanged. SciPy is used only
  by an optional test that is skipped where it is not installed (as in `.venv` and Docker).
- **Added cost.** Score intervals are vectorized (0.2 ms per 8,000 observations). Calibration runs
  only when replicates exist (about 13 ms per replicate pair). The plug-in comparison takes about
  0.2 s per pair, and its cost grows with depth. With many replicate pairs at very high depth this
  adds seconds, not minutes.
- **Gunicorn timeout.** It stays at 180 s. The constrained local run is in §10.
- **To verify on Render after review** (not done here):
  1. Upload `mmc2.xlsx` and confirm section 9 shows "Not calibrated". Confirm that
     `sprint5_vaf_observations.csv` and `sprint5_overdispersion_sensitivity.csv` download.
  2. Upload the replicate demo and confirm "Empirically calibrated" and the replicate figure.
  3. Watch the Render memory graph and the request time.
- **Nothing was pushed or deployed.** Render remains **unverified**.

## Acceptance criteria

| # | Mandatory criterion | Status | Evidence |
|---|---|---|---|
| 1 | Phase 1 and the current repository state audited | **Met** | §3 of the progress file: the git state, a 475-test baseline, the Phase 1 modules read, and both source papers read in full text |
| 2 | Numerically validated beta-binomial implementation | **Met** | `genomics/beta_binomial.py`; 23 tests against the mixture integral, the gamma form, closed forms, Newcombe, an analytic LR bound, and SciPy (§5, §11) |
| 3 | Defensible calibration procedure for replicated data | **Met** | Moment estimator with two identities verified exactly, site bootstrap, Tonkin-Hill site criteria, diagnostics; the published plug-in estimator evaluated and shown biased (§3, §5 C3–C4) |
| 4 | Never invents overdispersion when unsupported | **Met** | No ρ without a valid replicate pair and at least 10 measured sites. `mmc2` and the unreplicated examples have every `Rho` blank (web tests) |
| 5 | Validated interval procedure for supported conditions | **Met (synthetic)** | Score interval passes C1 for every ρ. The upper-bound rule passes C3 for every design with 10 or more sites (§5). Not validated on real replicates (§12) |
| 6 | Wilson intervals remain the baseline | **Met** | Phase 1 columns are identical: CLI comparison, and Wilson recomputed in the web test (§6) |
| 7 | Sensitivity distinguished from empirical calibration | **Met** | Separate columns, files, labels ("Hypothetical; not calibrated"), and the badge (§9); the sensitivity table shows wrong assumptions mis-cover (§5) |
| 8 | Simulations compare expected and observed coverage | **Met** | Exact versus simulated coverage (C1, C2), coverage versus nominal by design (C3), width, bias, RMSE, runtime (§5) |
| 9 | Phase 1 undercoverage investigated | **Met** | Same seed and identical draws: Wilson 0.52–0.55 at 300×, corrected to 0.95 (§5 C2) |
| 10 | Integrated into the existing website | **Met** | Section 9 panel, two figures, setting, status rows, downloads (§9); web tests |
| 11 | Unsupported datasets get explanations, not errors | **Met** | Required statement plus per-specimen reasons. Tested with no replicates, missing depth, missing dates, an invalid setting, and unusable pairs |
| 12 | `mmc2.xlsx` still runs | **Met** | Guided web path: HTTP 200, no figure errors, all downloads under 8 MB (§7, §10) |
| 13 | Sprints 2, 3, 4 and Phase 1 still work | **Met** | Ran 535 tests in 297.1 s — OK (skipped=1); CLI outputs and figures identical to Phase 1 (§11) |
| 14 | Memory and runtime evaluated | **Met** | A/B, four-run, and 512 MB / 10% CPU measurements (§10). Render not measured |
| 15 | New and existing tests pass without unexplained failures | **Met** | Ran 535 tests in 297.1 s — OK (skipped=1). The only skips are the optional NCBI files. No existing test was changed. Three expectations in the new tests were set from derived bounds (two numerical tolerances, one coverage pattern), explained in the progress file; no expected failures were added |
| 16 | Limitations and assumptions documented | **Met** | §2, §3, §4, §12 |
| 17 | No paid calls, commits, pushes, merges, or deployments | **Met** | `AI_ENRICHMENT_ENABLED=false` in every run; `git status` shows uncommitted changes only; nothing pushed |
| 18 | Implementation report completed | **Met** | This document |

All 18 mandatory criteria are met for the conditions tested. Criterion 5 is validated on synthetic
data that satisfy the beta-binomial assumptions. That is the standard the goal sets ("under simulated
datasets matching the stated beta-binomial assumptions"). It is **not** evidence that real wastewater
libraries are calibrated, and §12 says so explicitly. Render deployment remains unverified because
deploying was out of scope.

**SPRINT 5 PHASE 2 COMPLETE**
