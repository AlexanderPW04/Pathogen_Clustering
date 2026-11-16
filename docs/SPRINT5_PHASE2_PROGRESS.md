# Sprint 5 Phase 2 — Progress Checklist

Persistent record for the Sprint 5 Phase 2 goal loop (overdispersion-aware pathogen variant
analysis). Read this file first when resuming. Labels: **Verified** (run or read in this session),
**Inferred**, **Unverified**.

Scratchpad for diagnostics (not part of the project):
`C:\Users\ALEXAN~1\AppData\Local\Temp\claude\a7271a137b0bdb68\922255f3-0871-46d4-aedd-17d549d92b7a\scratchpad\p2`

## Milestone status

| Milestone | Status | Attempt |
|---|---|---|
| 0. Verify Phase 1, baseline | **Done** | 1 |
| 1. Beta-binomial foundation | **Done** | 1 |
| 2. Overdispersion calibration | **Done** | 1 |
| 3. Overdispersion-aware intervals | **Done** | 1 |
| 4. Simulation study | **Done** | 1 |
| 5. Website integration | **Done** | 1 (download-size regression found and fixed) |
| 6. Regression, performance | **Done** | 1 |
| Final report | **Done** (`docs/SPRINT5_PHASE2_IMPLEMENTATION_REPORT.md`) | 1 |

---

## Milestone 0 — Baseline (completed)

**Git state (verified, November 15, 2026):** branch `sprint4-integration`, HEAD `b641533`. The Phase 1
work is uncommitted: 17 modified files (+875/−141 lines) and untracked `docs/`, five
`genomics/` modules (`collection_dates`, `depth_aware_analysis`, `frequency_uncertainty`,
`vaf_trends`, `variant_evidence`), four test files, and `validation/simulate_vaf_intervals.py`.
Phase 2 adds files and edits only what integration needs; nothing is reset, stashed, or committed.

**Baseline tests (verified):** `python -m unittest discover -s tests -v` with
`AI_ENRICHMENT_ENABLED=false`, `GEONAMES_USERNAME=` (empty), `OPENBLAS_NUM_THREADS=1`, and
`FONTENELE_MMC2` pointing to a scratch copy of `mmc2.xlsx`: **Ran 475 tests in 330.5 s — OK
(skipped=2)**. The two skips are the optional NCBI files (as in Phase 1) and the real-MAFFT test
(`mafft` is not on this shell's PATH; Phase 1 reported 1 skip). This is an environment difference,
not a code change. System Python 3.11.9 (numpy 2.2.6, pandas 2.2.2) has SciPy 1.13.1; the project
`.venv` (numpy 2.4.6, pandas 2.3.3) and `requirements.txt` do **not**, so Phase 2 code must not
import SciPy. Tests use SciPy only as an optional independent check.

**Phase 1 implementation (verified by reading code):**

- Variant data model: `genomics/variant_evidence.py` (`build_variant_evidence`, 37
  `EVIDENCE_COLUMNS`). Read-count provenance is in `Alt_Count_Source`, with the values Supplied,
  Base counts, Reconstructed, Frequency only, and No depth. The other provenance fields are
  `Ref_Count_Source`, `Depth_Source`, and `Interval_Status`.
- Wilson interval: `genomics/frequency_uncertainty.py` (`wilson_interval`, `expected_interval_width`,
  `depth_dependence_table`).
- Time series: `genomics/vaf_trends.py` (`known_points`, `series_summary`, `choose_series`,
  `unknown_dates`).
- Orchestration: `genomics/depth_aware_analysis.py` (`run_depth_aware_analysis`), called from
  `genomics/integrated_analysis.py`. Its outputs are 4 tables and 3 figures, and it fills
  `result["depth_aware"]`.
- QC used: `Analysis_Included`, `Counted_Once`, `Duplicate_Status`, `Frequency_Check`. Repeated
  records are counted once.
- Reference verification: `genomics/quality_control.reference_verification`.
- Routes: `/review` → `/run` (`app.py`). Dashboard section 9 is built by `integrated_view.py`
  (`SECTION_SPRINTS["depth-aware"]`) and rendered by `templates/index.html`.
- Fixtures: `examples/integrated/` synthetic pooled SNVs, base counts, and reference.
  `tests/test_sprint5_depth_aware.py` holds the end-to-end Flask tests.
- Docker/Render: gunicorn uses 1 worker, 4 threads, and `--timeout 180`. The target is 512 MB and
  0.1 CPU. The Phase 1 measurements were 287 MB for one `mmc2` run, 307 MB after four runs, and
  107.5 s under the 512 MB / 10 % CPU cap.

**Source methods read (verified, Europe PMC full-text XML, saved in the scratchpad):**

- Tonkin-Hill et al. 2021 (eLife, PMC8363274), Methods "Beta-binomial modelling of replicate
  samples":
  - 1,181 samples were sequenced in duplicate. Each replicate had separate RT, PCR, and library
    preparation.
  - For each replicate pair, the variable sites used had coverage above 100× in both replicates
    and a mean VAF above 1 % and below 90 %.
  - One ρ is fitted per replicate pair by maximum likelihood, found by grid search. It uses the
    "approximate likelihood" Π_i Π_j BB(x_ij; n_ij, p = (x1j+x2j)/(n1j+n2j), ρ).
  - ρ ≤ 0.02 was used as an acceptable-discordance QC threshold; 58 % of samples met it.
  - ρ correlated strongly with Ct.
- MAGOS (Ahmadinejad et al. 2022, MBE, PMC9260306):
  - Each subclone cluster's VAFs are modelled as Beta(α, β), with α+β = mean depth ē and
    α/(α+β) = mean VAF v̄.
  - This gives var(v) = v̄(1−v̄)/(ē+1), which is binomial-like.
  - There is **no overdispersion parameter**. The variance model is used for clustering
    (later work; out of scope).
  - It is a tumor-subclone method, so for pathogen surveillance it is used only as context.

**Scientific assumptions recorded before implementation:**

1. One ρ is estimated per pair of technical replicates of one specimen, as in Tonkin-Hill. Pooling
   across loci assumes that the extra variation comes from the library process (input molecules,
   amplification), which is shared by every site in that library. Sites in the same amplicon may
   be correlated; this is a documented limitation.
2. Different dates, locations, or specimens are never treated as replicates. Each replicate needs
   its own `Sample_ID`; the replicates share a `Specimen_ID`.
3. Calibration uses only measured counts (Supplied or Base counts). Reconstructed counts are
   excluded, and so are sites whose count is unknown in either replicate. An unreported allele
   counts as 0 reads only when the A/C/G/T counts show it.
4. Multiallelic positions are excluded from calibration. For the intervals, an allele-versus-rest
   count at a multiallelic site is a marginal of a Dirichlet-multinomial, which is beta-binomial
   with the same ρ.
5. Datasets without replicates (including `mmc2.xlsx`) get no ρ estimate. Any ρ used for them is a
   user-specified hypothetical sensitivity value.

**Performance criteria, fixed BEFORE running the Phase 2 simulation study (not tuned afterwards):**

- **C1, known ρ, exact sums (no Monte Carlo error).** The grid is depth {20, 50, 100, 300} × p {0.01,
  0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95, 0.98, 0.99} × ρ {0, 0.001, 0.005, 0.02, 0.05}. Across each ρ
  the mean coverage must be 0.95 ± 0.01 and the minimum at least 0.90, the same criterion Phase 1
  applied to Wilson under its own model.
- **C2, Phase 1 scenario.** The Phase 1 setup (ρ = 0.02, depths 20–300, p ∈ {0.05, 0.2, 0.5},
  20,000 runs) must reproduce Wilson's under-coverage. With the true ρ, the simulated coverage of
  the new interval must lie within 4 Monte Carlo SE of its exact coverage, and that exact coverage
  must meet C1.
- **C3, estimated ρ (replicates).**
  - Designs: m ∈ {5, 10, 20, 50} calibration sites; true ρ ∈ {0, 0.002, 0.01, 0.02, 0.05}.
  - Each test observation is scored on depth {20, 50, 100, 300} × p {0.05, 0.2, 0.5}. Coverage given
    the estimated ρ is computed exactly, then averaged over 400 simulated calibration datasets.
  - A design is **supported** if mean coverage is at least 0.93 and every cell is at least 0.90.
    Designs that fail are disabled (no calibrated interval).
  - The relative bias of ρ̂ is reported, with a target of at most 10 % for ρ ≥ 0.01.
- **C4.** Every estimator, including the published plug-in likelihood, is evaluated on the same
  simulated data. An estimator that fails C3 is not used for calibrated intervals.
- A method is never called calibrated only because its intervals are wider.

---

## Milestone 1 — Beta-binomial foundation (completed)

**Code:** `genomics/beta_binomial.py` (new; no SciPy).

- Parameterization: mean p and overdispersion ρ (0 ≤ ρ < 1), with α = p(1−ρ)/ρ and
  β = (1−p)(1−ρ)/ρ.
- The pmf is computed as a product of rising factorials in which 1/ρ cancels. It is exact at
  ρ = 0 (the binomial), has no gamma-function cancellation at tiny ρ, and no overflow at 50,000×.
- Functions:
  - `log_pmf`, `log_pmf_all`, `pmf_all`;
  - `count_mean`, `count_variance`, `frequency_variance`, `design_effect`;
  - `count_prediction_interval`, which gives the range for a FUTURE count, not a CI;
  - `score_interval`, the main interval;
  - `likelihood_ratio_interval`, evaluated only;
  - `check_rho`.

**Tests:** `tests/test_beta_binomial.py`, 23 tests (verified OK, with SciPy and without it in `.venv`).
The independent references are:

- the Beta-mixture definition, by numerical integration;
- the gamma-function form;
- Beta(1,1), which gives a discrete uniform;
- n = 1, which gives a Bernoulli;
- moments summed from the pmf;
- Newcombe's 81/263 → 0.2553–0.3662;
- the analytic LR bound for x = 0;
- SciPy `betabinom`, where installed (max difference 1e-13).

Two tolerances were set from derived error bounds, not loosened blindly:

- 50,000 summed log terms carry about 1e-9 of rounding error;
- at ρ = 1e-12 the true gap from the binomial is about ρk²/(2p) ≈ 2e-7.

## Milestone 2 — Calibration from technical replicates (completed)

**Code:** `genomics/overdispersion.py` (new).

- Design: two samples that share a `Specimen_ID`, date, and location form one replicate pair.
  `Specimen_ID` was added as an optional SNV or metadata column. The loader checks it, it is
  passed through the provenance columns, and column detection recognizes its common synonyms.
- Sites follow the Tonkin-Hill criteria, applied verbatim:
  - depth above 100 in both replicates;
  - mean VAF above 1 % and below 90 %;
  - measured counts only, either Supplied or from A/C/G/T counts for an allele that was not
    reported;
  - multiallelic positions are excluded.
- **Finding (verified by simulation):**
  - The published plug-in likelihood plugs in p_j = (x1+x2)/(n1+n2) from the same two counts it
    is fitting. It underestimates ρ by about half: relative bias −46 % to −69 %.
  - Intervals built from that estimate cover only 0.90–0.93.
  - It is kept for comparison (`Plugin_Likelihood_Rho`) but is not used for intervals.
- **Estimator used:** a moment estimator built from two expectation identities, which were
  verified by exact enumeration in the tests. Its relative bias is within ±4.1 % for ρ ≥ 0.01
  (target ≤ 10 %). The site bootstrap uses 1,000 resamples and seed 20261115.
- **Statuses:** Estimated; Estimated, no extra-binomial variation; too few sites (< 10);
  unstable (> 5 % of resamples at ρ = 0.5); not a usable replicate pair (≠ 2 samples, different
  date, location, or reference).
- **Datasets without replicates:** no ρ is produced, and the panel states the required sentence.
  Sensitivity values are labelled "Hypothetical; not calibrated".

## Milestones 3–4 — Intervals and simulation study (completed)

**Interval:**

- It is a beta-binomial score interval, equivalent to Wilson at the effective depth
  n / [1 + (n−1)ρ] (Kish design effect).
- It is labelled "approximate (score)", not exact.
- Calibrated results use the **upper 95 % bootstrap bound of ρ**, so the uncertainty of the
  estimate is included.

**Study:** `validation/simulate_overdispersion.py`, seed 20261115, 90 s. Exact beta-binomial sums
were used wherever coverage given ρ could be computed exactly. Criteria C1–C4 were fixed above
before the run.

- C1, known ρ:
  - The score interval passes for every ρ: mean coverage 0.947–0.953, minimum 0.907–0.940.
  - Wilson fails for every ρ > 0; at ρ = 0.05 its minimum is 0.35.
  - The LR interval fails the minimum, with 0.83 near the boundaries, so it is **not promoted**.
- C2, Phase 1 scenario (same seed and identical draws, verified):
  - Wilson coverage at 300× is 0.553, 0.533, and 0.522.
  - The new interval gives 0.952, 0.951, and 0.948, within 4 Monte Carlo SE of its exact
    coverage.
- C3, ρ estimated from 5–50 sites:
  - The upper-bound rule is supported for every design: mean coverage 0.955–0.972 and minimum
    cell 0.924.
  - The point-estimate rule fails at 5 sites (0.924).
  - With 5 sites, the bootstrap upper bound is above the true ρ only 73–76 % of the time, so
    **10 sites are required**.
- Misspecification checks:
  - When ρ differs between sites (Gamma, CV 1), coverage is 0.96–0.97.
  - An assumed ρ that is too small undercovers (true 0.05, assumed 0.02, 300×: 0.81). One that is
    too large overcovers (≈ 1.0), so hypothetical values are not a calibration.
- Runtime:
  - 8,000 score intervals take 0.2 ms.
  - Calibration of 50 sites with 1,000 bootstrap resamples takes 13 ms.
  - The plug-in likelihood takes 0.2 s.

**Tests:** `tests/test_overdispersion.py` (19) and `tests/test_overdispersion_simulation.py` (8, fast
versions). One test expectation was corrected after inspection: Wilson at ρ = 0.02, 300× covers
0.85 at p = 0.01 and 0.99, not < 0.70. The assertion now states the real pattern; the code was not
changed.

## Milestone 5 — Website integration (completed)

- `genomics/overdispersion_analysis.py` (new) orchestrates the Phase 2 analysis. It is called from
  `run_depth_aware_analysis`.
- Section 9 starts with the panel **"Sequencing Uncertainty and Variant Frequency"**:
  - the required explanation sentence;
  - a calibrated or not-calibrated badge with a headline;
  - the model comparison table;
  - up to 6 example observations;
  - the reasons any result is unavailable.
- New figures: `uncertainty_models` and `replicate_concordance` (shown only when replicates exist).
  The time series overlays calibrated bars for replicated samples only. The three Phase 1 figures
  are kept.
- New setting: hypothetical ρ values (default "0.002, 0.02, 0.05"; empty means none; validated).
  It is added to the advanced settings, the toolbar re-run, and provenance.
- New status rows: "Overdispersion calibration (Sprint 5 Phase 2)" and "Overdispersion sensitivity
  (hypothetical)".
- Downloads:
  - `sprint5_vaf_observations.csv` gains 11 appended columns; the first 37 Phase 1 columns are
    unchanged.
  - New files: `sprint5_uncertainty_models.csv`, `sprint5_overdispersion_sensitivity.csv`,
    `sprint5_overdispersion_calibration.csv`, and `sprint5_calibration_sites.csv`.
- **Regression found and fixed:**
  - The first version made `mmc2`'s observation table 8.6 MB, over the 8 MB download limit, so the
    file disappeared.
  - Per-row labels were shortened; the full wording is kept on the page and in this report.
  - The sensitivity table no longer repeats Wilson columns.
  - Result: 6.27 MB and 3.86 MB.
- Synthetic demo: `examples/integrated/make_replicate_example.py`, which writes
  `synthetic_replicate_{snvs,base_counts}.csv` and `synthetic_replicate_reference.fasta`.
  - SPEC_A (true ρ 0.01) is estimated at 0.0105, and SPEC_B (true 0.002) at 0.0021.
  - The plug-in estimates are 0.0055 and 0.0009.
- Web tests: `tests/test_sprint5_phase2_web.py` (10), covering:
  - replicates;
  - no replicates;
  - user and empty settings;
  - an invalid setting;
  - missing depth;
  - replicates without dates;
  - Sprint 2+3 combined;
  - `mmc2`;
  - the loader.

## Milestone 6 — Regression and performance (completed)

- **Full suite, intermediate (system Python, no MAFFT on PATH):** Ran 535 tests in 301.9 s — OK
  (skipped=2). That is 475 existing tests plus 60 new ones.
- **Phase 2 tests in the project `.venv`** (no SciPy): Ran 60 — OK (skipped=2: the optional SciPy
  check, and mmc2 because no path was set).
- **Regression of established outputs (verified):**
  - `run_integrated_analysis.py --example` was rerun with the `.venv` interpreter and MAFFT, and
    compared with the Phase 1 output saved in the scratchpad (`s5/cli_out`).
  - Every Phase 1 CSV is identical. The exceptions are `analysis_status.csv` (+2 rows) and
    `provenance.csv` (+2 parameters).
  - In `sprint5_vaf_observations.csv` the first 37 columns are identical, with 11 appended.
  - All 23 Phase 1 figures are **byte-identical**.
  - New files: `sprint5_uncertainty_models.csv`, `sprint5_overdispersion_sensitivity.csv`, and
    `uncertainty_models.png`.
  - This covers MAFFT, Sprint 2 distances and screening, Sprint 3 QC, Yue–Clayton, PCoA,
    cross-source matching, Phase 1 VAF provenance, and missing-coverage logic. Sprint 4 code was
    not touched.
- **Performance** (`.venv`, guided web path, `mmc2.xlsx`):
  - One run peaked at 296 MB. The A/B run with Phase 2 disabled gave 291 and 309 MB, so there is
    no detectable regression.
  - Four runs peaked at 306 MB (Phase 1: 307).
  - With 512 MB and 10 % CPU: 285 MB in 110.4 s (Phase 1: 282 MB, 107.5 s).
- **Original data:** the `mmc2.xlsx` md5 is `c0a75b7fa5bde8db09746613354a9f64`, unchanged and read
  only.
- **Final full suite (MAFFT on PATH, final code):** Ran 535 tests in 297.1 s — OK (skipped=1, the
  optional NCBI files only).
- Late fixes before the final run:
  - the interval-width figure caption read "(; gray…" when nothing was calibrated;
  - the time-series caption now explains the calibrated bands;
  - a duplicate ρ curve is no longer drawn.

**Remaining work:** none for Phase 2. Render is unverified, because deploying was not allowed.
Recommended next steps are in §14 of the report: obtain real replicate libraries before any
mutation clustering.
