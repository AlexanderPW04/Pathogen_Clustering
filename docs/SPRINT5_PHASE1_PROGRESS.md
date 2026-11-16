# Sprint 5 Phase 1 — Progress Checklist

Persistent record for the Sprint 5 Phase 1 goal loop (depth-aware pathogen variant analysis).
Read this file first when resuming. Labels: **Verified** (run or read in this session),
**Inferred**, **Unverified**.

Scratchpad for diagnostics (not part of the project):
`C:\Users\ALEXAN~1\AppData\Local\Temp\claude\a7271a137b0bdb68\922255f3-0871-46d4-aedd-17d549d92b7a\scratchpad\s5`

## Milestone status

| Milestone | Status | Attempt |
|---|---|---|
| 0. Read audit, baseline | **Done** | 1 |
| 1. Repair Sprint 3 issues (1A–1F) | **Done** | 1 |
| 2. Sprint 5 variant data model | **Done** | 1 |
| 3. Depth-aware VAF + Wilson intervals | **Done** | 1 |
| 4. Mutation-frequency time series | **Done** | 1 |
| 5. Website integration | **Done** | 1 |
| 6. Validation, simulation, performance, regression | **Done** | 1 |
| Final report | **Done** (`docs/SPRINT5_PHASE1_IMPLEMENTATION_REPORT.md`) | 1 |

---

## Milestone 0 — Baseline (completed)

**Audit read:** `docs/SPRINT3_AUDIT_REPORT.md`, sections A–I and Appendices A–B (verified, read in full).

**Git state (verified, November 14, 2026):**

- Branch `sprint4-integration`, HEAD `b641533` ("Stabilize Sprint 3 Excel analyses and simplify the
  Sprint 4 page"), one commit ahead of `origin/sprint4-integration`. `main` = `origin/main` = `eb72013`.
- Working tree: only `?? docs/` (the audit report and the demonstration-guide PDF). No uncommitted
  application changes to preserve.
- The memory fix is in the tree: `observation_matrix` builds column lists, not one dict per row
  (`genomics/temporal_analysis.py:26-84`); compressed result storage and the analysis slot are in
  `app.py`; `MAX_OBSERVATION_CELLS` is enforced at review.

**Sprint 3 execution path (verified by reading code):**
`/review` (`app.review_upload`) → `input_inspection.review_dataset` → `/run` (`app.run_reviewed`) →
`input_inspection.prepare_inputs` → `app.run_analyses` → `app.run_genomic_analyses` →
`genomics.integrated_analysis.run_integrated_analysis` → `wastewater.processed_snv.load_processed_snvs`
→ `pooled_observations.pooled_variant_records` / `pooled_sample_records` →
`quality_control.apply_pooled_quality_control` → `temporal_analysis.observation_matrix` /
`frequency_trends` → `diversity_analysis.population_structure` → `integrated_view.build_dashboard`
/ `download_files` → `app.store_result_files`.

**Baseline tests (verified):**

```bash
FONTENELE_MMC2=<scratchpad>/data/mmc2.xlsx AI_ENRICHMENT_ENABLED=false PATH=/c/Tools/mafft-win:$PATH .venv/Scripts/python.exe -m unittest discover -s tests
```

`Ran 409 tests in 212.9s — OK (skipped=1)` (skip: `test_optional_real_ncbi_data`). Output:
`<scratchpad>/s5/baseline_tests.txt`.

**Baseline performance (verified, `repro_audit.py guided mmc2.xlsx`, Windows working set,
`OPENBLAS_NUM_THREADS=1`):** review 0.8 s (149 MB); run HTTP 200, page 1,899,584 bytes, 52 pooled
samples; workbook opened 4 times; total 10.6 s; **peak 286 MB**. Matches the audit (286–288 MB).
`mmc2.xlsx` md5 `c0a75b7fa5bde8db09746613354a9f64` (unchanged).

**New baseline observation (verified):** `mmc2` `ALT_FREQ` text has 1–6 significant digits
(6,858 rows with 6). `ALT_FREQ × TOTAL_DP` is within 0.0045 of an integer, not exactly integral:
the frequencies are rounded. With the column's 6-significant-digit precision, all 7,973 rows map to
one unique integer alternate count that reproduces the reported text. With each row's own text
precision, 204 rows (mostly `1`, `0.5`, `0.2` at high depth) would be ambiguous. The audit's
"whole number in every row" statement is therefore true only up to the rounding of the reported
frequency.

**Audit issues addressed in this phase:**

| Goal item | Audit ID | Code |
|---|---|---|
| 1A trends ignore non-detections | S3-02 | `temporal_analysis.frequency_trends` |
| 1B Yue–Clayton distance 0 without shared coverage | S3-03 | `diversity_analysis.paper_distance_matrix`, `population_structure` |
| 1C ambiguous dates | S3-05 | `processed_snv.load_processed_snvs` |
| 1D figure error without dates | S3-04 | `integrated_analysis.run_integrated_analysis`, `figures.trajectory_figure` |
| 1E pooled reference unverified | S3-07 | `quality_control.reference_check`, new reference status |
| 1F reported vs count frequency | S3-06 (VAF-1) | `processed_snv.load_processed_snvs` |
| 1F duplicates in a median | S3-13 (DUP-2) | `integrated_analysis.position_summary`, `temporal_analysis.location_date_summary` |

**Intended modules:** `wastewater/processed_snv.py`, `genomics/temporal_analysis.py`,
`genomics/diversity_analysis.py`, `genomics/integrated_analysis.py`, `genomics/quality_control.py`,
`genomics/pooled_observations.py`, `genomics/figures.py`, `integrated_view.py`, `dashboard_guide.py`,
`app.py` (date-order option only), `genomics/input_inspection.py` (date review), templates; new
`genomics/collection_dates.py`, `genomics/variant_evidence.py`, `genomics/vaf_trends.py`,
`genomics/depth_aware_analysis.py`; new tests.

**Scientific decisions recorded before coding:**

1. *Non-detection.* A value of 0 is used only when per-position A/C/G/T counts show zero alternative
   reads at a position with adequate depth. A depth-only table shows adequate coverage but not zero
   alternative reads (the caller may have a reporting threshold), so those positions are counted as
   "not reported at adequate depth (value below the caller's reporting limit, unknown)" and are not
   set to 0.
2. *Duplicates.* A pooled sample has one frequency per allele at a position. Repeated rows with the
   same sample, reference, position, Ref, and Alt are one observation; technical replicates must use
   separate sample IDs.
3. *Reconstruction.* Alternate counts are reconstructed only when no direct count exists and every
   frequency-only record in the dataset maps to a unique integer under the column's precision.

---

## Milestone 1 — Sprint 3 repairs (completed, attempt 1)

**Files changed:** `genomics/collection_dates.py` (new), `wastewater/processed_snv.py`,
`genomics/unified_schema.py`, `genomics/pooled_observations.py`, `genomics/temporal_analysis.py`,
`genomics/diversity_analysis.py`, `genomics/quality_control.py`, `genomics/integrated_analysis.py`,
`genomics/figures.py`, `genomics/input_inspection.py`, `integrated_view.py`, `dashboard_guide.py`,
`app.py`, `templates/index.html`, `templates/_toolbar.html`; tests `tests/test_sprint3_repairs.py`
(new, 22 tests), `tests/test_analyses.py` and `tests/test_integrated_example.py` (column rename only).

| Item | Repair | Test (independent expected value) |
|---|---|---|
| 1A S3-02 | `observation_matrix` adds `Callable_Frequency`/`Callable_Evidence`/`Alt_Read_Count`. Known values: reported frequency; A/C/G/T counts at adequately covered unreported positions (0 alt reads -> 0; some -> alt/depth). Depth-only non-reports stay unknown (not 0). `frequency_trends` uses callable values, counts each state, adds `Complete_Series` and `Trend_Basis`, and omits series with no known value at all (count in a note). | Slope 0.16/week with a base-count zero (hand OLS); 0.128571/week without counts; single point -> no trend |
| 1B S3-03 | `paper_distance_matrix` returns NaN below `minimum_shared_positions` (default 1); `population_structure` leaves out samples with no data, skips the PCoA when an undefined pair remains, and returns `yue_clayton_pairs.csv` (shared positions, per-position mean, reason). Jaccard of two empty sets is NaN. | No-overlap pair NaN (old 0); 0 and 1 hand values; PCoA skipped with pair named; sample with no reads left out |
| 1C S3-05 | `collection_dates.parse_collection_dates`: ISO, month names, unambiguous numerics; ambiguous, two-digit-year, and partial dates kept as text, flagged, NaT. Optional day/month order (review form "Pooled-sample collection dates"); contradictions are invalid. Original text and status kept in the variant table. Review step warns. | 14 parsing cases; selected order; loader flags and resolves |
| 1D S3-04 | Trajectory and cross-source timeline figures skip undated rows; with no usable dates the page shows "Temporal trends unavailable because collection dates were not provided…" and the status row says Skipped. | No figure errors; other figures still drawn |
| 1E S3-07 | `quality_control.reference_verification`: Verified compatible / with record-level mismatches / Verified mismatch / Supplied but unverified / No reference provided; label source (paper label is assigned by the app). Table `reference_verification.csv`, status row, quality note. | Four statuses; synthetic example 20 of 21 match (one deliberate mismatch) |
| 1F S3-06 | `add_frequency_checks`: `Count_Derived_Frequency`, `Frequency_Difference`, `Frequency_Check` (tolerance 0.01); reported value never changed. Precedence: Sprint 3 analyses use the reported frequency; Sprint 5 uses counts and withholds intervals on disagreement. | 0.30 vs 40/100 flagged; 0.333 vs 1/3 agrees |
| 1F S3-13 | `mark_repeated_records`: `Duplicate_Status`, `Counted_Once`; medians in `position_summary` and `location_date_summary` count each sample once. | Median 0.3 (old 0.4) |

**Tests run:** new file 22/22 OK; Sprint 3 files (71) OK after the documented rename; web/import
files (115) OK.

**Real workbook (direct call):** `mmc2.xlsx` 9.5 s, no figure errors; 3 identical repeats; dates
all parsed; reference "Supplied but unverified" (label assigned from the paper's methods).

**Expected output changes (documented, not test-fitting):**
- `temporal_trends.csv`: `Change_In_Reported_Frequency`, `First/Last_Reported_Frequency` replaced by
  callable columns; series with no known frequency omitted (mmc2: 130,387 -> fewer rows; the
  baseline file was 6.7 MB, the unfiltered new one 19.4 MB and would have exceeded the 8 MB download limit).
- `variant_positions.csv` and `site_date_summary.csv` medians count each sample once.
- New columns in `pooled_variants_with_qc.csv`; new tables `reference_verification.csv`,
  `yue_clayton_pairs.csv` (paper method only).
- Partial dates (year or year-month) are no longer read as the first of the month.

**Scientific assumptions:** a depth-only "adequate depth, not reported" is not 0 (caller threshold
unknown); the default minimum of 1 shared position only prevents an undefined distance and is
configurable; the 0.01 frequency tolerance absorbs rounding to two decimals.

---

## Milestones 2–5 — data model, VAF uncertainty, time series, website (completed, attempt 1)

**New modules:** `genomics/frequency_uncertainty.py` (Wilson interval, depth illustration),
`genomics/variant_evidence.py` (data model, count provenance, reconstruction gate, validation),
`genomics/vaf_trends.py` (three-state time series, vectorized summary, panel selection),
`genomics/depth_aware_analysis.py` (tables, figures, availability messages, summary).
**Changed:** `genomics/integrated_analysis.py` (one call after the Sprint 3 pooled analyses; two
status rows), `genomics/figures.py` (three figures), `integrated_view.py` (section 9 "Sprint 5 —
Depth-Aware Variant Analysis", finding line), `dashboard_guide.py` (three figure guides),
`templates/index.html` (visible availability list), `genomics/temporal_analysis.py` (one-byte
indicators; categorical comparisons).
**Tests:** `tests/test_variant_evidence.py` (20), `tests/test_sprint5_depth_aware.py` (20: 13 unit,
7+ Flask end-to-end including the published workbook) — all OK.

Key results (verified):
- Wilson implementation equals Newcombe (1998) Table I to 4 decimals and a bisection inversion of the
  score test to 1e-9. 4/20 -> 0.0807–0.4160 (width 0.335); 40/200 -> 0.1505–0.2609 (width 0.110).
- `mmc2.xlsx`: precision 6 significant digits; all 7,973 counts reconstructed uniquely; 7,970
  counted once with approximate intervals; 474 mutation-site series with 2+ dates (independently
  counted from the raw workbook in the test), 118 with 3+.
- Memory regression found and fixed: first version peaked at 336 MB (text copies of a 295,000-row
  categorical column and 64-bit indicator columns). After fixes: 287 MB, peak inside the existing
  Sprint 3 `observation_status_figure`.
- Web path: Sprint 5 section, availability messages, and all four Sprint 5 CSVs plus figures
  download for the synthetic example and `mmc2.xlsx`.

**Precedence (documented):** Sprint 3 tables and trends use the reported frequency; Sprint 5 uses
alternative reads / total depth when counts are valid and agree within 0.01, otherwise the reported
frequency without an interval.

---

## Milestone 6 — validation, simulation, performance, regression (attempt 1)

- **Simulation** (`validation/simulate_vaf_intervals.py`, seed 20261114, SYNTHETIC): bias <= 0.0005;
  exact Wilson coverage 0.911–0.983 (simulated 0.908–0.983; mean about 0.95); width at VAF 0.2: 0.326 / 0.215 / 0.154 / 0.090
  at 20 / 50 / 100 / 300x; beta-binomial (rho 0.02) coverage drops to about 0.52 at 300x;
  detections-only slopes attenuated (0.0319 vs 0.0333 rising, -0.0201 vs -0.0222 falling at 50x),
  supported non-detections remove the bias. Test `tests/test_vaf_simulation.py` (4).
- **Performance (mmc2, web path):** 287 MB single run; 307 MB after four runs (baseline 302);
  512 MB + 10% CPU: 107.5 s, 282 MB (baseline 131.0 s, 279 MB).
- **Scenario matrix (22 audit scenarios):** no crash; same accept/reject outcomes; +3 figures for
  pooled runs. Combined vs pooled-only: pooled and Sprint 5 tables byte-identical.
- **mmc2 + NC_045512.2:** Verified compatible, 7,973 / 7,973.
- **Full suite run 1:** 475 tests, 1 failure (`test_ui_navigation` pinned the 8-section list; the
  new section 9 is intended). Test updated to 9 sections with a comment.
- Section 9 checked in the built-in browser by reading the rendered page (pane hidden, so no
  screenshot): navigation lists section 9; messages, 7 statistics, 3 figures, 4 downloads.
- Original workbooks in Downloads unchanged (md5 c0a75b7fa5bde8db09746613354a9f64).
- **Full suite run 2 (final):** `Ran 475 tests in 305.6s — OK (skipped=1)`; all 3 real-workbook tests ran.
- **Status:** all mandatory criteria met (one goal-loop attempt). Not done, by instruction: commit, push,
  deploy, Render verification. Next action: user review.
