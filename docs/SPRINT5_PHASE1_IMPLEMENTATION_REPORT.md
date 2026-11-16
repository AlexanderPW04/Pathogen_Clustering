# Sprint 5 Phase 1 — Implementation Report

**Depth-Aware Pathogen Variant Analysis**

- **Branch:** `sprint4-integration`. The base is local commit `b641533`, the Sprint 3 memory fix.
- **State:** all changes are uncommitted, as the goal instructions required.
- **Date:** November 14, 2026.
- **Not done (as instructed):** nothing was pushed, merged, deployed, or committed, and no paid API request was made.
- **Unchanged inputs:** the original `mmc2.xlsx` is unchanged (md5 `c0a75b7fa5bde8db09746613354a9f64`).
- **Related documents:** the Sprint 3 audit is `docs/SPRINT3_AUDIT_REPORT.md`; the milestone checklist is `docs/SPRINT5_PHASE1_PROGRESS.md`.

**Labels used in this report:**

- **Verified:** run or read in this session.
- **Inferred:** a conclusion drawn from verified observations.
- **Unverified:** not tested.

---

## 1. Overview

Sprint 5 Phase 1 adds an optional section to the existing website: **"Sprint 5 — Depth-Aware Variant Analysis"**. It examines how strongly sequencing data support observed pathogen mutations and how their frequencies change over time.

It is built on the Sprint 3 pooled pipeline and does not duplicate it.

- **Inputs it reuses:** the quality-controlled pooled variant records, the per-position table, the reference check, and the observation matrix.
- **What it adds:**
  - a validated record per observation, keeping the source values and labelling every derived value;
  - a 95% Wilson read-sampling interval for each frequency that has valid read counts;
  - time series per mutation and site that keep three states apart: detected, supported non-detection, and unknown;
  - three figures, four downloadable tables, and plain statements of what the data cannot support.

Before the new work, the five scientific and reliability issues and the two numerical discrepancies from the Sprint 3 audit were repaired and tested (section 2).

Out of scope, as instructed:

- raw-read alignment or variant calling;
- MAGOS-style clustering;
- lineage assignment and functional annotation;
- new Claude or PubTator features.

The user-facing text uses no tumor or cancer terms.

## 2. Repairs completed from the Sprint 3 audit

| Audit item | Problem (audit) | Repair | Tests |
|---|---|---|---|
| S3-02 Trends | Slopes used detections only | **Known values** now feed trends: reported frequencies, plus base-count values at adequately covered unreported positions (0 alternative reads gives 0). A depth-only non-report is **not** set to 0, because the caller's reporting limit is unknown. Every state is counted. Series with no known value are omitted from `temporal_trends.csv` (the count is given in a note). | `test_sprint3_repairs.CallableTrendTest` (5) |
| S3-03 Yue–Clayton | Pairs with no shared coverage got distance 0 | A pair with fewer than `minimum_shared_positions` (default 1) shared covered positions gets **NaN**. Samples with no data are left out. If any undefined pair remains, the PCoA is skipped with the pairs named. `yue_clayton_pairs.csv` reports the shared positions, the positions per sample, a per-position mean, and the reason. Two empty Jaccard sets also give NaN. | `SharedCoverageDistanceTest` (5) |
| S3-05 Dates | `03/04/2021` was silently read month-first | New `genomics/collection_dates.py` (rules below). Ambiguous dates, two-digit years, and partial dates are kept as text, flagged, and not used. The user may choose day-first or month-first (Advanced settings). A date that contradicts the chosen order is invalid. The review step warns. | `CollectionDateTest` (5), web test |
| S3-04 Figures without dates | Raw matplotlib error on the page | Time figures drop undated rows. With no usable dates, the page says "Temporal trends unavailable because collection dates were not provided…" and the status row shows *Skipped*. | `MissingDateFigureTest` (2), web test |
| S3-07 Reference | Pooled-only runs trusted the labels | `reference_verification.csv` assigns one of five statuses per reference label: *Verified compatible*, *…with record-level mismatches*, *Verified mismatch*, *Supplied but unverified*, *No reference provided*. It also records where the label came from; the paper's NC_045512.2 label is assigned by the app. Pooled-only analyses still run, labelled unverified. | `ReferenceVerificationTest` (2) |
| S3-06 Frequency vs counts | The reported value was used and only a vague warning given | Per record: `Count_Derived_Frequency`, `Frequency_Difference`, and `Frequency_Check` (tolerance 0.01). The reported value is never changed. Precedence is documented in section 4.5. | `FrequencyAndDuplicateTest`, `test_variant_evidence` |
| S3-13 Duplicates | Medians counted repeated rows | `Duplicate_Status`/`Counted_Once` mark each repeat as identical or conflicting. Medians count each sample once, so the hand example's median moves from 0.4 to 0.3. Replicates must use separate sample IDs. | `FrequencyAndDuplicateTest` (3) |

**Date rules.**

- Year-first dates (`2021-04-03`, also with a spreadsheet time) are always read as year-month-day.
- Month names are unambiguous.
- A numeric day/month date is unambiguous when one part is above 12 or both parts are equal.
- Every other numeric day/month date is ambiguous.

## 3. New scientific methods

1. **Variant data model** (`genomics/variant_evidence.py`). The sources of read counts, in order of precedence:
   - *Supplied*: ALT_DP or `Alt_Depth`; a supplied REF_DP is also kept.
   - *Base counts*: the Ref and Alt counts in the A/C/G/T table, with that table's depth as the denominator.
   - *Reconstructed*: the reported frequency × total depth, accepted only under the rule in section 4.3.
   - *Frequency only* or *No depth*: no count is established.

   A reference count is **never** derived as depth − alternate reads.
2. **Depth-aware VAF with a 95% Wilson interval** (`genomics/frequency_uncertainty.py`).
3. **Three-state mutation time series per site** (`genomics/vaf_trends.py`), with intervals, depth range, and interpretation notes.
4. **Orchestration and availability messages** (`genomics/depth_aware_analysis.py`).

## 4. Formulas and assumptions

### 4.1 VAF

For count-based records, **VAF = x / n**:

- *x* is the number of alternative-allele reads;
- *n* is the total usable depth: the denominator of the source frequency, which includes any other alleles at the position.

**Checks applied:**

- n ≥ 1;
- x is a whole number with 0 ≤ x ≤ n;
- frequencies are fractions from 0 to 1 (percentages are rejected at import, as before);
- Ref and Alt are each one of A, C, G, T and are different;
- a missing depth stays missing and never becomes 0.

Depth 0, negative counts, non-whole counts, and counts above the depth are rejected at import.

**Multiallelic sites.** Each allele keeps its own binomial (that allele versus all other reads). If the counted frequencies at one sample and position sum above 1 + 0.01, no interval is made for those alleles.

### 4.2 Uncertainty: the 95% Wilson score interval

With p̂ = x/n and z = 1.959964:

```
centre = (p̂ + z²/2n) / (1 + z²/n)
half   = z/(1 + z²/n) · sqrt(p̂(1−p̂)/n + z²/4n²)
[lower, upper] = [centre − half, centre + half], with 0 exactly when x = 0 and 1 exactly when x = n
```

**Assumption:** reads are independent Bernoulli draws (a binomial model).

**What the interval does not include:**

- sequencing or PCR error;
- amplicon dropout or primer bias;
- mapping errors;
- RNA degradation;
- differences between wastewater samples.

It is **not** a calibrated measure of all wastewater uncertainty, and it is not the MAGOS error model. When counts are reconstructed, the interval is labelled **approximate**. With no valid counts, no interval is made.

### 4.3 Count reconstruction rule

The table's precision P is the largest number of significant digits written in any of its frequencies. Writers such as iVar print 6 significant digits and drop trailing zeros, so `0.14` means 0.140000.

For each frequency-only record:

1. Compute A = round(f·n).
2. Accept A if A/n, rounded to P significant digits, equals f, and neither A−1 nor A+1 does.
3. If a record has no reproducing A at all, the frequency and the depth may not share a denominator, so **no record in the table is reconstructed**.
4. A record whose A is not unique stays frequency-only.

For `mmc2.xlsx` (verified), f·n is not an exact integer: it lies within 0.0045 of one because the frequencies are rounded. At P = 6, all 7,973 records reconstruct uniquely. At each record's own written precision, 204 would be ambiguous. This is why the table-level precision is used, and the rule is documented.

**Arithmetic consistency does not prove that counts were measured.** Reconstructed counts are therefore labelled "derived, not measured" everywhere.

### 4.4 Time series

- **Identifier:** Reference_ID (one per genome segment) + position + Ref + Alt. Different references are never merged.
- **Per mutation and site:** samples, detections, supported non-detections, uncalled base-count reads, unknown samples, distinct dated time points, and the earliest and latest VAF with interval and depth.
- **Absolute change:** calculated with ≥ 2 distinct dates.
- **OLS slope per 7 days:** calculated with ≥ 3 distinct dates, from grouped sums equal to `numpy.polyfit`.
- **Notes per series:** depth spread (≥ 5-fold is flagged) and the share of unknown samples.
- These summaries are descriptive only. No significance test is made.

### 4.5 Precedence by analysis

| Analysis | Frequency used |
|---|---|
| Sprint 3 QC, matrix, trends, medians, Jaccard | Reported frequency (paper reproduction) |
| Sprint 5 VAF and time series | x/n when counts are valid and agree with the reported value within 0.01; otherwise the reported frequency, with no interval |
| Disagreement > 0.01 | Both values are kept; the interval is withheld; the record is flagged |

## 5. Changed files and major functions

**New files:**

- `genomics/collection_dates.py`: `parse_one`, `parse_collection_dates`, `date_warnings`.
- `genomics/frequency_uncertainty.py`: `wilson_interval`, `expected_interval_width`, `depth_dependence_table`.
- `genomics/variant_evidence.py`: `reconstruction_check`, `base_count_lookup`, `build_variant_evidence`, `validate_evidence`.
- `genomics/vaf_trends.py`: `known_points`, `series_summary`, `choose_series`, `unknown_dates`.
- `genomics/depth_aware_analysis.py`: `availability_messages`, `missing_information`, `run_depth_aware_analysis`.
- `validation/simulate_vaf_intervals.py`: synthetic validation.
- Tests: `tests/test_sprint3_repairs.py`, `tests/test_variant_evidence.py`, `tests/test_sprint5_depth_aware.py`, `tests/test_vaf_simulation.py`.

**Modified files:**

- `wastewater/processed_snv.py`: new date parsing; `REF_DP`→`Ref_Depth`; `significant_digits`; columns `Collection_Date_Original`, `Collection_Date_Status`, `Frequency_Significant_Digits`; clearer warnings.
- `genomics/unified_schema.py`: `FREQUENCY_TOLERANCE`.
- `genomics/pooled_observations.py`: `mark_repeated_records`, `add_frequency_checks`, provenance passthrough.
- `genomics/temporal_analysis.py`: `add_callable_evidence`, a rewritten `frequency_trends` with one-byte indicators, and a duplicate-aware `location_date_summary`.
- `genomics/diversity_analysis.py`: a guarded `paper_distance_matrix`, `pair_table`, and updates to `population_structure`.
- `genomics/quality_control.py`: `reference_verification`.
- `genomics/integrated_analysis.py`:
  - date guards and the `NO_DATES_MESSAGE`;
  - the reference table and the duplicate-aware `position_summary`;
  - the Sprint 5 call and new status rows.
- `genomics/figures.py`:
  - `trajectory_figure` and `frequency_change_figure` now plot known values;
  - new figures `vaf_time_series_figure`, `vaf_depth_figure`, `depth_dependence_figure`.
- `genomics/input_inspection.py`: `pooled_date_review`; the date order is passed to combined tables.
- `integrated_view.py`: section 9 and the reference and pair tables.
- `dashboard_guide.py`: three figure guides; updated wording on trends.
- `app.py`: the `pooled_date_order` setting (validated and passed on).
- `templates/index.html`: the date-order setting and a visible list of availability messages.
- `templates/_toolbar.html`: keeps the date-order setting.
- `tests/test_analyses.py`, `tests/test_integrated_example.py`: one column rename each (values unchanged).

**Not changed:** the Sprint 2 code (MAFFT, distances, clustering, geocoding) and the Sprint 4 code (`literature/`).

## 6. Dataset compatibility

| Input | Sprint 5 behaviour (verified by tests or runs) |
|---|---|
| iVar TSV with ALT_DP (and REF_DP) | Supplied counts; exact intervals; reference counts kept |
| Standardized CSV/XLSX with `Alt_Depth` | Supplied counts; disagreeing records flagged, no interval |
| Frequency + depth only, reproducible | Reconstructed counts; approximate intervals (labelled) |
| Frequency + depth, not reproducible | Frequency only; no intervals; message says why |
| No depth | "Read depth is unavailable, so frequency uncertainty cannot be calculated." |
| A/C/G/T position table | Base counts; supported non-detections (0 reads) enter time series |
| Depth-only position table | Adequate-depth non-reports are counted as unknown value (not 0) |
| No dates / ambiguous dates | "Collection dates are unavailable…" / flagged; selectable order |
| No reference FASTA | "Reference identity is unverified…"; analyses still run |
| Published `mmc2.xlsx` | Paper format; 52 samples; all counts reconstructed (approximate); detections-only series |
| Segmented references | Each segment is its own Reference_ID and is never merged. Uploading a multi-record reference FASTA is still not supported (S3-10, unchanged). |

## 7. Numerical validation (verified)

- **Wilson interval against published values.** It equals Newcombe (1998, Table I) to 4 decimals:
  - 81/263: 0.2553–0.3662
  - 15/148: 0.0624–0.1605
  - 0/20: 0–0.1611
  - 1/29: 0.0061–0.1718
- **Wilson interval against an independent method.** It equals a bisection inversion of the score test to 1e-9 for five more cases.
- **Depth dependence:**

| Reads | VAF | 95% interval | Width |
|---|---|---|---|
| 4/20 | 0.20 | 0.0807–0.4160 | 0.335 |
| 10/50 | 0.20 | 0.1124–0.3304 | 0.218 |
| 20/100 | 0.20 | 0.1334–0.2888 | 0.155 |
| 60/300 | 0.20 | 0.1587–0.2489 | 0.090 |
| 200/1000 | 0.20 | 0.1764–0.2259 | 0.050 |

- **Boundaries:** x = 0 gives [0, z²/(n+z²)]; x = n gives [n/(n+z²), 1].
- **Trend slopes:** hand OLS gives 0.16/week with a base-count zero and 0.128571/week without it. The synthetic South 15 T>C series gives 0.05/week, with a Wilson upper bound of z²/(120+z²) for 0/120.
- **Duplicate-aware median:** 0.3 (was 0.4).
- **`mmc2` reconstruction:** all 7,973 counts unique; |f·n − x| < 0.5 in every row.
- **Reference check:** `mmc2` with NC_045512.2 gives *Verified compatible*, 7,973 / 7,973 matching (agrees with the audit).
- **Number of eligible `mmc2` series:** 474 mutation–site pairs with two or more dates. This equals an independent count from the raw workbook in the test.

### Synthetic simulation (`validation/simulate_vaf_intervals.py`, seed 20261114)

All of these are **SYNTHETIC**. They check the interval under its own model; they do not validate real samples or lineage detection.

- **Bias:** |bias| of x/n ≤ 0.0005 at depths 20, 50, 100, and 300.
- **Coverage:** 0.911 to 0.983 by exact binomial sums; the simulation gives 0.908 to 0.983, in agreement within Monte Carlo error. The mean across the grid is about 0.95. This is the Wilson interval's known oscillation, lowest near VAF 0.01 or 0.99 at 50×.
- **Width at VAF 0.2:** 0.326 (20×), 0.215 (50×), 0.154 (100×), 0.090 (300×).
- **Overdispersion (beta-binomial, ρ = 0.02):** coverage falls to 0.52–0.53 at 300×. **Important limitation:** with amplicon-like extra variation, deeper sequencing does not make the binomial interval trustworthy.
- **Missing coverage:** a caller that reports only frequencies ≥ 3%, with linear true trends:

| Scenario | Depth | True slope / week | Detections only | With supported non-detections |
|---|---|---|---|---|
| Rising | 50× | 0.0333 | 0.0319 | 0.0333 |
| Rising | 300× | 0.0333 | 0.0331 | 0.0333 |
| Falling | 50× | −0.0222 | −0.0201 | −0.0222 |
| Falling | 300× | −0.0222 | −0.0219 | −0.0222 |

Using supported non-detections removes the bias; detections only flatten the slope.

## 8. Automated tests

**Full suite (verified):** `Ran 475 tests in 305.6s — OK (skipped=1)`. The suite was 409 tests at baseline; 66 were added.

- **Skipped:** `test_optional_real_ncbi_data`, as at the baseline.
- **Real-workbook tests (all ran and passed with `FONTENELE_MMC2` set):**
  - `test_optional_fontenele_table`;
  - `test_variant_evidence…test_published_workbook_counts_are_reconstructed`;
  - `test_sprint5_depth_aware…test_published_workbook`.
- **Environment:** MAFFT v7 on PATH, Python 3.11.9, Windows 11.

```bash
FONTENELE_MMC2=<path to mmc2.xlsx> AI_ENRICHMENT_ENABLED=false PATH=/c/Tools/mafft-win:$PATH .venv/Scripts/python.exe -m unittest discover -s tests
```

**First full run:** one failure. `test_ui_navigation` pinned the list of 8 dashboard sections, and the new section 9 is intended. The test now expects 9 sections, with a comment.

**New tests:**

| File | Tests | What they cover |
|---|---|---|
| `test_sprint3_repairs` | 22 | Repairs 1A–1F |
| `test_variant_evidence` | 20 | Wilson, reconstruction, direct/reconstructed/mismatch/missing/zero/invalid/multiallelic/duplicate/base-count records, real workbook |
| `test_sprint5_depth_aware` | 20 | Time series, availability, dashboard, and Flask end to end (details below) |
| `test_vaf_simulation` | 4 | Synthetic properties |

The Flask end-to-end tests in `test_sprint5_depth_aware` cover:

- Sprint 3 alone, and Sprint 2 + 3 combined;
- direct counts with no dates;
- a sparse frequency-only table;
- coverage gaps;
- ambiguous dates resolved by the user's choice;
- invalid and unsupported inputs;
- other routes;
- the published workbook.

**Documented expected-output changes** (bug fixes, not test-fitting):

- `temporal_trends.csv` columns changed: `Change_In_Reported_Frequency` became `Change_In_Callable_Frequency`, and the same numbers are asserted for the unchanged cases.
- Medians are duplicate-aware.
- Partial dates are no longer read as the first of the month.

## 9. Real workbook execution (verified, web path `/review` → `/run`)

- `mmc2.xlsx` returns HTTP 200 with 52 pooled samples and 7,973 records (7,970 counted once; 3 identical repeats).
- No figure error.
- Sprint 5:
  - 7,970 of 7,970 observations have an approximate interval;
  - 474 series have two or more dated values, 118 have three or more, and 6,770 are single time points;
  - all four Sprint 5 CSVs download.
- `temporal_trends.csv` (1.56 MB) and `sprint5_vaf_observations.csv` (5.08 MB) are under the 8 MB download limit. `sprint5_temporal_summary.csv` is 2.47 MB.
- The 22-scenario upload matrix from the audit: no crash, the same accept and reject outcomes as before, and 3 more figures wherever pooled data run.
- With and without Sprint 2 data, the pooled and Sprint 5 tables are byte-identical. Only `variant_positions.csv` gains individual carriers, as intended.

## 10. Performance and memory (local Windows working set, `OPENBLAS_NUM_THREADS=1`)

| Run (`mmc2.xlsx`, guided web path) | Audit baseline (`b641533`) | This phase |
|---|---|---|
| One run, peak memory | 286 MB (10.6 s, re-measured at start) | **287 MB** (10.7 s); highest point is the existing Sprint 3 `observation_status_figure` |
| Four runs in one process | 302 MB | **307 MB** (timing not comparable: the test suite ran at the same time) |
| 512 MB cap and 10% of one CPU core | 279 MB, 131.0 s | **282 MB, 107.5 s** (job committed 284 MB) |
| Workbook opened per run | 4 | 4 (unchanged) |
| Results page | 1.90 MB | 2.50 MB (three more figures and section 9) |
| Samples / observations / series with 2+ dates | 52 / 7,973 / — | 52 / 7,973 (7,970 counted once) / 474 |

**Regression found and fixed.** The first version peaked at 336 MB. The causes were text conversion of a 295,000-row categorical column and 64-bit indicator columns. Both were removed, and the peak is back at the baseline.

**Not measured:** Linux and Render memory, and Render's proxy timeout (unverified).

## 11. Known limitations

- The intervals cover binomial read sampling only. Under realistic overdispersion they are too narrow (section 7).
- `mmc2.xlsx` has no measured read counts, so its intervals are approximate. It has no per-position coverage, so its time series use detections only, and unknown samples are counted but not valued.
- A depth-only coverage table cannot show 0 alternative reads, so a non-report there is "value unknown", not 0.
- Earliest and latest values and slopes are descriptive. With unknown samples left out, a series may over-represent detections; each series carries a note saying so.
- Strand bias, base quality, and mapping quality cannot be checked from processed tables.
- Multi-record (segmented) reference FASTA upload remains unsupported (S3-10).
- The default minimum of one shared position for Yue–Clayton only prevents an undefined distance. Summed distances still grow with shared coverage; per-position means are offered for comparison.
- The time-series figure shows six series, chosen by the number of dated values. All series are in the CSV.
- No lineage, transmission, outbreak, or clustering claims are made. Similar frequencies do not show that mutations are linked on one genome.

## 12. How to use Sprint 5 locally

1. Start the app as described in the README, for example with `.venv/Scripts/python.exe app.py`, then open the local address.
2. Choose **Pooled sequencing** (or **Combined**) and upload a processed SNV table, for example `mmc2.xlsx` or `examples/integrated/synthetic_pooled_snvs.csv`. Optional files:
   - an A/C/G/T position table;
   - a reference FASTA.
3. In Step 2, read the review. If it warns about ambiguous dates, open **Advanced settings → Pooled-sample collection dates** and choose the order, but only if you know it.
4. Run the analysis and open section **9. Sprint 5 — Depth-Aware Variant Analysis**. The section shows:
   - the statements of what the data allow;
   - the summary statistics;
   - three figures: time series, depth and uncertainty, and the same-VAF illustration;
   - tables and downloads.
5. To run the synthetic validation:

```bash
.venv/Scripts/python.exe validation/simulate_vaf_intervals.py
```

## 13. Testing on Render (after review; not done here)

1. After you approve, commit and push the branch Render deploys, then rebuild. No new packages are needed: Sprint 5 uses numpy, pandas, and matplotlib only.
2. Upload `mmc2.xlsx` twice. In Render **Metrics**, check that the peak stays well below 512 MB (about 290–310 MB locally). In **Events**, check there is no "Ran out of memory".
3. Confirm that the request finishes within the proxy timeout. Locally it took 107.5 s under a 512 MB / 10% CPU cap; Render's proxy timeout was not measured.
4. Check that section 9 appears and that `sprint5_vaf_observations.csv` downloads.
5. Upload the synthetic pooled CSV with the A/C/G/T table, and confirm the open (non-detection) points appear in the time-series figure.

## 14. Recommended next development stage

1. **Overdispersion first.** Estimate the extra variation from technical replicates or from the spread between neighbouring samples. Replace or supplement the Wilson interval with a beta-binomial interval, and validate it by simulation, as section 7 shows is needed.
2. **Offline regeneration** of per-position A/C/G/T tables for the paper's samples (SRA PRJNA662596). This enables supported non-detections and the paper's Yue–Clayton PCoA.
3. **Only then:** an exploratory, depth-aware grouping of co-varying mutations (MAGOS-inspired, without tumor assumptions), validated on simulated mixtures of known genomes (audit section G.1).
4. **In parallel:** GFF-based gene and amino-acid annotation, replacing the published Sy/Ns labels (S3-18).

---

## Acceptance criteria

| # | Criterion | Met | Evidence |
|---|---|---|---|
| 1 | Audit read and referenced | Yes | Section 2; progress file Milestone 0 |
| 2 | Five Sprint 3 issues repaired or safeguarded | Yes | Section 2; `test_sprint3_repairs` |
| 3 | Frequency and duplicate discrepancies handled and tested | Yes | Sections 2 and 4.5 |
| 4 | Memory fix intact | Yes | 287 MB peak (baseline 286); 282 MB under a 512 MB cap |
| 5 | Validated Sprint 5 data model | Yes | `variant_evidence.py`, `validate_evidence`, 20 tests |
| 6 | Original values kept; observed vs reconstructed counts kept apart | Yes | `Alt_Count_Source`, `Reported_Frequency` |
| 7 | Uncertainty numerically validated | Yes | Newcombe values; bisection; exact coverage |
| 8 | Missing coverage never shown as 0 VAF | Yes | Tests for unknown, gap, and depth-only cases |
| 9 | Temporal summaries work with suitable data | Yes | Synthetic hand checks; `mmc2` 474 series |
| 10 | Results and downloads in the website | Yes | Section 9; Flask tests; rendered-page check |
| 11 | `mmc2.xlsx` completes without a crash | Yes | Web path, HTTP 200, no figure errors |
| 12 | Existing functionality operational | Yes | 475 tests; 22 scenarios; CLI runner |
| 13 | All new tests pass | Yes | 66 new tests |
| 14 | No unexplained regression failures | Yes | One expected UI test change, explained |
| 15 | Performance recorded | Yes | Section 10 |
| 16 | No unsupported biological claims | Yes | Section 11; dashboard limitations |
| 17 | No paid calls, pushes, merges, or deployments | Yes | AI enrichment disabled in every run; nothing committed or pushed |
| 18 | Final report exists | Yes | This file |

Behaviour on Render itself is **unverified**: its memory, its proxy timeout, and the deployed branch were not checked. Section 13 lists the checks to run after your review.

**SPRINT 5 PHASE 1 COMPLETE**
