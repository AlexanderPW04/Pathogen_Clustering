# Sprint 3 Audit Report — Pooled Sequencing Pipeline

**Scope:** the pooled (wastewater) sequencing pipeline on branch `sprint4-integration`, local commit `b641533` (one commit ahead of `origin/sprint4-integration`; `origin/main` is at `eb72013`).
**Date:** November 14, 2026.
**Type:** audit only. No application file was changed. No paid API call was made. All diagnostic inputs, outputs, and scripts are in the session scratchpad (listed in Appendix B). The original `mmc2.xlsx` was not modified (md5 `c0a75b7fa5bde8db09746613354a9f64`, unchanged).

Labels used throughout:

- **Verified:** observed in this audit by running code or reading it.
- **Inferred:** a conclusion drawn from verified observations, not directly observed.
- **Unverified:** an assumption that was not tested.

---

## A. Executive summary

**What Sprint 3 does.** It takes a processed pooled-sequencing table: one row per sample and single-nucleotide variant (SNV), with reference position, reference and alternative base, alternative-allele frequency, and read depth. Optional inputs are a reference FASTA, per-position depth or A/C/G/T counts, pooled metadata, a masked-site list, and variant signatures. From these it:

- applies the Fontenele et al. (2021) quality filters only where the supplied data allow them;
- builds a sample-by-variant observation table that separates "not reported" from "0%";
- summarizes reported frequencies over time and by site;
- draws a population-structure plot (PCoA). It uses the paper's Yue–Clayton method when A/C/G/T counts exist. Otherwise it uses a clearly labelled exploratory Jaccard ordination;
- compares pooled and individual-sequence variants on one verified reference.

**What works (verified):**

- The published workbook `mmc2.xlsx` imports completely: 7,973 records, 52 samples, 10 states. It runs end-to-end through the web upload path with no error in the current code.
- All 7,973 published reference bases match the RefSeq NC_045512.2 sequence, so the reference the importer assigns is correct for this workbook.
- 44 of 46 hand-calculated numerical checks agree with the program. The 2 differences are a deliberate design choice and a minor inconsistency (section D).
- The full test suite passes: 409 tests, 1 skipped, 242 s.
- Adding Sprint 2 individual sequences does not change any Sprint 3 output: 10 of 11 pooled tables are byte-identical. The 11th only gains individual-carrier columns.

**What is broken or weak (verified):**

- No open crash on valid input.
- **Bug:** with no collection dates, the "Variant frequencies over time" figure fails and its raw matplotlib error text is shown on the page.
- **Scientific gaps (none are crashes):**
  - Frequency trends use only *reported* frequencies, so adequate-depth non-detections never lower a slope.
  - The Yue–Clayton sum lets samples with no shared covered positions appear identical (distance 0).
  - Ambiguous day/month dates are read month-first without a warning.
  - Reported frequencies are trusted even when read counts disagree; only a warning is given.
- No VAF uncertainty or depth weighting exists anywhere. This is expected for Sprint 3, but it is the core of Sprint 5.

**The Excel crash.**

- *Reproduced on the pre-fix code (commit `6e122ff`):* it was resource exhaustion, not a parsing error.
  - Locally, the published workbook peaked at 499–500 MB per run. Repeated runs in one process reached 537 MB.
  - Under a 440 MB cap, figures were silently dropped.
  - On a 600,000-cell table under a 512 MB cap, `observation_matrix` raised `MemoryError`. The user saw the misleading message "Check the reference and pooled-sample files."
- *Current code:* about 286 MB per run, 302 MB after four runs, and 131 s at 10% of one CPU core under a 512 MB cap.
- *Not confirmed on Render:* no Render logs were available, so whether Render's out-of-memory killer ended the original hosted run remains **unverified**.

**Scientific validation.** VAF handling, depth QC, the Yue–Clayton formula, PCoA, Jaccard, the depth regression, temporal slopes, and cross-source categories were checked against independent hand calculations (section D). The formulas are correct. The concerns are about interpretation and missing guards, not arithmetic.

**Sprint 5 readiness: READY AFTER SPECIFIC FIXES** (section I).

- The published table supports per-sample VAF and total depth. Alternate read counts can be recovered exactly.
- It does **not** contain reference counts, per-position coverage, or A/C/G/T counts. Depth-aware absence and the paper's PCoA therefore need offline-regenerated tables.

---

## B. Execution architecture (verified call chain)

### B.1 Main guided path (`/` → `/review` → `/run`)

| Stage | Function (file:line) | What happens | Notes |
|---|---|---|---|
| Upload | `review_upload` ([app.py:1306](../app.py#L1306)) | Files are staged in a temp folder by `stage_uploads` ([app.py:1202](../app.py#L1202)); `MAX_CONTENT_LENGTH` is 10 MB ([app.py:70](../app.py#L70)) | Staged for 30 min, at most 6 sets ([app.py:208-209](../app.py#L208)) |
| File inspection | `inspect_file` ([input_inspection.py:316](../genomics/input_inspection.py#L316)) → `inspect_table` ([:201](../genomics/input_inspection.py#L201)) | Detects the format (paper table, iVar, standard, VCF), IDs, reference labels, variant count, 3-row preview | |
| Worksheet choice | `choose_sheet` ([:184](../genomics/input_inspection.py#L184)) → `describe_sheets` ([excel_tables.py:73](../genomics/excel_tables.py#L73)) | Scores the first 15 rows of each visible sheet; asks the user when several sheets look like data | Workbook opened `read_only=True, data_only=True` ([excel_tables.py:69-70](../genomics/excel_tables.py#L69)) |
| Header detection and parse | `read_sheet` ([excel_tables.py:119](../genomics/excel_tables.py#L119)) via `processed_snv.read_table` ([processed_snv.py:78](../wastewater/processed_snv.py#L78)) | Single streaming pass; every cell as text; header row is the row with most recognized names; Fontenele table recognized from its title row ([processed_snv.py:95-100](../wastewater/processed_snv.py#L95)) | |
| Dataset review | `review_dataset` ([input_inspection.py:499](../genomics/input_inspection.py#L499)), `check_analysis_size` ([:750](../genomics/input_inspection.py#L750)), `compatibility_review` ([:816](../genomics/input_inspection.py#L816)) | Role assignment, column mapping, value warnings, metadata readiness, reference-label compatibility, size limit (`MAX_OBSERVATION_CELLS`) | Re-runs on `/run` |
| Run | `run_reviewed` ([app.py:1592](../app.py#L1592)) → `prepare_inputs` ([input_inspection.py:1044](../genomics/input_inspection.py#L1044)) → `run_analyses` ([app.py:832](../app.py#L832)) → `run_genomic_analyses` ([app.py:908](../app.py#L908)) | Writes standardized inputs; one analysis at a time (`ANALYSIS_SLOT`, [app.py:195](../app.py#L195)) | `MemoryError` → `OUT_OF_MEMORY` message; other errors logged ([app.py:851-857](../app.py#L851), [953-963](../app.py#L953)) |
| Loading and validation | `load_processed_snvs` ([processed_snv.py:105](../wastewater/processed_snv.py#L105)) | Renames columns; validates bases, positions, frequency in [0,1], depth; sets aside 0% rows; warns on duplicates and count/frequency disagreement | Paper format gets reference `NC_045512.2` ([:135](../wastewater/processed_snv.py#L135)) |
| Shared schema | `pooled_variant_records` ([pooled_observations.py:200](../genomics/pooled_observations.py#L200)), `pooled_sample_records` ([:294](../genomics/pooled_observations.py#L294)) | One variant row per record; one sample row per sample, with coverage from metadata or a complete depth table | |
| QC | `apply_pooled_quality_control` ([quality_control.py:110](../genomics/quality_control.py#L110)) | Depth ≥ 20, Fisher p < 0.05, sample coverage > 90%, masked sites, Ref-base agreement; "Not evaluable" when values are missing | Paper mode excludes failures; exploratory mode flags them |
| Analyses | `run_integrated_analysis` ([integrated_analysis.py:160](../genomics/integrated_analysis.py#L160)) | Per-sample counts, depth regression, position summary, `observation_matrix` ([temporal_analysis.py:26](../genomics/temporal_analysis.py#L26)), `frequency_trends` ([:99](../genomics/temporal_analysis.py#L99)), `population_structure` ([diversity_analysis.py:159](../genomics/diversity_analysis.py#L159)), `compare_evidence` ([variant_matching.py:122](../genomics/variant_matching.py#L122)) | Each figure goes through `add_figure` ([integrated_analysis.py:147](../genomics/integrated_analysis.py#L147)); a drawing error becomes a note |
| Results | `build_dashboard` ([integrated_view.py:197](../integrated_view.py#L197)), `download_files` ([:454](../integrated_view.py#L454)), `store_result_files` ([app.py:618](../app.py#L618)) | CSV and PNG kept zlib-compressed in memory for 10 min, at most 4 result sets | Tables larger than 8 MB are skipped and listed |
| Download | `result_download` ([app.py:638](../app.py#L638)) | Unguessable 16-byte token URL | |

### B.2 Legacy and separate routes (verified working on `mmc2.xlsx`)

- `/wastewater` (`wastewater_page`, [app.py:445](../app.py#L445)): import and summary only (`summarize_snvs`). 7,973 records, 0.6 s, peak 151 MB.
- `/wastewater-data` and `/environmental` (`environmental_page`, [app.py:326](../app.py#L326)): for measurement tables, not SNVs. Given `mmc2.xlsx`, it correctly answers "This file looks like variant calls, not measurements."
- `run_integrated_analysis.py`: command-line runner of the same `run_integrated_analysis`.
- `offline_workflow/`: template for turning raw reads into tables offline. **Not run** in this audit; the tools are not installed.

### B.3 Pipeline behaviour checks

| Question | Finding | Evidence |
|---|---|---|
| Stages executed more than once | The workbook is opened **4 times** per guided run: sheet description at review, a read at review, sheet description inside `read_sheet` at run, and a read at run. Each pass is about 0.3 s; peak memory is unaffected. `review_dataset` also re-reads mapped tables on `/run`. | Instrumented `openpyxl.load_workbook` (Appendix B, `repro_audit.py`) |
| Large objects in session | No Flask cookie session. Staged files stay on disk for 30 min; result files stay in process memory (compressed) for 10 min. | [app.py:194-209](../app.py#L194) |
| Columns discarded | Only `Published_*`, `Source_*`, and `Alt_Depth` extras are kept ([processed_snv.py:330-340](../wastewater/processed_snv.py#L330)). Any other user column is dropped. | Code read |
| Silent data changes | Reference bases are uppercased. 0% rows are set aside *with* a warning. Duplicates are kept with a warning, and analyses use the first row. Ambiguous dates are parsed month-first **without** a warning (issue S3-05). | Section D; `audit_numeric_checks.py` |
| Upload mode inconsistency | Pooled-only and combined runs give byte-identical pooled tables. | `audit_combined_vs_pooled.py` |
| Failures without a clear message | Fixed: `MemoryError` and unexpected errors now have clear messages. Remaining: a figure-drawing failure shows raw library text (S3-04). | Scenario 12 |

---

## C. Root-cause analysis

### C.1 The `mmc2.xlsx` failure (reproduced on pre-fix code)

**Reproduction.**

1. Extract the pre-fix code with `git archive 6e122ff`.
2. Run `repro_audit.py guided mmc2.xlsx`, which drives the Flask test client through `/review` then `/run`.
3. Run it unconstrained, under a Windows Job Object memory cap, and four times in one process.

| Run (pre-fix `6e122ff`) | Result | Peak memory (process working set) |
|---|---|---|
| Unconstrained, 1 run | Completed (HTTP 200, 52 samples) | **499–500 MB** |
| 4 runs in one process | Completed | 499 → 512 → 524 → **537 MB** (about 37 MB retained per run) |
| 512 MB cap | Completed | 499 MB (job committed 503 MB) — 9 MB headroom |
| 440 MB cap | HTTP 200 but page **1.26 MB instead of 1.90 MB**: figures silently lost to `MemoryError` caught in `add_figure` | 454 MB |
| 600,000-cell generated table, 512 MB cap | `numpy._core._exceptions._ArrayMemoryError` in `observation_matrix` (`genomics/temporal_analysis.py:77` at `6e122ff`) ← `run_integrated_analysis:287` ← `app.run_genomic_analyses:887`. The user saw "The integrated analysis could not be completed. Check the reference and pooled-sample files." | 481 MB |

**Root cause (verified locally).**

- The workbook was read correctly. No parsing, header, numeric-conversion, or worksheet error occurred on any route.
- The analysis stage used almost all of a 512 MB instance:
  - `observation_matrix` built one Python dict per variant × sample: 5,677 × 52 ≈ 295,000 rows.
  - The heatmap converted the whole matrix again to draw 40 variants.
  - Figures and result files were retained between runs.
- Memory grew with each repeated run, and a single synchronous worker with a 180 s timeout ran the request at 0.1 CPU.

**What remains unverified.** That Render's out-of-memory killer, rather than its worker timeout, ended the original hosted run. That needs Render's **Events** ("Ran out of memory") and **Metrics** for the failed deploy. On Render, a killed worker usually appears to the browser as HTTP 502 (**inferred**). Linux memory use may differ from the Windows working-set figures above.

**Status in the current code (verified):**

| Run (current `b641533`) | Result |
|---|---|
| Unconstrained / 512 MB cap / 350 MB cap | Completed; peak **286–287 MB** (job 290 MB) |
| 4 runs in one process | 286 → 290 → 292 → **302 MB** (about 4–5 MB per run) |
| 512 MB cap and 10% of one core | Completed in **131.0 s** (review 9.9 s), peak 279 MB |
| 1.5 M-cell table with `MAX_OBSERVATION_CELLS=700000` (Dockerfile value) | Refused at review: "15,000 variants x 100 samples is 1,500,000 combinations, and this server can process 700,000. No samples or variants are removed automatically." |

**Smallest input that triggers the problem.**

- *Current code:* the review step refuses anything over 700,000 combinations. Measured peaks were 394 MB at 600,000 combinations and 665 MB at 1,500,000. The analysis limit is therefore enforced before memory runs out (verified).
- *Pre-fix code:* the published 295,000-cell workbook already sat within about 10 MB of 512 MB. The exact smallest failing size under Linux was not measured.

### C.2 Other reproduced failures (current code)

| ID | Reproduction | Error | Cause | Impact |
|---|---|---|---|---|
| F-1 | Upload `synthetic_pooled_snvs.csv` with every `Collection_Date` blank; run pooled mode | Results page note: "Variant frequencies over time could not be drawn: Failed to convert value(s) to axis units: array([NaT, NaT, NaT, NaT], dtype=object)" | `integrated_analysis.py:312-319` draws `trajectories` whenever a variant is reported twice, without checking dates; `figures.trajectory_figure` ([figures.py:220](../genomics/figures.py#L220)) plots NaT dates | One figure missing; technical message shown; the rest of the run is fine |
| F-2 | Same CSV with one `Alt` = `R`, or one blank `Sample_ID` | Review page shows **no** problem. The run then says "Processed SNV validation failed: 1 record(s) have invalid Alt; expected one A, C, G, or T (first at row 4: 'R')." | Nucleotide and blank-ID validation exists only in `load_processed_snvs`, not in the review step | User learns of the problem one step late; message is clear |

---

## D. Scientific validation

All "Verified = Yes" rows come from `audit_numeric_checks.py` (46 checks, 44 agree) or the scenario runs, unless stated otherwise.

| Method | Implementation | Verified? | Problems | Evidence |
|---|---|---|---|---|
| **VAF** | Reported `Alt_Frequency` is used as given ([processed_snv.py:246-250](../wastewater/processed_snv.py#L246), [305](../wastewater/processed_snv.py#L305)). It is never computed from counts. When `Alt_Depth` exists, `Alt_Depth/Read_Depth` is compared, and only a difference **> 0.05** produces a warning ([:296-302](../wastewater/processed_snv.py#L296)). | Yes (accepts reported VAF only) | (1) 0.30 reported vs 40/100 = 0.40 → 0.30 is used, with a warning (VAF-1). (2) A 0.03 disagreement passes silently (VAF-2). (3) Percentages > 1 are rejected (VAF-3), but "0.5 meaning 0.5%" cannot be detected (VAF-4). (4) Missing frequency is rejected, not set to 0 (VAF-5). (5) Two ALTs at one position may sum above 1 without a check (VAF-7); in `mmc2` none do (0 of 110 multi-ALT sites). (6) The denominator is whatever the source used: for iVar and the paper it is `TOTAL_DP` (all reads at the position), not ref + alt. | VAF-1…8; workbook inspection |
| **Read depth** | `Read_Depth` must be a whole number ≥ 1 or blank. Blank means "not evaluable" and the record is **kept** ([processed_snv.py:251-262](../wastewater/processed_snv.py#L251); [quality_control.py:65-70](../genomics/quality_control.py#L65)). `Alt_Depth` > depth is rejected. | Yes | Depth is used for QC and display only. It is never used to weight frequencies or estimate uncertainty. A table with all depths blank runs with no depth filtering, reported as "Not evaluable" (scenario 7). | DEPTH-1, DEPTH-2, scenario 07 |
| **QC** | `apply_pooled_quality_control` ([quality_control.py:110](../genomics/quality_control.py#L110)) with `PAPER_SETTINGS` ([:29](../genomics/quality_control.py#L29)): depth ≥ 20, Fisher p < 0.05 (only if supplied), coverage > 90% at 10× (only if supplied), masked sites (only if supplied), Ref vs reference (only with a FASTA). | Yes | Thresholds come from the paper's methods and are documented in code; they are not configurable from the UI except paper vs exploratory mode. Excluded records stay in the downloadable table with a reason. Base quality, mapping quality, and strand bias **cannot** be checked: processed tables do not carry them, and this is stated. For `mmc2`: 12 records < 20× excluded (7,961 kept); coverage, Fisher, and masking are not evaluable; all 52 samples kept. | QC-1…4; scenario 18 |
| **Temporal trends** | `observation_matrix` ([temporal_analysis.py:26](../genomics/temporal_analysis.py#L26)) marks every included sample × variant as observed, excluded by QC, adequate depth not reported, low or no coverage, or coverage unknown. `frequency_trends` ([:99](../genomics/temporal_analysis.py#L99)) gives first/last reported frequency and an OLS slope per week ([:87](../genomics/temporal_analysis.py#L87)) from **reported** values only, with ≥ 3 points. | Yes (arithmetic) | (1) **Selection bias:** an adequate-depth non-detection is counted (TIME-3) but ignored by the slope and change (TIME-4), so series can only be computed from detections. (2) No uncertainty, no depth weighting, and no test. (3) In `mmc2`, only 6 of 23 sites have ≥ 3 dates. (4) Dates come from the collection-date field, never the upload date. Ambiguous dd/mm values parse month-first silently (S3-05). (5) Wording says "not estimates of community prevalence", which is appropriate. | TIME-1…5; workbook inspection |
| **Yue–Clayton** | `yue_clayton_dissimilarity` ([diversity_analysis.py:35](../genomics/diversity_analysis.py#L35)): `1 − Σab / (Σa² + Σb² − Σab)` on A/C/G/T frequencies (count / sum of the four counts, [:47-66](../genomics/diversity_analysis.py#L47)). The pairwise distance is the **sum** over positions; positions without reads in either sample contribute 0 ([:69-87](../genomics/diversity_analysis.py#L69)). | Yes | The formula matches the hand calculations (0.5, 0, 1, a random pair). **However:** two samples with no covered position in common get distance 0, so they look identical (YC-7). `compared_positions` is computed but never shown or thresholded (S3-03). The summed scale also grows with coverage. Both follow the paper's stated method, but are risky without a guard. Requires A/C/G/T counts: never run on SNV-only tables (verified). | YC-1…8 |
| **PCoA** | `principal_coordinates` ([diversity_analysis.py:104](../genomics/diversity_analysis.py#L104)): classical Gower centering `B = −½ J D² J`, `eigh`, axes from positive eigenvalues; negative eigenvalues counted and reported, not corrected. Needs ≥ 3 samples on one reference. | Yes | Euclidean distances are recovered exactly. Eigenvalues equal squared singular values of the centered points. Explained fractions sum to 1. A non-metric matrix reports negative eigenvalues. For `mmc2` the method is **Exploratory Jaccard** on reported SNV sets, labelled "EXPLORATORY – not the paper's Yue–Clayton PCoA" on the figure. The figure guide's "measures" line still says "nucleotide frequencies" (S3-15). | PCOA-1…5, JAC-1…2 |
| **Cross-source integration** | `harmonize_pooled_variants` ([variant_matching.py:53](../genomics/variant_matching.py#L53)) and `compare_evidence` ([:122](../genomics/variant_matching.py#L122)): records are comparable only when the reference ID is identical or checksums match, the position lies within the reference, and the Ref base matches. Categories are at [:232-237](../genomics/variant_matching.py#L232). | Yes | All 7 hand-built cases match: shared, individual-only (adequate pooled depth), pooled-only (individual reference call), both unresolved directions, Ref-mismatch incomparable, and version-mismatch incomparable. "Individual only" requires at least one adequately covered pooled sample; QC-excluded calls count as unknown. | XS-* checks; scenarios 16, 17, 22 |
| Depth regression (paper Fig 2B) | `depth_regression` ([temporal_analysis.py:195](../genomics/temporal_analysis.py#L195)): OLS, bootstrap band, permutation p, seed 0 | Yes | Not available for `mmc2` (no mean depth) | REG-1 |

**Mathematical inconsistencies found:**

1. `position_summary` ([integrated_analysis.py:100](../genomics/integrated_analysis.py#L100)) takes the median over **all** rows, including duplicate sample/variant rows. Everywhere else uses the first row. Example: 0.2 and 0.6 for S1 plus 0.4 for S2 gives a median of 0.4 here, versus 0.3 if each sample counts once (DUP-2). This affects `variant_positions.csv`. In `mmc2` only 3 records are repeats.
2. No other inconsistency was found.

**Information present in `mmc2.xlsx` (verified by inspection):**

- One sheet, "All detected SNVs per sample": title row 1, headers row 2, 7,973 data rows, 20 columns, no merged cells.
- The `BIN` column holds 7,973 `COUNTIFS` formulas. Their cached values are read because of `data_only=True`.
- Sample IDs are mixed numbers and text (1,410 int, 6,563 str). `ALT_FREQ` and `TOTAL_DP` are mixed: 2,005 numeric and 5,968 stored as text. All parse correctly.
- REF bases are lowercase (uppercased on import).
- 29 collection dates from 2020-04-23 to 2020-06-16, one date per sample.
- Depth 12–12,317 (median 622). `ALT_FREQ` ranges from 0.00018 to 1.0, with 5,249 values below 0.05.
- `ALT_FREQ × TOTAL_DP` is a whole number in every row, so **alternate read counts are recoverable exactly**. Reference counts and other-allele counts are not.
- 3 exact duplicate rows: 2 true repeats, and 1 position annotated under both ORF7a and ORF7b. 110 sample positions carry two ALTs.
- No per-sample coverage, no Fisher p, no A/C/G/T counts, no position depth.
- Published annotation: every UTR row (144) is labelled "Syn". Its synonymous/non-synonymous labels are not reliable for non-coding regions.

`mmc3.xlsx` (Supplementary Table 3) has a different structure: 3 header rows with merged group headers; per-SNV novelty and clinical matches; no sample-level frequency or depth. It is rejected with a specific explanation (scenario 20), which is correct.

---

## E. Prioritized issue register

**No P0 issue is open** (verified): no crash on valid input, no invalid formula, no data corruption, no security exposure found. The original crash (P0 at the time) is fixed in `b641533` but **not deployed or pushed**.

| ID | Priority | File / function | Description | Evidence / reproduction | Impact | Recommended repair | Before Sprint 5? |
|---|---|---|---|---|---|---|---|
| S3-01 | **P1** | Deployment (`Dockerfile`, branch state) | The memory fix exists only in local commit `b641533` (ahead 1 of origin). It has not been tested on Render or in a rebuilt Docker image. Which branch Render deploys was not verified. | `git status -sb`; section C.1 | The hosted site may still run code that crashes on `mmc2.xlsx` | After review: push, build the image, run `mmc2.xlsx` twice on Render, read Metrics/Events | Yes |
| S3-02 | **P2** | `temporal_analysis.frequency_trends` (lines 126-142) | Slope and first→last change use reported frequencies only; adequate-depth non-detections are counted but ignored | TIME-3/TIME-4 | Trend direction can be biased toward detections | Treat adequate-depth non-detections as below-detection observations, or show them on the trend; add depth-aware uncertainty (Sprint 5) | Yes (Sprint 5 builds on it) |
| S3-03 | **P2** | `diversity_analysis.paper_distance_matrix` / `population_structure` | Positions without shared reads contribute 0, so samples with no common coverage are distance 0; `compared_positions` is not reported or thresholded | YC-7/YC-8 | Low-coverage samples can appear similar in the paper-method PCoA | Report shared callable positions per pair; set a minimum, or flag pairs; consider normalizing by compared positions as a labelled variant | Yes |
| S3-04 | **P2** | `integrated_analysis.run_integrated_analysis` lines 312-319; `figures.trajectory_figure` | Trajectory figure drawn without dates → matplotlib conversion error shown to the user | Scenario 12 | One missing figure; confusing message | Skip the figure with a plain note when no included sample has a date | No (but cheap) |
| S3-05 | **P2** | `processed_snv.load_processed_snvs` (lines 153, 169) | Ambiguous dates such as `03/04/2021` parse month-first silently; mixed formats are rejected | Date test: `03/04/2021` → 2021-03-04 | Wrong dates feed temporal analysis | Accept ISO dates, or require the user to confirm day/month order at review | Yes |
| S3-06 | **P2** | `processed_snv.load_processed_snvs` (lines 296-302) | Reported frequency used when counts disagree; tolerance 0.05 is absolute; no option to compute VAF from counts | VAF-1/VAF-2 | Inconsistent tables pass with a warning only | Sprint 5: prefer counts when present; report disagreement per record | Yes (for Sprint 5) |
| S3-07 | **P2** | `quality_control.reference_check`; pooled-only runs | Without a reference FASTA, Ref bases and the reference label are trusted ("Not checked") | Code read; scenario 18 vs 19 | Wrong coordinate systems could pass in pooled-only mode | Recommend or require the reference FASTA for pooled analyses that compare positions | Yes |
| S3-08 | **P2** | `app.py` in-process state (`STAGED_LIMIT=6`, `RESULT_DOWNLOAD_LIMIT=4`, `ANALYSIS_SLOT`) | A 7th concurrent upload evicts the oldest staged set; a 5th result evicts the oldest downloads; state requires one worker | Code read ([app.py:194-209](../app.py#L194)) | Class demos with many users can lose uploads or links | Document; raise limits within memory, or move staging to disk-only metadata | No |
| S3-09 | **P2** | Hosting | 131 s at 10% of one core; the Render proxy's request timeout and cold start were not measured | Section C.1 | Users may give up or hit a proxy timeout (unverified) | Verify on Render; consider a background job with a status page | No |
| S3-10 | **P2** | `references.load_reference_fasta` (exactly one sequence) | One reference sequence per run; segmented genomes (for example influenza) can use one segment only | Code read ([references.py:41-44](../genomics/references.py#L41)) | Limits non-SARS-CoV-2 pooled data | Multi-record reference registry keyed by segment | Only if Sprint 5 targets segmented pathogens |
| S3-11 | **P3** | Review step (`input_inspection`) | Invalid nucleotides and blank sample IDs are caught only at run | Scenarios 13, 14 | One extra round trip | Reuse loader checks in `value_warnings` | No |
| S3-12 | **P3** | Review messages | When one column is missing, "Required information missing" lists **every** pooled field | Scenarios 6, 8 | Misleading | List only fields not mapped | No |
| S3-13 | **P3** | `integrated_analysis.position_summary` (line 100) | Median includes duplicate rows (DUP-2) | DUP-2 | Small numeric inconsistency | Drop duplicates by sample before the median | No |
| S3-14 | **P3** | Excel reading | Workbook opened 4 times per guided run; secondary XLSX tables read via `pd.read_excel` (first sheet, header row 1) in `read_any_table` ([pooled_observations.py:82-93](../genomics/pooled_observations.py#L82)) outside the guided path | Instrumented loads | Small time cost; two Excel readers | Use `excel_tables.read_sheet` everywhere; cache the sheet choice | No |
| S3-15 | **P3** | `dashboard_guide.py` PCoA entry (line 143) | "Measures … nucleotide frequencies" also shown for the Jaccard ordination | Code read | Slight overstatement (the figure title is correct) | Method-specific text | No |
| S3-16 | **P3** | `validation/check_fontenele_2021.py` line 48 | Row "Samples with >90% coverage (analyzed)" in the "Application calculates" column only counts sample IDs; coverage is not in the table | Code read | Implies a calculation that is not made | Relabel "Samples in table (paper: >90% coverage)" | No |
| S3-17 | **P3** | Paper-format detection ([processed_snv.py:95-100](../wastewater/processed_snv.py#L95)) | Recognition depends on exact title-row phrases | Code read | A re-saved copy without the title row falls back to column mapping | Detect by header set plus reference | No |
| S3-18 | **P3** | Data caveat | Published `Sy/Ns` labels UTRs as "Syn"; no own annotation exists | Workbook inspection | Sprint 5 must not reuse these labels for non-coding sites | Annotate from a GFF instead | Yes (for annotation) |
| S3-19 | **P3** | Tests | No automated test runs `mmc2.xlsx` through `/review` → `/run`; the optional test checks counts only. The memory limit has no regression test. | Test inventory | A future regression could pass CI | Add an optional web-path test plus a memory-budget test | Recommended |

---

## F. Test and performance findings

### F.1 Automated tests (run in this audit)

```bash
FONTENELE_MMC2=<path to mmc2.xlsx> AI_ENRICHMENT_ENABLED=false PATH=/c/Tools/mafft-win:$PATH .venv/Scripts/python.exe -m unittest discover -s tests -v
```

- **Result:** `Ran 409 tests in 242.258s — OK (skipped=1)`.
- **Skipped:** `test_optional_real_ncbi_data` (needs NCBI files in the project root).
- **Real data:** the real-workbook test `test_optional_fontenele_table` **ran and passed**.
- **Environment:** Python 3.11.9 (`.venv`), MAFFT v7 on PATH, Windows 11.

Sprint 3-related test files (test count in brackets):

- `test_wastewater` (9), `test_quality_control` (6), `test_analyses` (10: Yue–Clayton, PCoA, Jaccard, trends, regression);
- `test_excel_import` (22: header row, worksheets, size limits, web review), `test_genomics_schema` (7), `test_input_detection` (27), `test_mode_upload` (12), `test_guided_workflow` (12);
- `test_integrated_example` (13), `test_audit_regressions` (20), `test_offline_workflow` (3), `test_environmental_analysis` (14), `test_dashboard` (6), `test_optional_fontenele_table` (1).

**What the tests use:**

- Synthetic fixtures and hand-made tables, plus one optional real-data test (counts only).
- Web routes are exercised with the Flask test client in 16 test files.
- MAFFT is real when installed, otherwise `tests/fake_mafft.py`.
- GeoNames, NCBI, and Claude are mocked.

### F.2 Independent checks (this audit)

- `audit_numeric_checks.py`: 46 checks, 44 agree. The 2 that differ are VAF-1 (reported VAF is used by design) and DUP-2 (S3-13).
- `audit_scenarios.py`: 22 web-path scenarios (section F.4).
- `audit_combined_vs_pooled.py`: 10 of 11 pooled tables identical; 1 differs as intended.
- `audit_workbook.py`: workbook structure; REF agreement with NC_045512.2 (7,973 / 7,973).

### F.3 Performance (measured locally)

All figures are measured in this audit on Windows 11, `.venv` Python 3.11.9, `OPENBLAS_NUM_THREADS=1`, Flask test client. Peak is the process working set; caps use a Windows Job Object.

| Dataset | Cells (variants × samples) | Code | Peak | Time |
|---|---|---|---|---|
| `mmc2.xlsx` | 295,000 | pre-fix | 499–500 MB | 16.6–17.2 s |
| `mmc2.xlsx` ×4 in one process | — | pre-fix | 537 MB | 68 s |
| `mmc2.xlsx` | 295,000 | current | 286–288 MB | 12.5–13.8 s |
| `mmc2.xlsx` ×4 in one process | — | current | 302 MB | 55 s |
| `mmc2.xlsx`, 512 MB and 10% of one core | — | current | 279 MB | 131 s |
| Generated 10,000 × 60 | 600,000 | current | 394 MB | 14 s |
| Generated 15,000 × 100 | 1,500,000 | current (no limit) | 665 MB | 26 s |
| Generated 10,000 × 60, 512 MB cap | 600,000 | pre-fix | `MemoryError` | 5.6 s |

**What each stage scales with (from the code; inferred):**

| Stage | Scales with |
|---|---|
| Observation matrix, heatmap input, trends | variants × samples |
| Paper Yue–Clayton distances | samples² × reference length (time); samples × length × 4 (memory, float32: about 25 MB for 52 samples on SARS-CoV-2) |
| Jaccard | samples² × variants per sample |
| PCoA | samples³ (eigen-decomposition), samples² memory: negligible below a few hundred samples |
| Cross-source comparison | variant keys × (individual + pooled samples) |
| Excel loading | rows × columns (streaming) |
| Downloads | total CSV size; files > 8 MB are skipped and listed |

The dominant term for realistic data is variants × samples, which `MAX_OBSERVATION_CELLS` guards. The quadratic terms are not limiting at current sizes.

### F.4 Diagnostic scenario matrix (web path: `/review` → `/run`)

| # | Scenario | Expected | Actual | Crash? | Science correct? | Limitation explained? |
|---|---|---|---|---|---|---|
| 1 | Small valid CSV (9 samples) | Runs | Runs; 13 figures; exploratory ordination | No | Yes | Yes (coverage and PCoA warnings) |
| 2 | Same as XLSX | Runs | Identical outcome | No | Yes | Yes |
| 3 | Two data worksheets | Ask user | "Multiple worksheets were detected… Select the worksheet"; run blocked until chosen | No | n/a | Yes (message shown twice) |
| 4 | Empty worksheet | Reject | "Worksheet 'empty' is empty" | No | n/a | Yes |
| 5 | Lowercase / spaced headers | Map | Mapped; runs | No | Yes | — |
| 6 | No `Read_Depth` column | Reject with name | "no column was recognized for read depth (Read_Depth)" plus an over-long missing list (S3-12) | No | n/a | Partly |
| 7 | `Read_Depth` blank | Run, flag | Runs; depth QC "Not evaluable"; all retained | No | Yes (no invented depth) | Yes, in the method table |
| 8 | No `Alt_Frequency` column | Reject | Specific message | No | n/a | Yes |
| 9 | Text in `Position` | Reject | "1 value(s) are not numbers (for example 'about 80')" | No | n/a | Yes |
| 10 | Frequencies as percent | Reject | "values above 1, so it may be a percentage… not rescaled automatically" | No | n/a | Yes |
| 11 | Duplicate rows | Warn | Runs; warning; first row used | No | Yes (S3-13 aside) | Yes |
| 12 | No collection dates | Run without temporal | Runs; **trajectory figure error text** (S3-04) | No | Yes | Partly |
| 13 | ALT = `R` | Reject at review | Rejected at run, clear message (S3-11) | No | n/a | Yes |
| 14 | Blank Sample_ID | Reject at review | Rejected at run (S3-11) | No | n/a | Yes |
| 15 | Two pooled samples | No PCoA | "Population structure: Not available" | No | Yes | Yes |
| 16 | SNVs + A/C/G/T counts + reference | Paper PCoA | "Paper method: Yue-Clayton PCoA"; 8 of 9 samples kept (coverage QC); off-reference records incomparable | No | Yes | Yes |
| 17 | SNV reference ≠ uploaded FASTA | Incomparable | Warning at review; records incomparable; pooled analysis still runs | No | Yes | Yes |
| 18 | `mmc2.xlsx` | Runs, exploratory | 52 samples; 13 figures; exploratory Jaccard | No | Yes | Yes |
| 19 | `mmc2.xlsx` + NC_045512.2 FASTA | Runs, Ref verified | 52 samples; 13 figures; Ref check applied | No | Yes | Yes |
| 20 | `mmc3.xlsx` | Reject | "contains variant positions but not the processed SNV values… no column was recognized for pooled sample ID…" | No | n/a | Yes |
| 21 | Malformed CSV | Reject | "The columns do not match any supported table" | No | n/a | Yes |
| 22 | Sprint 2 + 3 combined (synthetic) | Runs both | 21 figures; paper PCoA; cross-source categories | No | Yes | Yes |

**Not performed:** any live Render test; a Docker build; the offline FASTQ workflow; Linux memory measurement.

---

## G. Sprint 5 readiness matrix

| Planned feature | Existing code | Missing | Required input fields | Scientific validation needed | Approach | Depends on |
|---|---|---|---|---|---|---|
| 1. Reference-based SNV comparison | `variant_key`, `references_compatible`, `harmonize_pooled_variants`, `compare_evidence` (verified) | Multi-segment references; indels | Reference FASTA; `Reference_ID`, `Position`, `Ref`, `Alt` | Already hand-checked; add segment tests | Reuse | S3-07, S3-10 |
| 2. Read-depth and VAF analysis | Reported VAF and depth validated; QC depth filter | VAF from counts; ref/alt/other counts in the schema | `Alt_Depth` (or ALT_DP), `Read_Depth`; ideally `Ref_Depth` | Agreement checks; percentage/fraction detection | Extend schema with count columns; compute VAF when counts exist | S3-06 |
| 3. Uncertainty for measured allele frequencies | None | Binomial / beta-binomial intervals | Alt and total counts per record (`mmc2`: derivable as ALT_FREQ × TOTAL_DP) | Coverage of intervals by simulation; overdispersion from replicates | New module `genomics/frequency_uncertainty.py` | Item 2 |
| 4. Mutation-frequency trends over time | Observation matrix, statuses, OLS slope | Depth-aware model; non-detections as data; uncertainty | Dates (ISO), depth at non-detected positions (position table) | Selection-bias test (TIME-4); simulation; minimum dates per series | New trend module reading the existing matrix | S3-02, S3-05; position-depth tables |
| 5. Depth-aware mutation grouping (MAGOS-inspired) | None (Sprint 3 has no clustering; PCoA is visualization only) | Entire method | Alt/total counts per variant per sample; multiple samples | Simulation with known mixtures; stability; see section G.1 | Separate optional module | Items 2–3 |
| 6. Gene and protein annotation | Passthrough of published/iVar annotation (`Published_*`, `Source_*`) | Own annotation from GFF | Reference GFF matching the FASTA | Compare with published annotations; overlaps (ORF7a/7b) | New module using reference + GFF | S3-18 |
| 7. Synonymous / non-synonymous | Published labels only (UTRs mislabelled "Syn") | Codon translation, frame handling, ribosomal slippage (ORF1ab) | Reference FASTA + GFF | Reproduce published labels in coding regions; flag UTR | Same module as item 6 | Item 6 |
| 8. Individual vs pooled comparison | Implemented and verified | Frequency-aware summaries | As now | Already checked | Reuse | — |
| 9. Public-health interpretation | Cautious wording throughout | Guidance tied to uncertainty | — | Expert review of wording | Text plus limits in the dashboard guide | Items 3–4 |

**Capability assessment (what exists today):**

| Capability | Status |
|---|---|
| Reference read counts | **Unsupported** (not in schema; not in `mmc2`) |
| Alternate read counts | **Partial** (`Alt_Depth` optional; `mmc2` derivable, not stored) |
| Total usable depth | **Implemented** (`Read_Depth`; source definition) |
| Per-position VAF | **Implemented** (reported values only) |
| Per-position uncertainty | **Unsupported** |
| Multiple time points | **Implemented** (one date per sample) |
| Comparable mutation IDs across samples | **Implemented** (`REF:POS:REF>ALT`) |
| Full coverage / callability matrices | **Dependent on input** (position table); not available for `mmc2` |
| Reference FASTA mapping | **Implemented** (one sequence) |
| Gene annotations | **Partial** (passthrough only) |
| Amino-acid consequences | **Partial** (passthrough only) |

### G.1 MAGOS: what does and does not transfer

MAGOS (Ahmadinejad et al., *Mol Biol Evol* 2022, 39(7):msac136) groups somatic mutations into tumor subclones.

**What the method does** (from the verbatim part of the article):

- It models each cluster's VAFs with a moment-matched beta distribution whose variance falls with mean depth: `var = v̄(1−v̄)/(ē+1)`.
- It merges variants by agglomerative clustering on a weighted negative log-likelihood.
- It splits the resulting tree adaptively with a variance test (p < 0.01).
- With several samples, a group must hold in every sample.

The article was read verbatim only for about its first 40%. The Discussion, Methods details, and supplement came from a summary and were not verified.

**Assumptions that do not transfer to pathogens or wastewater:**

1. **Cellular prevalence = 2 × VAF.** This depends on heterozygous mutations in a diploid genome. Viral genomes are haploid, and wastewater mixes many genomes from many people, so VAF is a mixture proportion of genomes, not a cell fraction.
2. **Copy-number correction and tumor purity.** These have no meaning for viral RNA. Amplicon dropout and primer bias are the relevant distortions instead.
3. **Shared VAF within a group means shared lineage.** In wastewater, mutations at similar frequency may come from different lineages whose abundances happen to be similar. Similar VAFs **do not** show that mutations sit on the same genome. Short reads rarely span two distant sites, so linkage cannot be confirmed.
4. **Concordant change across samples.** In tumors, samples come from one patient. Wastewater samples from different dates or sites have different source populations, so concordance is weaker evidence.
5. **Sampling variance.** The beta variance assumes reads are binomial draws. Amplicon PCR adds overdispersion, and RNA degradation and low input create jackpot effects. The variance–depth relationship must be re-estimated from replicate or pathogen data.
6. **Validation.** The article validated with simulated tumors at 30–500× and with AML data. A pathogen version needs its own simulations: mixtures of known lineage genomes at known proportions, amplicon-level depth noise, and dropouts.

**Recommended validation before any pathogen grouping is shown:**

- Simulate mixtures of 2–5 known genomes (for example published SARS-CoV-2 consensus sequences) at set proportions and realistic per-amplicon depth.
- Measure how often grouped mutations recover the true lineage-defining sets, across depth (30–1,000×) and proportion gaps.
- Estimate overdispersion from technical replicates, if available, or from the published depth–frequency spread.
- Test stability under bootstrap resampling of reads.
- Report groups as "co-varying mutation sets (exploratory)", never as lineages or haplotypes, with per-group uncertainty and the minimum depth used.

### G.2 Architecture for expansion

**Recommendation: add separate optional modules that read the existing shared tables.** Do not extend the Sprint 3 functions.

- **Keep as is:** the existing variant and sample schema (`unified_schema`), the QC, `observation_matrix`, and the reference registry. They are a sound base.
- **Add:**
  - optional `Ref_Depth`, `Alt_Depth`, and `Other_Depth` columns in the schema and loader. Today `Alt_Depth` is accepted;
  - new modules, for example `genomics/frequency_uncertainty.py`, `genomics/depth_aware_trends.py`, `genomics/mutation_grouping.py`, and `genomics/annotation.py`, each with its own tests and a clear "requires" check;
  - entries in `analysis_status` that say *why* each is unavailable, as Sprint 3 already does.
- **Run Sprint 5 analyses only for compatible datasets:**
  - counts present for uncertainty and grouping;
  - a position-depth table for absence-aware trends;
  - a GFF for annotation.

**Architecture risks:**

- `run_integrated_analysis` is already 300 lines in one function ([integrated_analysis.py:160-472](../genomics/integrated_analysis.py#L160)). Adding Sprint 5 there would make it hard to review. A small registry of analysis steps would help.
- The observation matrix is variants × samples. Adding per-cell intervals increases memory, so keep `MAX_OBSERVATION_CELLS` and measure.
- Single-reference assumptions run through the code: one FASTA, and the figures pick one reference.

---

## H. Recommended repair sequence (not implemented)

1. **Critical before deployment**
   - S3-01: after your review, push `b641533`, rebuild the Docker image, and verify on Render. Upload `mmc2.xlsx` twice; Metrics should show a peak well below 512 MB, and Events should show no "Ran out of memory".
   - Record which branch Render deploys.
2. **Scientific corrections before Sprint 5**
   - S3-05: dates — prefer ISO and confirm day/month order.
   - S3-02: trends that account for adequate-depth non-detections, or at least show them.
   - S3-03: report and threshold shared callable positions in the Yue–Clayton distance.
   - S3-06: VAF from counts when counts exist; per-record disagreement column.
   - S3-07: recommend or require the reference FASTA for pooled position work.
   - S3-18: own annotation rather than published Sy/Ns.
3. **Performance**
   - S3-14: open the workbook once per step.
   - Optional background job for long runs (S3-09).
   - Add a memory regression test (S3-19).
4. **User interface**
   - S3-04: skip the trajectory figure without dates.
   - S3-11: validate bases and IDs at review.
   - S3-12: list only missing fields.
   - S3-15, S3-16: wording.
5. **Optional**
   - S3-10: multi-segment references.
   - S3-08: concurrency limits.
   - S3-17: sturdier paper-format detection.
   - S3-13: duplicate-aware median.

---

## I. Final verdict

**READY AFTER SPECIFIC FIXES.**

Why not "Ready":

- The crash fix has not been deployed or tested on Render (S3-01).
- Several parts of the analysis need guards before depth-aware methods are built on them: dates (S3-05), trend selection (S3-02), shared-coverage-blind Yue–Clayton sums (S3-03), and count-free VAF (S3-06).
- The published workbook does not contain the reference counts, per-position depth, or A/C/G/T counts that depth-aware absence, uncertainty for reference alleles, and the paper's PCoA need. Those tables must be regenerated offline from the raw reads (SRA PRJNA662596) or come from another dataset.

Why not "Not ready":

- No P0 is open.
- The formulas that exist are correct by independent check.
- The real workbook runs end to end in about 286 MB.
- Combining with Sprint 2 does not change Sprint 3 results.
- The code separates "not reported", "not covered", and "0%" consistently, which is the foundation Sprint 5 needs.

---

## Appendix A. Observed, inferred, and unverified

| Claim | Status |
|---|---|
| Pre-fix code peaks at about 500 MB on `mmc2.xlsx` and loses figures under 440 MB | Verified (local) |
| Pre-fix code raises `MemoryError` at 600,000 cells under 512 MB, with a misleading message | Verified (local) |
| Render killed the original run for memory | **Unverified** (no Render logs) |
| Current code finishes in 131 s at 0.1 CPU equivalent | Verified (local Windows Job Object) |
| Render's proxy allows a 131 s request | **Unverified** |
| Linux memory use is close to the Windows figures | **Unverified** (Inferred: same Python libraries) |
| `mmc2` coordinates are NC_045512.2 | Verified (7,973 / 7,973 REF bases match) |
| `ALT_FREQ × TOTAL_DP` gives exact alternate counts | Verified (all rows are whole numbers) |
| MAGOS details beyond the first 40% of the article | **Unverified** (summary only) |

## Appendix B. Reproducible commands and files

Scratchpad: `C:\Users\ALEXAN~1\AppData\Local\Temp\claude\a7271a137b0bdb68\922255f3-0871-46d4-aedd-17d549d92b7a\scratchpad`.

Scripts (in its `scripts` folder):

- `repro_audit.py`: one route, counts workbook loads, reports peak memory. Uses `AUDIT_PROJECT` to choose the code.
- `repro_repeat.py`: the same with `REPEAT=n` runs in one process.
- `render_like.py`: runs a script under a memory cap (MB) and CPU cap (% of one core).
- `audit_numeric_checks.py`: hand-calculated checks (output `audit_numeric.txt`).
- `audit_workbook.py`: workbook structure and REF agreement (output `audit_workbook.txt`).
- `audit_scenarios.py`: the 22-scenario web matrix (output `audit_scenarios.txt` and `.csv`).
- `audit_combined_vs_pooled.py`: Sprint 3 outputs alone versus combined.

Data:

- `prefix_6e122ff\` — `git archive 6e122ff`, the pre-fix code.
- `refdata\NC_045512.2.fasta` — from NCBI E-utilities.

Examples (from the project folder, with `.venv`):

```bash
.venv/Scripts/python.exe <scratchpad>/scripts/render_like.py 512 10 <scratchpad>/scripts/repro_audit.py guided <scratchpad>/data/mmc2.xlsx
```
```bash
.venv/Scripts/python.exe <scratchpad>/scripts/audit_numeric_checks.py
```

Set `OPENBLAS_NUM_THREADS=1` for capped runs: without it, OpenBLAS fails to allocate per-core buffers under a 512 MB cap on a 32-core machine. The Dockerfile already sets this.
