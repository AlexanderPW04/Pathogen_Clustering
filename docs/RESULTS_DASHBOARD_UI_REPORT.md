# Results dashboard UI — report

**Branch:** `feature/homepage-sprint4-ui-20261116` (local only). It still holds the earlier homepage + Sprint 4 UI work, which is uncommitted and was kept; a snapshot of it is in `stash@{0}`.

**What was not done:**
- Nothing was pushed, merged, or deployed.
- The recovery branch is unchanged.
- No paid API call was made.

**Day-by-day log:** [RESULTS_DASHBOARD_UI_PROGRESS.md](RESULTS_DASHBOARD_UI_PROGRESS.md).

## 1. Original dashboard problems

**Where the uploaded 30-page PDF came from**
- It is a *browser print* of the results page, not the application's report.
- The page had no print styles, so the print included the step pills, toolbar buttons, the "Jump to" bar, every closed `<details>` summary, the system-status table, and large empty areas.

**Problems on the page itself (pooled-only run of the published mmc2.xlsx)**
- **Same sections for every run:** nine numbered sections, always rendered and titled "Sprint 2 + Sprint 3", even with no individual sequences.
- **Empty panels:**
  - "Cross-source evidence" said only that it was unavailable;
  - "Individual sequence relationships" had no analysis;
  - a genome-coverage chart drew every sample as "unknown";
  - a geography chart drew one gray bar of 52 named sites.
- **Overlapping figures:**
  - two frequency-trend figures;
  - five temporal figures in total;
  - a genome-position chart whose 516 multiple-base rings and 3,698 novelty stars buried the points near zero.
- **Unexplained counts:** 7,973 / 7,961 / 7,970 / 5,669 / 5,677 appeared in different places with no explanation.
- **Downloads and navigation:**
  - a download list under every section, plus a separate list at the bottom;
  - a 13-link navigation bar.

## 2. New conditional-rendering rules

**One central record**
- `genomics/result_availability.py` builds one record per result.
- The website and the PDF report both use it.
- Every rule reads the *outputs* of the run (tables, figures, and the existing `result["analyses"]` status record), never the mode alone.

| Decision | Rule |
|---|---|
| Mode and title | pooled only → **Pooled Pathogen Sequencing — Sprint 3 + 5**; individual only → **Individual Pathogen Sequences — Sprint 2**; both on one verified reference → **Combined Analysis — Sprint 2 + 3 + 5**; both without a shared verified reference → "Individual and Pooled Sequencing, analyzed separately". "+ 5" appears only when the depth-aware analysis ran. |
| Cross-source category | only when the cross-source comparison actually ran (both sources on one verified reference) |
| Individual relationships | only when pairwise MAFFT distances exist (or individual results exist) |
| Genomic variation | only when reference-position variants or annotations exist |
| Coverage chart | main chart only when genome coverage was measured; with pooled data and no coverage, it is replaced by one sentence; with individual sequences only, the callable-position bars stay |
| Geography chart | only when at least one sample has source or GeoNames coordinates; otherwise a short precision summary |
| PCoA | only when it was computed; "Paper method" or "Exploratory — not the paper's method" is shown as before |
| Experimental grouping | always collapsed, with its "Profiles only / Experimental" badge and validation text |
| Sprint 4 | a category only after an enrichment report (approved changes); a pre-analysis publication check adds one line to the findings; otherwise one small link under Downloads |
| Unavailable or limited analyses | a short list in the overview, filtered to the run's mode. Individual-only rows are hidden in a pooled run and replaced by one "Not part of this analysis" line; if individual sequences were uploaded but could not be compared, it says so instead |

**Where the full record is kept:** the complete status table remains under "Files, connection details, and analyses that ran" and in `analysis_status.csv`.

## 3. Sections hidden or reorganized for each mode

| Mode | Categories shown (in order) |
|---|---|
| Pooled | Overview · Genomic Variation · Temporal and Geographic Patterns · Sequencing Quality and Uncertainty · Advanced Analysis · Downloads and Provenance |
| Individual | Overview · Individual Sequence Relationships · (Genomic Variation and Sequencing Quality only with a reference) · Temporal and Geographic Patterns · Downloads and Provenance |
| Combined | Overview · Cross-Source Variant Evidence · Individual Sequence Relationships · Genomic Variation · Temporal and Geographic Patterns · Sequencing Quality and Uncertainty · Advanced Analysis · Downloads and Provenance |

**Navigation**
- The navigation lists exactly the categories that render, with short labels (for example "Time and place" and "Quality and uncertainty").
- On phones it scrolls sideways.

**Removed**
- the "Sprint 2 + Sprint 3" tags;
- the section numbers;
- the per-section download lines;
- "Back to results overview" under the overview itself.

**Moved**
- "How this works" is now below the results on results pages;
- "Analyses predicted in Step 2" is now under Downloads and Provenance.

**Overview**
- the title and status ("Completed" or "Completed with limitations");
- the headline counts, all calculated;
- the main findings, without repeating the headline counts;
- an "Important limitations" box;
- the open "Unavailable or limited analyses" list;
- links to each category and a **Download results** button;
- collapsed: "What each count means", "Files, connection details, and analyses that ran", and "Input completeness".

## 4. Charts retained as primary

**Primary charts by category**

| Category | Primary charts |
|---|---|
| Cross-source | evidence categories |
| Individual relationships | genetic-distance distribution; relationship network |
| Genomic variation | variants along the reference genome; mutation annotation overview, inside the annotation panel |
| Temporal / geographic | mutation frequencies over time with uncertainty (or the plain trajectories when no uncertainty series exists); sampling dates and sites; geographic precision (only with coordinates) |
| Quality / uncertainty | read-depth distribution; sequencing depth and frequency uncertainty (the one main VAF / depth / uncertainty view); coverage bars (only when measured) |
| Advanced | PCoA (only when calculated) |

**"About this chart"** (closed by default) under every chart:
- **What this shows**
- **How to interpret it** (axes, colors, symbols)
- **Important limitation**

## 5. Charts moved to "Additional figures" (one click away)

| Category | Additional figures |
|---|---|
| Cross-source | evidence matrix, positions by category, timeline |
| Temporal | variant frequencies over time (duplicates the uncertainty series), changes between dates, missing / unknown observations, sample-by-variant heatmap |
| Quality | QC checks, SNVs per sample, depth vs SNVs, callable positions, "same VAF at different depths" |
| Advanced | binomial vs beta-binomial width (hypothetical rho); replicate concordance |
| Overview | metadata completeness, under "Input completeness" |

- **Model comparison:** the comparison table, examples, and hypothetical-rho scenarios are behind "Model comparison and hypothetical overdispersion scenarios". It is labelled "calibrated and hypothetical" when real replicates exist.
- **Annotation table:** the detailed mutation table is collapsed; its summary says "first 25 of N mutations; the CSV has every row".

## 6. Charts suppressed because required data were missing

- **Which charts:**
  - "Sequencing quality by source", when no coverage was supplied and there are no individual sequences;
  - "Geographic precision", when no sample has coordinates. The reason says either "named places, no coordinates" or "no locations".
- **Where the reason is shown:**
  - in its category (for example "Genome-wide coverage was not supplied, so coverage percentages could not be calculated.");
  - in the unavailable list ("Chart not shown");
  - next to the PNG in Downloads, which still offers it.
- **Calculation is unchanged:** the figures are still calculated exactly as before.

**Drawing-only figure changes**

| Figure | Change |
|---|---|
| Crowded genome-position figures (>300 variants) | smaller translucent points; multiple-base positions and published novelty flags drawn as two tick rows under zero; no negative axis ticks; a shorter x-axis label that names annotated and unannotated regions |
| Sampling timeline | the legend lists only the sources drawn |
| Profile heatmap | titled "profiles, not grouped" when there are no groups |
| Read-depth distribution | the caption mentions a right panel only when one exists |

## 7. Count and denominator explanations

**Investigated in code and confirmed on mmc2 + FASTA + GenBank. No calculation error was found.**

| Label now used | mmc2 | Definition |
|---|---|---|
| SNV observations in the upload | 7,973 | every uploaded record |
| Variant records after QC | 7,961 | passed the QC preset (12 records below 20x) |
| Unique mutation observations after duplicate handling | 7,970 | 3 identical repeats counted once; intervals are calculated for these, including QC-excluded records (their QC status is kept) |
| Distinct mutations after QC | 5,669 | different position:Ref>Alt changes in records passing QC (Genomic Variation, grouping eligibility) |
| Distinct mutations evaluated for annotation | 5,677 | every reported change, including 8 seen only in QC-excluded records |

**Where these labels appear**
- the "What each count means" table;
- the stat cards;
- the PDF report.

**Tests:** `test_displayed_counts_match_their_denominators` builds a synthetic table in which all five stages differ, then recomputes each value from the downloaded CSVs.

## 8. Download organization

**One area: "Downloads and Provenance"**
- **Download Results:** the XLSX workbook and the PDF report. Both are prepared on demand, as before.
- **Grouped files:** CSV and PNG files grouped by category. The first group is open; every stored file appears exactly once.
- **Provenance and input mapping:**
  - Step 2 predictions;
  - Step 3 readiness;
  - the provenance record;
  - "How each uploaded file was read (input mapping and import provenance)";
  - metadata matching.
- **Sprint 4:** the optional link.

**Omitted large files**
- They are stated explicitly, for example: "`pooled_observation_matrix.csv` is too large for the hosted download, so it is in neither the file list nor the workbook."
- The workbook is not claimed to contain them.

**What is unchanged**
- Figures keep their own "Download PNG" link.
- Tables keep their "Download …csv" link inside "Detailed tables".
- No download route, cache limit, or memory setting was changed.

## 9. PDF and print improvements

**Kept separate:** the browser print uses CSS in `templates/index.html`; the application report is built by `genomics/report_exports.py`.

**Browser print** (`@media print`)
- **Hidden:**
  - the step pills;
  - the toolbar;
  - the sticky navigation;
  - buttons, forms, and back links;
  - "How this works";
  - the system status;
  - the chart guides;
  - closed `<details>` and their summary lines;
  - the download lists. A one-line note replaces them.
- **Layout:**
  - figures, stat cards, conclusions, and table rows do not split across pages;
  - headings stay with the content that follows;
  - wide figures shrink to the page;
  - report-like typography.
- **Page counts (Chrome print-to-PDF):**

| Run | Pages |
|---|---|
| mmc2 pooled | about 10 (previously 30) |
| Pooled example | 7 |
| Individual | 4 |
| Combined | 13 |

**Application report** (`analysis_report.pdf`)
- **Page 1:**
  - the title and status;
  - the executive summary (the same headline counts as the website);
  - the main findings;
  - the important limitations;
  - the unavailable or limited analyses;
  - "What each count means".
- **Then:**
  - the dataset and methods;
  - short limitations by analysis;
  - a note that full tables are in the CSV and XLSX files.
- **Figures:** only the website's primary figures, so no empty coverage or geography chart.
- **Page count for mmc2:** 10 (3 text pages and 7 figures).
- **Older stored summaries** keep the previous layout.

## 10. Exact files changed (relative to the homepage snapshot `stash@{0}`)

**New files**
- `genomics/result_availability.py` — the central availability record;
- `results_layout.py` — categories and grouped downloads;
- `tests/test_results_dashboard.py` — 18 tests;
- `docs/RESULTS_DASHBOARD_UI_PROGRESS.md`, this report, and `docs/results_dashboard_screenshots/`.

**Changed files**
- **`integrated_view.py`:**
  - count labels;
  - the pooled-only page no longer shows the individual-variant stat;
  - attaches availability, figure views, categories, and download groups;
  - findings without repeated counts.
- **`templates/index.html`:**
  - the new results markup and macros;
  - "How this works" placement;
  - results and print CSS.
- **`genomics/report_exports.py`:**
  - the summary JSON carries the availability record;
  - the executive-summary PDF;
  - the report figure order.
- **`genomics/figures.py`:** the drawing-only changes listed in §6.
- **`genomics/integrated_analysis.py`:** the read-depth caption text only.
- **`validation/local_ui_preview.py`:** the results preview states.
- **Updated tests, layout and wording only; no scientific value changed:**
  - `test_dashboard.py`, `test_ui_navigation.py`, `test_mode_upload.py`;
  - `test_sprint5_depth_aware.py`, `test_sprint5_phase2_web.py`, `test_sprint5_phase3_web.py`, `test_sprint5_phase4_web.py`;
  - `test_guided_workflow.py`, `test_final_polish_navigation.py`, `test_phase4a_web.py`;
  - `test_excel_import.py`, `test_upload_enrichment.py`.
- **What the test updates did:** they replaced old layout strings ("1. Overview", "Pooled samples: 9", "Distinct SNVs evaluated", "Sprint 5 — Depth-Aware Variant Analysis") with the new headline sentences and category headings. Each assertion still checks the same underlying value.

## 11. Screenshots and local preview

**Run the preview**

```
python validation/local_ui_preview.py          # then open http://127.0.0.1:5058/preview
```

**Results states:**
- individual;
- individual + reference;
- pooled;
- combined;
- pooled without dates or locations;
- individual without coordinates;
- the Sprint 4 example.

**Real data:**
- `fontenele` and `fontenele-no-fasta` run the published mmc2.xlsx only when `FONTENELE_MMC2` points to a local copy, which is never in Git.
- Each run takes about a minute; pages are cached for 8 minutes.

**Screenshots in `docs/results_dashboard_screenshots/`:**
- mmc2 pooled at 1280×800, at 1024×768, as a full overview at 1280, and on mobile at 375 (overview);
- the combined overview;
- individual at 1024×768;
- pooled without dates;
- page 1 of the generated PDF report.

**Browser checks:**
- At 1280, 1024, and 375 px no page scrolls sideways.
- On mobile, the category bar scrolls on its own and the overview follows directly below it.
- Every state renders without a figure error, and the navigation matches the rendered categories.

## 12. Full test results

**Command:** `AI_ENRICHMENT_ENABLED=false python -m unittest discover -s tests`

**Result: Ran 738 tests — OK (skipped=16).**

- **Count:** 720 tests before this work, plus 18 new tests in `tests/test_results_dashboard.py`.
- **Skips:** the 16 skipped tests are the same optional ones as in the baseline. They need the real-data variables (`FONTENELE_MMC2`, `CRITS_CHRISTOPH_DIR`, `SCHUMANN_DIR`, `ANNOTATION_MPOX_GENBANK`) or MAFFT on `PATH`.
- **Real data:** the mmc2 configurations were exercised separately through the preview and the CSV comparison in §13.

**How the earlier failures were resolved**
- **33 old assertions** failed on wording or layout during development and were updated (§10). No scientific expected value was changed.
- **Two real wording bugs** were found and fixed:
  - the "no individual sequences were uploaded" line appeared when MAFFT was unavailable;
  - the geography reason was wrong when there were no locations.

**Configurations covered**

| Configuration | Covered by |
|---|---|
| Individual only, pooled only, combined | new tests and preview |
| With Sprint 4 | the existing enrichment end-to-end test and the Sprint 4 preview state |
| Without Sprint 4 | the small-link test |
| Missing metadata (dates and coordinates) | missing-data preview states and the existing missing-date tests |
| Missing reference | pooled without a FASTA ("supplied but unverified"); individual without a reference |
| Missing dates | new test |
| Missing coverage | new test |
| No paper-method PCoA | the exploratory-label test |
| No annotation | the pooled example (annotation "Skipped", listed in the overview) |
| Experimental grouping off | the collapsed-panel test; the default stays "profiles" |

## 13. Scientific-output equivalence

**Method:**
- every CSV download was compared by SHA-256 between a clean worktree of `145dc4c` (code equal to staging `8690561`) and this branch;
- configurations: the integrated example, the Sprint 2 example, mmc2 + GenBank, and mmc2 + FASTA + GenBank.

**Result: all 120 CSV files are byte-identical.**

**PNG figures:** only the four drawing changes in §6 alter PNG bytes. No calculated value changed.

## 14. Remaining usability issues

- **Long unavailable list:** in missing-data runs it can reach 10–12 items; it is open by default and could be collapsed past about 6 items.
- **Annotation statistics:** the annotation panel still shows nine stat cards; a smaller consequence summary could replace some of them.
- **Synthetic example notice:** on example runs it still sits above the navigation.
- **Mobile screenshots:** only the overview was captured at 375 px; deeper sections were checked by measurement (no horizontal overflow).
- **Not yet done:** there has been no human usability test and no staging deployment of this layout.
- **Combined report length:** the combined example's generated report is 11 pages, against 10 before, because the executive summary adds text before the 8 figures.

**RESULTS DASHBOARD UI READY FOR REVIEW**
