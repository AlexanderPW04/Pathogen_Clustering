# Results dashboard UI — progress log

**Branch:** `feature/homepage-sprint4-ui-20261116` (local only).
- The homepage + Sprint 4 UI work from earlier today is still uncommitted on this branch and was kept.
- Before any change, it was saved as `stash@{0}` ("snapshot: homepage + Sprint 4 UI uncommitted work"). The working tree was not touched by that snapshot.
- A copy of its untracked files was kept outside the repository.

**What was not done:**
- Nothing was pushed, merged, or deployed.
- The recovery branch `backup/sprint5-staging-validated-20261115` (`24f81fe`) is unchanged.
- No paid API call was made; `AI_ENRICHMENT_ENABLED=false` was used throughout.

## 2026-11-16

### 1. Inventory (no code changed)
- **The uploaded 30-page PDF** ("Integrated Genomic Surveillance.pdf") is a browser print of the results web page, not the application's own report (`analysis_report.pdf`). It contains:
  - the header step pills and toolbar buttons;
  - the "Jump to" bar;
  - every closed `<details>` summary line;
  - the system-status table.
- **The two PDFs come from different code:**
  - the application's report is built separately by `genomics/report_exports.py` (matplotlib `PdfPages`);
  - the browser print came from `templates/index.html`, which had no print styles.
- **How the page decided what to show:**
  - `integrated_view.build_dashboard` always built the same nine numbered sections, titled "Sprint 2 + Sprint 3".
  - The template rendered all of them; a section with nothing to show printed a gray "needs …" box.
  - Every section ended in its own download list.
  - Analysis status already existed as `result["analyses"]`, but the page did not use it to hide anything.
- **Real-data check** (pooled-only run of the published mmc2.xlsx with NC_045512.2 FASTA + GenBank, local copy only):
  - the empty genome-coverage chart (every sample "unknown");
  - a geography chart with one gray bar (52 named sites, no coordinates);
  - an empty "Individual sequence relationships" section;
  - a "Cross-source evidence" section that only says it is unavailable;
  - five temporal figures, two of which show the same frequency trends;
  - a genome-position chart whose 516 multiple-base rings and 3,698 novelty stars cover the points near zero.

### 2. Baseline
- **Clean worktree:** a separate worktree of `145dc4c` (code identical to staging `8690561`).
- **Configurations run there:**
  - the integrated example;
  - the Sprint 2 example;
  - mmc2 + GenBank;
  - mmc2 + FASTA + GenBank.
- **Baseline saved:** 120 CSV downloads (SHA-256) plus the page HTML and figure PNGs.

### 3. Counts and denominators (investigated in code and on the real data)

| Number | What it counts | Where it comes from |
|---|---|---|
| 7,973 | every uploaded SNV record | `len(pooled["variants"])` |
| 7,961 | records passing QC | `Analysis_Included` (12 records below the 20x depth limit) |
| 7,970 | records after 3 identical repeats are counted once | `depth_aware["observations_counted"]` (`Counted_Once`) |
| 5,669 | distinct mutations in records passing QC | `variant_positions.csv`, grouping eligibility |
| 5,677 | distinct mutations in all reported records | annotation (`distinct_variants` uses every record) |

- **Uncertainty intervals:** the 7,970 intervals include the 12 QC-excluded records; their QC status is kept.
- **Annotation:** of the 5,677 distinct mutations, 8 occur only in QC-excluded records. Annotation intentionally covers every reported change.
- **No calculation error was found.** The difference is a deliberate scope choice, now explained on the page and in the report.

### 4. Central availability record
- **New `genomics/result_availability.py`:** one JSON-serializable record per result, used by both the page and the PDF report. It holds:
  - the mode and title;
  - which analyses actually produced output;
  - the headline counts;
  - the important limitations;
  - the short "unavailable or limited" list;
  - the figures not shown and why;
  - the primary/secondary role of each figure;
  - the order of the report figures;
  - the count definitions.
- **Decided from outputs, not from the mode alone.** Examples:
  - coverage counts only when coverage was measured;
  - coordinates count only when a sample has precise ones;
  - the paper's PCoA counts only when the method is not the exploratory one.

### 5. Page layout
- **New `results_layout.py`:** arranges the existing section content into categories that fit the run:
  - Overview;
  - Cross-Source Variant Evidence (only with a verified shared reference);
  - Individual Sequence Relationships;
  - Genomic Variation;
  - Temporal and Geographic Patterns;
  - Sequencing Quality and Uncertainty;
  - Advanced Analysis;
  - Downloads and Provenance.
- **Empty categories are not rendered.**
- **Template (`templates/index.html`):**
  - a compact "Results:" navigation;
  - a new overview with the title and status, headline counts, findings, a limitations box, the open "Unavailable or limited analyses" list, links, and a collapsed "What each count means";
  - one figure-list pattern: main figures, then "Additional figures";
  - the uncertainty model table and hypothetical rho scenarios behind a labelled disclosure;
  - experimental grouping collapsed and visually de-emphasized;
  - the mutation table collapsed;
  - one grouped Downloads area;
  - provenance moved under Downloads;
  - Sprint 4 shown as a small link unless it was run;
  - "How this works" moved below the results on results pages;
  - print styles.
- **Chart explanations:** each chart's help is now "About this chart" with three parts: What this shows, How to interpret it, Important limitation.

### 6. Figures (drawing only; no calculation)
- **Genome positions (crowded figures only):**
  - smaller, slightly transparent points;
  - multiple-base positions and published novelty flags drawn as two tick rows below zero instead of rings and stars on the points;
  - no negative y-axis ticks;
  - a shorter x-axis label that names annotated and unannotated regions.
- **Sampling timeline:** the legend lists only the sources that are drawn.
- **Mutation-profile heatmap:** titled "profiles, not grouped" when no groups exist.
- **Read-depth distribution caption:** mentions a right panel only when there is one.

### 7. PDF report and print
- **`genomics/report_exports.py`:**
  - executive summary: title and status, headline counts, findings, important limitations, unavailable analyses, and count definitions;
  - then the dataset, methods, and short limitations;
  - then only the figures the page shows first;
  - older stored summaries without the availability record still use the earlier layout.
- **Mmc2 + FASTA + GenBank:**
  - the report has 10 pages: 3 text pages and 7 figure pages;
  - a browser print of the page has about 10 pages, against 30 before.

### 8. Preview and tests
- **`validation/local_ui_preview.py`:** gained results states:
  - individual;
  - individual + reference;
  - pooled;
  - combined;
  - pooled without dates or locations;
  - individual without coordinates;
  - Sprint 4 example;
  - mmc2 with and without a FASTA (only when `FONTENELE_MMC2` is set).
- **New tests:** `tests/test_results_dashboard.py`, 18 tests (see the report for results).

### 9. Problems found while testing, and fixed
- **Run with MAFFT unavailable:** individual sequences uploaded but not compared made the overview say "No individual sequences were uploaded". It now says the sequences were uploaded but could not be compared, both in the list and in the limitations.
- **Pooled table without locations:** the geography reason claimed the samples had place names. It now says no sample has a location.
- **Old test assertions:** 33 failed because they checked the earlier layout wording. They were updated to the new headings and headline sentences; the values they check did not change.
- **Wording kept from before:** "Sprint 4: …", "Sprint 4 metadata enrichment was enabled: …", and "Revisit metadata enrichment" were restored after a test showed I had changed them without need.
- **Preview banner:** it appeared twice on results pages; fixed.
