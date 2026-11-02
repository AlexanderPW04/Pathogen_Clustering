# Sprint 3 Change Log

This file records work added **after the Sprint 2 prototype**. Use it as the source
for a later Sprint 3 report or presentation. Add an entry whenever a Sprint 3
change is completed. Do not list an existing Sprint 2 file as a new Sprint 3 file
just because it is used again.

## Starting point: Sprint 2 baseline (October 19, 2026)

The following files and capabilities already existed before Sprint 3 work began:

| Existing file or folder | What Sprint 2 already provided |
|---|---|
| `validate_inputs.py` | FASTA plus CSV/XLSX metadata loading, NCBI column aliases, accession matching, and validation. |
| `compare_genetic_distances.py` | Pairwise comparison of aligned sequences within the same pathogen and gene/segment group. |
| `add_temporal_geographic.py` | Days apart, Haversine kilometers when coordinates exist, and normalized location matching. |
| `detect_potential_clusters.py` | Threshold filtering and connected-component detection for potential clusters. |
| `run_pipeline.py` | One command that runs the four backend stages in order. |
| `app.py` and `templates/index.html` | Local Flask upload form, threshold inputs, and cluster display. |
| `requirements.txt` | pandas, openpyxl, Biopython, and Flask dependencies. |
| `tests/` | Six test files covering input validation, all analysis stages, the pipeline, and Flask. The last reported Sprint 2 run passed 27 tests. |
| `outputs/` | Sprint 2 demonstration inputs and generated results, including a 10-record NCBI-style example. |
| `README.md` | Cross-platform setup and usage instructions. |
| `PYTHON_LINE_BY_LINE_GUIDE.md` | Beginner-friendly explanation of the six main Python files. |
| `output/pdf/BMI461_Sprint2_Code_Review_Guide.md` and `.pdf` | Sprint 2 conceptual documentation. |

The folder was not a Git repository when this written baseline was created.
It is now versioned, and the entries below remain the human-readable change
record. The existing folder name `sprint_2` does not make later work part of
Sprint 2.

## Sprint 3 changes so far

### November 1, 2026 — Guided upload, more formats, flexible metadata, and richer integrated dashboard

This responds to the professor's review of the Sprint 2 + Sprint 3
integration (more meaningful visuals, a friendlier site, more file formats,
flexible metadata, a clearer integration, and more of Fontenele et al. 2021).
It is still Sprint 2 + Sprint 3 work, not Sprint 4.

- **Added files:** `genomics/file_formats.py` (gzip, FASTA/table/VCF/FASTQ
  recognition, VCF-to-SNV conversion), `genomics/column_detection.py`
  (role-specific header aliases, ambiguity, unit checks),
  `genomics/input_inspection.py` (file roles with evidence, dataset review,
  analysis plan, standardized copies), `validation/check_fontenele_2021.py`,
  `ARTICLE_METHOD_VALIDATION.md`, and the tests `test_input_detection.py`,
  `test_guided_workflow.py`, and `test_optional_fontenele_table.py`.
- **Changed files:** `app.py` adds the `/review` and `/run` routes and a
  30-minute temporary staging folder; the named-field `POST /` still works.
  `templates/index.html` is a three-step page (upload, review, results) with a
  prominent Run Integrated Example card, collapsible advanced settings, and one
  layout for all eight sections. `integrated_view.py` gives every section a
  purpose, statistics ("Not available" instead of 0), figures, tables,
  limitations, and downloads, and puts the integrated sections first.
  `genomics/figures.py` adds 13 figures (metadata completeness, variant-overlap
  matrix, cross-source positions, sampling timeline, frequency changes,
  missing observations, geographic precision, QC outcomes, depth distribution,
  callable positions, genetic distances, relationship network) and improves
  the genome-position and PCoA figures. `genomics/integrated_analysis.py`
  records per-sample cross-source evidence, which analyses ran, and how the
  sources were connected. `genomics/temporal_analysis.py` computes trends with
  grouped sums (the real 7,973-record table went from 106 s to about 8 s; same
  values). `validate_inputs.py` and `run_pipeline.py` relax individual
  metadata requirements (see below). `wastewater/processed_snv.py` reads SNV
  VCF files. `run_integrated_analysis.py` accepts `.gz` inputs and draws the
  Sprint 2 distance figure. `geocode_locations.py` labels a blank location
  "No location given".
- **Removed files:** None.
- **Scientific rule changes (documented and tested):** individual metadata now
  requires only `Accession` and `Pathogen`; missing dates, locations, or hosts
  are notes, not errors, and disable only the affected screening. A year-only
  or year-month date is no longer turned into a specific day. Units are never
  rescaled automatically. VCF frequencies are taken only from allele-specific
  values; genotype-only entries are not converted to frequencies.
- **Before:** one long form with eight separate file fields; exact column
  names; eight sections with few figures; pooled-only or individual-only
  sections were mostly blank.
- **After:** files are uploaded together, recognized by content, and reviewed
  with stated evidence before running; 21 figures across the sections; the
  overview explains how the two sources were connected and which analyses ran.
- **Verification:** 204 tests run, 202 pass, 2 optional real-data tests skipped
  in the default run (both passed when their real files were supplied). Real
  MAFFT (Windows, `C:\Tools\mafft-win`) reproduced the synthetic example's
  categories and the real H5N1 baseline (45 comparisons, 22 candidate pairs, 0
  strict clusters). The real Supplementary Table 2 matched the paper's record,
  sample, state, per-sample, and D614G counts exactly. GeoNames and the Docker
  build were not checked (no username; Docker daemon not running).

### November 1, 2026 — Integrate Sprint 2 and Sprint 3 into one reference-based analysis

- **Added files:** `genomics/` package (`references.py`, `unified_schema.py`,
  `pooled_observations.py`, `sequence_variants.py`, `quality_control.py`,
  `temporal_analysis.py`, `diversity_analysis.py`, `variant_matching.py`,
  `figures.py`, `integrated_analysis.py`); `integrated_view.py`;
  `run_integrated_analysis.py`; the synthetic paired demo in
  `examples/integrated/`; the optional `offline_workflow/`; `tests/fixtures/`,
  `tests/fake_mafft.py`, and seven new test files;
  `INTEGRATION_IMPLEMENTATION_LOG.md`.
- **Changed files:** `app.py` and `templates/index.html` now form one dashboard
  for individual, reference, and pooled inputs (the Sprint 3 pages remain
  linked). `wastewater/processed_snv.py` keeps the paper's annotation columns
  and reads iVar TSV files. `run_pipeline.py` can hand its validated sample
  table to the integrated analysis. `align_sequences.py` factors out the
  unchanged MAFFT command. `validate_inputs.py` and
  `compare_genetic_distances.py` accept IUPAC ambiguity codes and treat them
  like N. `system_diagnostics.py` reuses its MAFFT version lookup.
  `Dockerfile` copies the new package; `requirements.txt` lists numpy.
  Two regression tests now use committed synthetic fixtures; the real-NCBI
  check became an optional, clearly skipped test.
- **Removed files:** None. The Sprint 3 branch's generated `tmp/pdfs/` files
  were left out of the merge.
- **Before:** Sprint 2 compared individual sequences pairwise on a
  sample-to-sample alignment; Sprint 3 summarized pooled SNVs and processed
  feature tables on separate pages. The three data formats were unrelated.
- **After:** Both evidence types share one reference-coordinate variant
  table. Individual sequences receive true 1-based reference positions;
  pooled records are checked with the Fontenele et al. (2021) settings where
  the data allow; trends, a Yue–Clayton PCoA, and a cross-source comparison
  (shared, source-specific, unresolved, incomparable) appear on one dashboard
  with the Sprint 2 cluster screening as one section.
- **Why:** The instructor's feedback was that the Sprint 3 addition was not
  integrated with Sprint 2 and should use more of the article's methods.
- **How it works:** See the README section "Integrated analysis (Sprint 2 and
  Sprint 3)" and `INTEGRATION_IMPLEMENTATION_LOG.md`.
- **Verification:** 164 tests run, 162 pass, 2 skipped (real MAFFT not
  installed; optional real-NCBI data). MAFFT and GeoNames are mocked in tests.
- **Limits or next steps:** Not yet checked with real MAFFT, a Docker build,
  or the real `mmc2.xlsx`. The exact PCoA needs per-position A/C/G/T counts.

### October 19, 2026 — Complete local setup README refresh

- **Changed files:** `README.md` now begins with Windows setup and repeat-run
  commands, with PowerShell and macOS/Linux alternatives. It documents all
  current pages, inputs, optional MAFFT/GeoNames setup, diagnostics, downloads,
  local troubleshooting, and optional Docker/Render hosting. The AI tools
  acknowledgment remains at the end.
- **Accuracy updates:** The test count is 110 on Python 3.11.9. Full-suite
  regression fixtures excluded from GitHub are listed separately from the
  portable measurement-page tests. Stale deployment and Python-version claims
  were removed.
- **Verification:** All five documented local GET pages returned HTTP 200
  through the Flask test client on Python 3.11. No analysis logic changed.

### October 19, 2026 — Wastewater Data Analysis for processed feature tables

- **Changed files:** `app.py` exposes `/wastewater-data` and its organized CSV
  download using the existing environmental-measurement module and template.
  `templates/index.html` links to the requested page. Earlier draft routes and
  files are retained. `environmental_analysis.py` now accepts missing optional
  dates, locations, frequencies, and coverage, and counts every listed feature.
  `templates/environmental.html` displays availability messages, sample/date/
  location counts, and descriptive charts. README and focused tests document
  and verify this first table-summary addition. `.dockerignore` keeps its
  bundled synthetic CSV available in the Docker image.
- **Missing values:** Supplied zeros remain zero. Missing numerical values
  remain blank in tables and CSV downloads, and all-missing dated measurements
  remain missing in summaries and plots. Frequency and coverage are never
  derived from counts.
- **Scope:** Descriptive pandas operations on already processed tables only.
  Sprint 2 validation, alignment, geocoding, genetic comparison, temporal/
  geographic rules, and cluster detection are unchanged.
- **Verification:** All 110 automated tests passed on Python 3.11.9, including
  uploads in both formats, missing measurements, supplied zeros, summary
  calculations, PNG charts, temporary downloads, and Sprint 2 regressions.
  Hosted deployment was not performed during this change.

### October 19, 2026 — Processed wastewater SNV surveillance mode

- **Added files:** `wastewater/__init__.py` and `wastewater/processed_snv.py`
  load, validate, and summarize processed calls separately from individual
  genomes. `templates/wastewater.html` presents the upload and results page.
  `tests/test_wastewater.py` covers published-column mapping, validation,
  metadata joins, counts, the webpage, and CSV download.
- **Changed files:** `app.py` adds `/wastewater` and a brief in-memory CSV
  download. `templates/index.html` links to the new mode. `Dockerfile`
  includes the new Python package. `README.md` documents the exact published
  columns, mapping, scientific limits, and local/hosted usage.
- **Before (Sprint 2):** The site screened comparable individual pathogen
  sequences for potential outbreak clusters. Pooled wastewater SNV calls were
  not an accepted input.
- **After (Sprint 3):** A separate page accepts Fontenele et al. Supplementary
  Table 2 or standardized processed SNV calls and summarizes observed variants.
  It retains and reports repeated keys instead of treating a missing call as
  confirmed mutation absence.
- **Source:** Fontenele et al. (2021), Water Research, Supplementary Table 2
  (`mmc2.xlsx`), read directly from the publisher's published file. Table 3
  (`mmc3.xlsx`) is a novel-SNV subset without the full frequency/depth fields.
- **Scope:** No FASTQ processing, BWA, iVar, Freyja, lineage assignment, or
  wastewater clustering. Existing Sprint 2 calculations and routes remain.
- **Verification:** The published workbook loaded as 7,973 retained records
  across 52 samples, with three repeated keys reported. Seven focused loader
  tests pass; full Flask/regression tests could not run in this workspace
  because the test interpreter lacks Flask/Biopython and dependency
  installation was not approved. They remain to be run locally or in Docker.

### October 19, 2026 — Environmental processed-measurement exploration

- **Added files:** `environmental_analysis.py` loads and validates processed
  CSV/XLSX measurements and builds sample, time, location, missing-data, and
  quality summaries. `templates/environmental.html` provides a separate upload,
  example, chart, and results page. `examples/synthetic_environmental_measurements.csv`
  contains clearly labeled invented data. `tests/test_environmental_analysis.py`
  covers file formats, validation, summaries, charts, and the Flask download.
- **Changed files:** `app.py` adds `/environmental` and its short-lived
  organized-CSV download. `templates/index.html` links to the separate mode.
  `requirements.txt` adds Matplotlib for descriptive plots, and `Dockerfile`
  includes the analysis module. `README.md` describes the schema, safeguards,
  synthetic example, and local/hosted usage.
- **Before:** The Sprint 2 pathogen pipeline handled FASTA plus metadata, while
  environmental processed measurements had no dedicated analysis page.
- **After:** Visitors can explore sample-level processed measurements by
  feature, collection date, and location, including missing-data and quality
  summaries. Results are descriptive only and do not establish a diagnosis,
  transmission, or outbreak.
- **Scope:** No raw reads, sequence alignment, Sprint 2 calculation changes,
  statistical significance testing, or diagnostic/outbreak classification.
- **Verification:** All 105 automated tests passed, including CSV and XLSX
  loading, validation, summaries, chart generation, Flask results, CSV download,
  and the existing Sprint 2 and other Sprint 3 regression tests. The tests use
  the installed Python 3.11 environment; hosted Render deployment was not
  performed from this workspace.

### October 19, 2026 — On-demand MAFFT and GeoNames status check

- **Added files:** `system_diagnostics.py` runs a small real MAFFT alignment of
  invented sequences. `tests/test_diagnostics.py` tests success and failure
  responses without external calls.
- **Changed files:** `geocode_locations.py` adds a separate one-request
  GeoNames diagnostic using its existing secure URL and SSL checks; the
  analysis geocoder still follows its original rules. `app.py` adds a
  five-minute in-process cache and a click-only status route.
  `templates/index.html` shows MAFFT and GeoNames statuses near the top.
  `Dockerfile` includes the new module, while `README.md` and Flask tests
  explain and verify local and hosted use.
- **Scope:** No changes to genomic distances, alignment grouping, geographic
  screening, thresholds, clusters, or wastewater analysis. `/health` remains
  separate from the dependency checks.
- **Verification:** All 87 automated tests passed with MAFFT and GeoNames
  mocked for the new diagnostics. The deployed Render service has not been
  updated or checked from this workspace.

### October 19, 2026 — Prepare the Sprint 2 website for Docker hosting

- **Added files:** `Dockerfile` installs Python 3.11, Linux MAFFT, Python
  requirements, and runs the existing Flask app with Gunicorn. `.dockerignore`
  and `.gitignore` keep local inputs, generated outputs, caches, and credential
  files out of the image or new Git commits.
- **Changed files:** `requirements.txt` adds Gunicorn. `app.py` limits uploads,
  adds `/health`, and keeps internal errors off public response pages.
  `templates/index.html` explains zero strict clusters when geography is
  unresolved. `README.md` adds Windows-to-GitHub-to-Render Free instructions;
  alignment and Flask tests cover Linux executable discovery, upload limits,
  health checks, and error presentation.
- **Verification:** All 74 automated tests passed. Docker CLI is present, but
  its daemon is not running/reachable here (and the local Docker config denied
  access), so no image build or Render deployment was verified.
- **Scope:** Deployment preparation for the existing six-stage Sprint 2
  website only. Scientific calculations, MAFFT grouping, GeoNames behavior,
  and strict-versus-candidate result meanings are unchanged. No wastewater
  analysis was added.

### October 19, 2026 — Clearer website counts and candidate-pair inspection

- **Added files:** None.
- **Changed files:** `app.py` counts distinct sample IDs in the pipeline's
  region-assignment output and in unresolved-geography candidate pairs.
  `templates/index.html` highlights sample and pair counts, keeps the full
  candidate list collapsed until selected, rounds displayed p-distances to
  six decimals, and labels example-run results as synthetic. `README.md` and
  `tests/test_app.py` explain and verify the presentation.
- **Before:** Candidate pairs were all shown at once, and a count of 22 pairs
  could be mistaken for 22 additional samples.
- **After:** The NCBI-style test shows 10 unique samples, 45 pairwise
  comparisons, zero strict related pairs/clusters, and 22 unresolved-geography
  candidate pairs involving 9 of the original samples.
- **Verification:** All 70 automated tests passed with MAFFT mocked where
  needed. No pipeline calculations or CSV precision changed.
- **Scope:** Website presentation only; no wastewater analysis.

### October 19, 2026 — Clearer Flask choices and geography-aware results

- **Added files:** `examples/synthetic_sequences.fasta` and
  `examples/synthetic_metadata.csv` provide an explicitly invented demo,
  including one name-only location and unequal sequence lengths for MAFFT.
- **Changed files:** `templates/index.html` now offers automatic detection or
  explicitly confirmed complete genomes, a Run Example Analysis button, a
  compact run summary, and a separate unresolved-geography candidate section.
  `app.py` handles those controls and reads the existing enriched CSV without
  changing strict cluster outputs. `validate_inputs.py` checks the new
  complete-genome assertion against metadata, FASTA, and pathogen identity.
  `run_pipeline.py` passes that option through and optionally supplies run
  counts to Flask. `README.md` and relevant tests were updated.
- **Removed files:** None.
- **Before:** The website had a free-text default region and a general
  same-region checkbox. It showed only strict clusters, so state-only data
  produced an empty result with little context.
- **After:** Unknown regions are not guessed in automatic mode. A confirmed
  complete-genome choice cannot override conflicting evidence. Pairs passing
  genomic/time thresholds but lacking usable geography are shown separately
  from all-criteria potential clusters. The synthetic demo uses the same
  six-stage backend.
- **Verification:** All 70 automated tests passed. The 10-record state-only
  NCBI-style upload produced zero strict clusters and nonzero unresolved-
  geography candidates with MAFFT mocked. The synthetic demo produced two
  strict clusters in the mocked-alignment test. A real MAFFT run was not
  possible here because `mafft` is not installed on `PATH`; the command-line
  run stopped with the expected installation message.
- **Scope:** No wastewater analysis, new scientific distance formula, or
  change to the existing strict related-pair and cluster CSV meanings.

### October 19, 2026 — Secure GeoNames endpoint and location precision

- **Added files:** None.
- **Changed files:** `geocode_locations.py` now calls the secure GeoNames URL
  with SSL verification enabled and labels state/administrative locations as
  too coarse. `add_temporal_geographic.py` keeps coordinates and sources in
  its output but leaves kilometer distance and same-location status blank for
  pairs involving coarse locations. `templates/index.html` and `README.md`
  explain the limit; geocoding and temporal/geographic tests were extended.
- **Removed files:** None.
- **Before:** The geocoder used an older URL, and a state centroid or identical
  state name could satisfy a local geographic cutoff.
- **After:** Named places can still be geocoded and cached, while state-level
  descriptions such as `USA: ND` cannot create a local-distance connection.
  The enriched CSV reports each sample's `Geographic_Precision`.
- **Verification:** All 57 automated tests passed. They check the secure URL,
  SSL verification, caching, state-name handling, and exclusion from the 50 km
  criterion. A live request was not run here because this shell has no GeoNames
  username; the README gives a command to run with a configured account.
- **Scope:** This corrects the inherited Sprint 2 backbone; it does not add
  wastewater or other Sprint 3 scientific analysis.

### October 19, 2026 — Optional genomic-region input and provenance

- **Added files:** `genomic_regions.py` reads explicit segment, gene, and
  whole-genome wording from metadata and FASTA descriptions without using
  sequence length or regular expressions.
- **Changed files:** `validate_inputs.py` now accepts flexible metadata headers,
  makes `Gene_Segment` optional, checks metadata/FASTA disagreements, and records
  each label's source. `run_pipeline.py` accepts a default or a same-region
  confirmation and writes `*_region_assignments.csv` after success.
  `compare_genetic_distances.py` accepts the same options when run alone.
  `app.py` and `templates/index.html` add equivalent controls and show sources
  beside clustered samples. `README.md` and input, pipeline, and Flask tests
  were updated.
- **Removed files:** None.
- **Before:** An exactly spelled `Segment` mapped to required `Gene_Segment`,
  but other spacing or capitalization did not. FASTA descriptions were ignored.
- **After:** The project resolves each region from metadata, then the FASTA
  header, then a user default, then explicit one-pathogen confirmation. It
  rejects contradictory evidence and does not merge different segments for
  MAFFT alignment.
- **Verification:** The real `sequences.csv` and `sequences.fasta` matched all
  10 records and resolved `Segment 4` from metadata. 53 automated tests cover
  the requested cases, including separate MAFFT calls for mixed segments;
  MAFFT itself remains unavailable locally, so that external call is mocked.
- **Scope:** This is another update to the inherited Sprint 2 backbone. It does
  not add Sprint 3 wastewater or SNV analysis.

### October 19, 2026 — MAFFT and GeoNames added to the inherited backbone

- **Added files:** `align_sequences.py` aligns each pathogen/gene group with the
  external MAFFT program. `geocode_locations.py` optionally fills missing
  coordinates from GeoNames and records their source. New alignment and
  geocoding tests were added under `tests/`.
- **Changed files:** `run_pipeline.py` now validates, geocodes, aligns, compares,
  enriches, and clusters in that order. `compare_genetic_distances.py` uses
  alignment when run directly. `add_temporal_geographic.py` includes coordinates
  and their provenance in the enriched CSV. `app.py` and `templates/index.html` show
  GeoNames warnings and short user guidance. `requirements.txt` adds `requests`.
  Existing integration tests, `README.md`, and `PYTHON_LINE_BY_LINE_GUIDE.md`
  were updated.
- **Removed files:** None.
- **Before (Sprint 2 baseline):** Comparable sequences had to arrive already
  aligned; missing coordinates used only the named-location fallback.
- **After (backbone update):** The pipeline invokes MAFFT separately for each
  pathogen/gene group and can estimate missing coordinates from GeoNames when
  `GEONAMES_USERNAME` is configured. Existing p-distance, temporal,
  geographic-distance, threshold, and connected-component formulas remain.
- **Why:** Accept unequal-length, unaligned sequences within a comparable
  group, and make location metadata more useful when coordinates are absent.
- **How it works:** MAFFT is found on `PATH` and called through `subprocess` on
  a temporary FASTA for each group. GeoNames searches a location once per run,
  returns approximate place coordinates when available, and labels them
  `GeoNames`; complete original coordinates are labeled `Provided`.
- **Verification:** 40 automated tests pass with Python 3.12 locally, including
  both existing example datasets through the full pipeline with a mocked MAFFT
  pass-through. The actual MAFFT executable is not installed on this laptop, so
  a real run remains to be checked after installation. GeoNames network calls
  are mocked in tests.
- **Limits or next steps:** The user requested this as a Sprint 2 backbone
  improvement, not as the Sprint 3 wastewater/SNV pipeline. No wastewater,
  FASTQ, BWA, iVAR, PANGOLIN, or SNV analysis was added. This entry records the
  post-baseline change for the future Sprint 3 report only.

### October 19, 2026 — Start the Sprint 3 record

- **Added:** `SPRINT_3_CHANGELOG.md` (this file).
- **Changed existing files:** None.
- **Removed files:** None.
- **Behavior added or changed:** None yet. The application still has the Sprint 2
  capabilities listed above.
- **Reason:** Establish a clear boundary so future additions can be explained as
  Sprint 3 work without accidentally presenting Sprint 2 features as new.
- **Verification:** Checked the current project files and existing documentation.
  No code was changed, so tests were not rerun for this entry.

## Add each future change here

Copy this small template for each completed Sprint 3 task, placing the newest
entry above older entries:

```markdown
### YYYY-MM-DD — Short name of change

- **Added files:** Names and purpose, or `None`.
- **Changed files:** Names and exactly what changed, or `None`.
- **Removed files:** Names and reason, or `None`.
- **Before (Sprint 2):** What the project could do previously.
- **After (Sprint 3):** What users can do now.
- **Why:** The problem or requirement this change addresses.
- **How it works:** A few plain-English sentences about the new code or data flow.
- **Verification:** Tests, sample run, and observed result; note if not tested.
- **Limits or next steps:** Any important unfinished work.
```

When a future change modifies an existing file, list that file under **Changed
files** rather than **Added files**. Record documentation changes too, so the
later Sprint 3 write-up can distinguish implementation from explanation.
