# Project handoff — Codex → Claude Code

**Purpose:** Transfer working context for a biomedical informatics class coding project. This handoff was prepared from discussion, Codex-reported summaries, and screenshots. It is **not** a code audit or an assertion that every file exists in the present repository. Verify the current files before editing.

## Sprint 4 on real NCBI records: preparation failure and evidence eligibility (November 2, 2026)

- **Bug and root cause:** on the real ten-record NCBI example (QB228092.1-QB228101.1), "Analyze Literature with Claude AI" showed only "The analysis could not be prepared."
  - The project's `.venv` lacks the `anthropic` package: it is in `requirements.txt` but was never installed there.
  - `ClaudeProvider.count_tokens` imports it during preparation. The `ModuleNotFoundError` is not an `anthropic.APIError`, so it reached the generic handler in `app.py`.
  - No request was ever sent; key, model, and request format were fine.
- **Fix for the error:**
  - The missing SDK is reported before any request: `sdk_installed` in `literature/config.py`, shown by `choose_provider` and the status panel.
  - Token-counting errors the paid request would also hit (rejected key, unknown model, invalid request) stop with a specific message; temporary outages fall back to an approximation.
  - Any other preparation error names its type, says nothing was billed, and logs the details on the server.
  - **To use Claude locally, install the requirements in `.venv`.**
- **Evidence eligibility (`literature/discovery.py`, new `literature/chronology.py`):**
  - A publication is analyzed for a sample only when it was not published before the sample's collection.
  - A study-level link (shared BioProject or dataset identifier) also needs the sample named in the text, for example in a sample-to-site table.
  - Each candidate row records `Chronology` and `Eligibility_Reason`; the page separates unique, direct, study-only, naming, usable, and unusable publications.
- **No evidence:** when no article is usable, nothing is sent to Claude. The page shows "One study-associated publication was found, but it does not provide verifiable collection information for these samples. No metadata changes are proposed…" and the analysis can continue with the original metadata.
- **Discovery:** nuccore → BioSample → PubMed links count as sample-specific. BioSample accessions from NCBI records are added as sample aliases (the real records gained, e.g., SAMN63447205).
- **Synthetic literature example:** SYNPUB-007 (study-level, names no sample) is no longer sent, so 19 instead of 21 facts are produced. Association scoring in `literature/evaluation.py` still counts it as correctly found.
- **Tests:** `tests/test_literature_eligibility.py` (16).

## Metadata readiness on the review step (November 2, 2026)

- **What it shows:** `genomics/metadata_readiness.py` adds a "Metadata readiness for: <analysis>" table to Step 2. Requirements depend on the chosen or detected mode (individual, pooled, combined). Each field is labeled:
  - **Required** (the analysis cannot run without it);
  - **Recommended** (one named step is limited or unavailable);
  - **Optional** (context or extra precision).
- **Fields by mode:**
  - *Individual:* required are the identifier, the pathogen (can be typed on the page), and the genomic region (or confirm complete genomes); dates and locations are recommended.
  - *Pooled:* the variant-call columns are required; dates, sites, and genome coverage are recommended.
  - *Combined:* adds the shared reference FASTA, which is required only for the cross-source comparison.
- **What each row shows:** status (Detected / Partly filled / Imprecise / Missing / Needs your choice) with counts, the calculation it feeds, and the source columns kept unchanged. Each gap gets a specific sentence about what it affects, for example: "Collection locations are only available at the state level. Genetic comparisons can proceed, but precise geographic-distance screening cannot be performed".
- **Blocking and labels (display only; no calculation changed):**
  - The run is still blocked only by the existing review problems (missing identifier, pathogen, or variant columns; column choices).
  - The review plan now says "Unavailable" instead of "Skipped", and gives specific temporal and geographic reasons.
  - Analyses of other modes without files are hidden.
- **Sprint 4:** with the toggle on, the table adds "Could Sprint 4 help?" (dates, locations, coordinates, sample source). After enrichment, `gap_resolution` in `literature/preanalysis.py` re-lists the gaps from the derived metadata: resolved, improved but still limited, or unresolved.
- **Tests:** `tests/test_metadata_readiness.py` (9), using complete, partly filled, and sparse versions of the synthetic NCBI fixture.

## NCBI Virus geographic columns (November 2, 2026)

- **Bug:** NCBI Virus metadata (Geo_Location "USA: SD", Country "USA", and a column named "USA" holding the state) stopped at "Choose which column is Location_Name (Geo_Location, Country)". The cause was that `country` was listed as a Location_Name alias in `genomics/column_detection.py`, so two candidates looked ambiguous.
- **Fix:** Country and State_Province are now separate detected fields. `genomics/location_hierarchy.py` treats the columns as one hierarchy:
  - Geo_Location is used as given and never overwritten;
  - blank rows are built from country/state and labeled with their level and source;
  - contradictory rows are listed under "Location values to review".
- **Precision:** country-, state-, and conflict-level rows are labeled too coarse for the kilometer cutoff and are not sent to GeoNames. Genetic and date comparisons still run, and such pairs are reported as unresolved geography.
- **No change** to distance, date, or cluster calculations.
- **Tests:** `tests/test_ncbi_location.py` (15) with a synthetic fixture copying the 17-column structure: `tests/fixtures/ncbi_virus_structure_*`.

## Sprint 4 from the upload page, before the analysis (November 2, 2026)

- **Toggle.** The upload page has a toggle, off by default: "Use PubTator + Claude AI to find missing or more precise sample information". There is an optional field for dataset-wide publication or study identifiers. It works in every analysis mode.
- **When it is on:** the review step lists metadata gaps without sending anything. "Continue to metadata enrichment" then leads through literature discovery, cost-approved Claude analysis, and review. "Apply approved changes and run the analysis" runs the genomic analysis once, with the approved values:
  - individual values are written into a copy of the prepared metadata, so the pipeline's own GeoNames step resolves place names;
  - pooled values reach the integrated analysis through `pooled_overrides` (site names, dates, publication coordinates).
- **Skipping:** available at every step; the analysis then runs with the original metadata.
- **Code:** `literature/preanalysis.py`, routes `/enrich/start` and `/enrich/run`, `templates/_enrichment_report.html`. Discovery now accepts dataset identifiers as study-level links.
- **Results:** a Sprint 4 report (counts, changed values with evidence, Sprint 2 screening before and after, evidence downloads) and an overview line. Without enrichment there is only a "revisit" link, no empty section.
- **Offline demo:** "Run Integrated Example" with the Sprint 4 box ticked uses `examples/integrated/synthetic_integrated_articles.json` and `synthetic_integrated_mock_responses.json`. IND_D gets Tempe; eight pooled samples get publication site coordinates; one fabricated quote is rejected. Result: one more pair meets the screening, and the shared variant gains a 2.5 km distance.
- **Tests:** `tests/test_upload_enrichment.py` (13). Checkpoint tag before this change: `checkpoint/before-upload-enrichment`.

## Mode-specific upload fields (November 2, 2026)

- The upload card now shows the fields of the chosen mode (`UPLOAD_FIELDS` in `app.py`, a small script in `templates/_upload.html`). Automatic detection keeps the single multi-file box.
- **Individual (Sprint 2):** sequences FASTA and sample metadata (Accession + Pathogen) are required; a reference FASTA is optional.
- **Pooled (Sprint 3):** a processed SNV table or VCF is required. Pooled metadata, per-position depth or A/C/G/T counts, and a reference are optional. Masked sites and variant profiles sit under Advanced inputs.
- **Combined:** the Sprint 2 fields plus the SNV results and a reference are required; the other pooled files are optional.
- Switching modes keeps chosen files but uploads only the current mode's fields, and the page names any kept files that will not be uploaded.
- A role declared through a field is kept when the content fits (`apply_declared_role`). Otherwise the review shows "Conflict — please choose" and runs nothing until a role is chosen. Ambiguous files in the automatic box still need a choice.
- Tests: `tests/test_mode_upload.py` (12). Checkpoint tag before this change: `checkpoint/before-mode-upload`.

## Usability redesign (November 2, 2026, after Sprint 4)

- **Interface only; no scientific calculation changed.** Checkpoint tag before it: `checkpoint/before-ui-redesign`.
- **Home page:** a "How this works" diagram with one lane per sprint. The shared reference genome is shown as the connection, and the steps done before upload (sequencing reads, variant calling) are drawn dashed. Also an optional "What do you want to analyze?" choice (Detect / Sprint 2 / Sprint 3 / Combined), per-analysis file guides with downloadable example files (`/examples/<name>`, an allow-list), an "unlocks" table, and a glossary.
- **Review step:** shows the detected analysis mode and, when a mode was chosen, which files that mode would still need. File-role detection is unchanged.
- **Navigation:** every review, results, and Sprint 4 page has a top toolbar: Back to upload, Back to file review, Results overview, Start new analysis. `POST /back` restores a review only while its staged files still exist. After an analysis, files are already deleted (unchanged privacy rule), so the page explains this and lists the earlier file names to select again.
- **Results:** section 1 is the results overview (mode, major findings from the actual outputs, data quality, collapsed details and limitations). Each section shows its sprint, its key figure, a "How to read this chart" guide (`dashboard_guide.py`), and an "In these results" sentence. Further figures, tables (25 rows, up to 200 more on demand, full CSV download), and limitations are collapsible. The cross-source section follows one shared variant from the reference base to both sources.
- **Sprint 4:** renamed **Sprint 4 — AI Metadata Enrichment**. The article "surveillance context" statements were removed (prompt version `literature-extraction-v2`).
- **System status:** a table of each external dependency, what it does, and what still works without it. Legacy Sprint 3 pages moved into a collapsed footer.
- **Tests:** 296 run; with real MAFFT 294 pass, 2 optional skip. New: `tests/test_ui_navigation.py` (15). The 33-check real-MAFFT matrix still passes. Desktop and 375-px phone layouts were checked in the browser; no horizontal page scrolling.

## Sprint 4 update: AI literature and metadata enrichment (November 2, 2026)

- **What was built:** an optional section 9 of the same dashboard (and optional workflow step 4) that finds publications linked to uploaded samples (NCBI E-utilities record links, metadata publication fields, exact identifiers in PubTator3 text), sends a few relevant passages to a replaceable extraction provider (Claude API with structured outputs; an offline mock; a documented local-model HTTP contract for a future ASU SOL model), verifies every quotation against the retrieved text, and applies only person-approved values to a derived metadata table with a before/after rescreen of dates and distances. Code: `literature/`, `literature_view.py`, `templates/_literature.html`, routes at the end of `app.py`. Full details: `SPRINT4_AI_NLP_INTEGRATION_REPORT.md`.
- **Paid API:** disabled by default (`AI_ENRICHMENT_ENABLED`), no key in the repository or browser, explicit two-step approval with a cost estimate, per-run and per-day caps, and an access code or local-only rule for paid runs. **No live Claude request has been made**; all Claude tests use a fake client or the mock provider. Settings template: `ai_enrichment.example.env` (copy to the gitignored `.env`).
- **Live services checked:** NCBI E-utilities and PubTator3 (free) worked for MN908947.3 → PMID 32015508 with PMC full text. GeoNames and the Docker build are still unverified.
- **Synthetic demonstration and evaluation:** `examples/literature/` (13 invented samples, 11 invented articles, hand-written expected labels, mock responses with planted mistakes); `python validation/evaluate_literature_enrichment.py`. These numbers test the pipeline's checks, not Claude's accuracy.
- **Tests:** 281 run; with real MAFFT 279 pass and 2 optional skip (without MAFFT 278 pass, 3 skip). The 33-check real-MAFFT end-to-end matrix still passes. Checkpoint tag before Sprint 4: `checkpoint/before-sprint4`. Commits are local on `sprint4-integration`; nothing was pushed, merged, or deployed.

## Current state after the Sprint 2 + Sprint 3 integration (November 1, 2026)

This section supersedes the "not yet confirmed" items further below; the rest of the file is the original handoff, kept for history.

- **Branches:** work is on `sprint4-integration` (the branch name predates the work; this is not Sprint 4 work). The branch is pushed and pull request #1 into `main` is open, not merged. `origin/sprint-3` (1c2d59b) was merged in without its generated `tmp/pdfs/` files. A local backup branch `backup/sprint4-before-integration` points at the pre-merge commit `eb72013`. `main` and the Render deployment were not changed.
- **Sprint 3 files that actually existed:** `wastewater/processed_snv.py` (Fontenele Supplementary Table 2 / standardized SNV importer), `environmental_analysis.py` (generic sample-by-feature page), `templates/wastewater.html`, `templates/environmental.html`, and their tests. They remain available as legacy pages.
- **Integration implementation:** `genomics/` package (shared reference-based model, reference-anchored individual variants, paper-based QC, temporal analysis, Yue–Clayton PCoA, cross-source matching, figures, orchestration with provenance), unified dashboard in `templates/index.html`, `integrated_view.py`, `run_integrated_analysis.py`, synthetic paired demo in `examples/integrated/`, optional `offline_workflow/`. See `INTEGRATION_IMPLEMENTATION_LOG.md`, the README section "Integrated analysis (Sprint 2 and Sprint 3)", and the November 1 entry in `SPRINT_3_CHANGELOG.md`.
- **Paper facts verified from the full text** (Europe PMC, PMC8464352) are recorded in `INTEGRATION_IMPLEMENTATION_LOG.md`, including the exact PCoA method (per-position Yue–Clayton summed over positions; zero-depth positions contribute 0).
- **Professor-feedback update (same day):** guided upload → review → run workflow; gzip, TSV, and SNV VCF inputs; FASTQ refused; flexible column detection with explicit ambiguity; optional individual metadata (only `Accession` + `Pathogen` required); 21 dashboard figures with integrated sections first; `ARTICLE_METHOD_VALIDATION.md`. See the newest entry in `SPRINT_3_CHANGELOG.md` and Phase 10 in `INTEGRATION_IMPLEMENTATION_LOG.md`. Local checkpoint tag: `checkpoint/before-guided-workflow`.
- **Data reality:** the local Sprint 2 real dataset is H5N1 influenza HA (segment 4); the paper is SARS-CoV-2 wastewater. They cannot be compared. The real `mmc2.xlsx` (obtained through Europe PMC; not in Git) was imported and matched the paper's reported counts; see `validation/check_fontenele_2021.py`.
- **Tests:** 204 run, 202 pass, 2 optional real-data tests skipped by default. Real MAFFT (`C:\Tools\mafft-win`) was exercised on Windows and passed; GeoNames was not called live and the Docker image was not built.

## Project objective
Build **one integrated genomic surveillance system** from the foundations of two research sprints, then extend that system in Sprint 4. Avoid ending with two unrelated analyses on adjacent web pages.

- **Sprint 2 origin:** Individual/consensus sequence analysis — load FASTA + metadata; validate/accession-match; detect compatible genomic regions; MAFFT-align comparable sequences; measure nucleotide differences/p-distance; combine collection dates with geospatial metadata and threshold-screen potential sample relationships.
- **Sprint 3 origin:** Pooled/wastewater sequencing evidence — initially processed variant/SNV tables and metadata; eventually reflect more of the assigned study's coverage, read depth, allele frequencies, variant-calling quality, temporal variation, and population-structure analyses.
- **Sprint 4 integration target:** One reference-aware variant and metadata system that can compare compatible variants from individual consensus sequences with allele-frequency information from pooled samples. Outputs should include article-grounded quantitative analyses, temporal and geographic visualizations, population structure (where data permit), and transparent uncertainty.

Use general scientific language in the interface and communications where possible. Keep precise genomic grouping and metadata logic internally.

## Important research references
- Sprint 2 conceptual article: Gardy & Loman (2018), *Towards a genomics-informed, real-time, global pathogen surveillance system*, Nature Reviews Genetics. DOI https://www.nature.com/articles/nrg.2017.88
- Sprint 3 study: Fontenele et al. (2021), *High-throughput sequencing of SARS-CoV-2 in wastewater provides insights into circulating variants*, Water Research: https://www.sciencedirect.com/science/article/pii/S0043135421009040
- Related repository shared by instructor/user: https://github.com/ZooPhy/radx . Treat as a reference implementation, not guaranteed current or directly deployable.
- Before reproducing a figure or analysis, review the study's methods and supplementary tables; verify actual columns, reference coordinates, and data coverage. Do not invent paper-derived results.

## Last confirmed Sprint 2 architecture (from discussion)
A Flask app served through Docker + Gunicorn on Render Free, with `GEONAMES_USERNAME` configured in the host environment and MAFFT installed as a Linux executable.

Likely modules:
- `validate_inputs.py`: FASTA/metadata loading, standardize headers/accessions, data validation
- `genomic_regions.py`: identify region from metadata or FASTA header; explicit complete-genome option with safeguards
- `align_sequences.py`: MAFFT per compatible grouping
- `compare_genetic_distances.py`: pairwise p-distance, percent identity, ignoring ambiguous/gap positions
- `geocode_locations.py`: GeoNames lookups when usable locations lack coordinates
- `add_temporal_geographic.py`: days apart, Haversine distance where justified, precision flags
- `detect_potential_clusters.py`: thresholded related pairs and connected components
- `run_pipeline.py`: orchestrator
- `app.py`, `templates/index.html`: upload UI, summaries, synthetic examples, diagnostics
- `system_diagnostics.py`: on-demand two-sequence MAFFT test + GeoNames city lookup, cached briefly
- `Dockerfile`, `requirements.txt`, `tests/`, `examples/`, `README.md`

**Previously observed validation:** Windows local MAFFT v7.526 aligned real sequences; on hosted Render diagnostic MAFFT v7.505 and GeoNames both returned *Operational*. The hosted application ran a real NCBI-style FASTA/CSV test and a synthetic example. This validates the previous version, not necessarily the new checkout.

**Example baselines from previous screenshots (do not hard-code as universal expected results):**
- Publicly sourced individual-sequence demonstration: 10 sequences; 45 pairwise comparisons; 22 genetic+temporal candidate pairs involving 9 samples; 0 strict geographically qualified clusters because all available sampling places were state-level.
- Explicitly synthetic example: 8 sequences; 28 comparisons; one GeoNames-resolved location; 4 strict related pairs; 2 potential clusters.
- Earlier Codex snapshot reported 87 automated tests passing before the later Sprint 3 additions; tests used mocks for certain external dependencies.

## GeoNames and MAFFT details
- GeoNames endpoint: `https://secure.geonames.org/searchJSON` (avoid previously problematic `api.geonames.org` endpoint). Use TLS certificate verification. Username comes from `GEONAMES_USERNAME`; never commit secrets.
- State-only geocoding is too coarse for a 50-km sampling-location threshold; do not interpret state centroid or equal state label as a 0-km collection separation. Keep genomic+temporal candidates separate from three-criteria clusters.
- MAFFT aligns homologous consensus/individual sequences, grouped by compatible region. Do not use MAFFT as a general variant caller for mixed short reads. Windows `mafft.bat` and Linux `mafft` have distinct invocation details.
- Current web prototype reportedly limits upload requests to roughly 10 MB; deployed free service is resource-constrained. Confirm code and hosting settings.

## What is NOT yet confirmed
- Exact latest Sprint 3 file names/modules, whether they were pushed to GitHub, and extent of their processed wastewater SNV validation/plots.
- Which study supplementary tables have actually been loaded versus just mentioned.
- Whether current files support consistent reference genome/variant coordinates across evidence sources.
- Whether the new Claude checkout has all synthetic fixtures, tests, scripts, and unpushed working changes.
- Whether the later application state differs from previous recorded screenshots/results.

**Do not fabricate these details; inspect the repository and ask for absent files only when necessary.** A GitHub download/clone may omit ignored local sample FASTA/CSV, generated outputs, environment configuration, and work that was never pushed.

## Professor feedback / priority for Sprint 4
The current temporary Sprint 3 addition is not sufficiently integrated. The instructor recommended incorporating **more of the actual research methods and evidence**, including variant calling / variant-level quality, temporal structure, PCoA/ordination as appropriate, and richer scientific graphics. Reuse the existing backbone rather than replacing it. The project must ultimately tell one story: how genomic observations from two kinds of sampling evidence relate across compatible reference positions, time, and geography.

Important distinctions:
- The study's population-structure analysis uses **PCoA** on nucleotide-frequency-derived distances, not necessarily standard PCA. Check the method and underlying data before promising exact reproduction.
- Pooled-sample alternate allele fraction is not directly comparable as a measurement to a binary individual-consensus allele call; overlap can be assessed only under compatible reference/quality assumptions.
- A missing variant entry may mean no variant call or insufficient coverage; it must not automatically be treated as 0%.
- Different reference builds, segments, nonoverlapping regions, or incompatible alignment coordinates cannot be overlapped without explicit harmonization.
- More figures should be driven by genuine inputs and transparent denominators, not visually impressive but scientifically invalid results.

## Recommended staged implementation
**First session — audit only:** Identify the actual modules, dependencies, routes, data schemas, test status, workflow gaps, and uncommitted changes. Produce a precise small-step integration plan without editing.

**Next feature 1 — common reference-aware data model:** Create robust import/adapters to carry metadata and standardized variant events for both single-sequence and pooled evidence. Preserve original analysis and explicit source/provenance. Build a tiny compatible synthetic paired demonstration set. Test reference mismatch behavior.

**Next feature 2 — research-grounded analyses:** Work toward variant quality/depth, genome-position frequency, time-series visualization, coverage-aware sample comparisons and PCoA when available. Validate against paper methods/available supplements and label approximations.

**Next feature 3 — unified dashboard:** Replace/augment disconnected pages with a shared workflow: one analysis submission/session; genomic evidence-by-source; compatible variant overlap; temporal and geographic patterns; qualified clustering and caveats. Do not silently convert one data type into another.

**Optional later — raw read ingestion:** BWA/SAMtools/iVar or other supported toolchain can run locally/on research compute; parse resulting QC/variant tables in the web app. Keep large processing out of Render Free unless explicitly benchmarked.

## Branch / deployment practice
- Keep the existing working Render deployment on `main`; develop on a branch such as `sprint4-integration`.
- Never push or deploy to `main` without approval. Run tests first, compare synthetic baseline results, and review staged files for secrets or real sample data.
- Clone the real GitHub repo rather than using Download ZIP if you want to preserve commit history and push new commits.
- Keep local external binaries, API usernames, private data, and bulky sequencing files outside git. Recreate local dependencies from `requirements.txt`, `Dockerfile`, and README.

## Onboarding prompt (suggested)
"Read CLAUDE.md, PROJECT_HANDOFF.md, README.md, the Python modules, tests, Dockerfile and Sprint 3 implementation. First perform a read-only audit. Explain current architecture and input/output schemas, identify what exists and what is missing from the planned shared variant model, compare existing methods to the assigned study, and propose three small, testable integration tasks in priority order. Run non-destructive inspection/tests only if necessary. Do not modify files, push commits, or deploy yet. Report discrepancies between handoff notes and actual code."
