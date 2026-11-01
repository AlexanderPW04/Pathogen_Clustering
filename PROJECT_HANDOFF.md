# Project handoff — Codex → Claude Code

**Purpose:** Transfer working context for a biomedical informatics class coding project. This handoff was prepared from discussion, Codex-reported summaries, and screenshots. It is **not** a code audit or an assertion that every file exists in the present repository. Verify the current files before editing.

## Current state after the Sprint 4 integration (November 1, 2026)

This section supersedes the "not yet confirmed" items further below; the rest of the file is the original handoff, kept for history.

- **Branches:** work is on `sprint4-integration` (local commits only, not pushed). `origin/sprint-3` (1c2d59b) was merged in without its generated `tmp/pdfs/` files. A local backup branch `backup/sprint4-before-integration` points at the pre-merge commit `eb72013`. `main` and the Render deployment were not changed.
- **Sprint 3 files that actually existed:** `wastewater/processed_snv.py` (Fontenele Supplementary Table 2 / standardized SNV importer), `environmental_analysis.py` (generic sample-by-feature page), `templates/wastewater.html`, `templates/environmental.html`, and their tests. They remain available as legacy pages.
- **Sprint 4 implementation:** `genomics/` package (shared reference-based model, reference-anchored individual variants, paper-based QC, temporal analysis, Yue–Clayton PCoA, cross-source matching, figures, orchestration with provenance), unified dashboard in `templates/index.html`, `integrated_view.py`, `run_integrated_analysis.py`, synthetic paired demo in `examples/integrated/`, optional `offline_workflow/`. See `SPRINT_4_IMPLEMENTATION_LOG.md` and the README section "Integrated analysis (Sprint 4)".
- **Paper facts verified from the full text** (Europe PMC, PMC8464352) are recorded in `SPRINT_4_IMPLEMENTATION_LOG.md`, including the exact PCoA method (per-position Yue–Clayton summed over positions; zero-depth positions contribute 0).
- **Data reality:** the local Sprint 2 real dataset is H5N1 influenza HA (segment 4); the paper is SARS-CoV-2 wastewater. They cannot be compared. `mmc2.xlsx` is not in the repository and was not available in this environment, so the real-table import was verified only with synthetic replicas of its documented headers.
- **Tests:** 164 run, 162 pass, 2 skipped (real MAFFT not installed here; optional real-NCBI test needs untracked files). MAFFT and GeoNames are mocked in automated tests; neither was exercised for real in this environment.

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
