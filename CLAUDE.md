# Project instructions — Unified Genomic Surveillance (BMI 461)

## What this repository is
An undergraduate biomedical informatics class project, currently a Flask application hosted with Docker/Render. It began with an individual-sequence comparison workflow (Sprint 2) and gained a preliminary wastewater/pooled-sequencing variant workflow (Sprint 3). The objective for Sprint 4 is **one integrated analysis**, not two unrelated applications or tabs.

Read `PROJECT_HANDOFF.md`, `README.md`, and the actual code/tests before planning changes. The handoff is historical context, **not** a verified description of the current checkout. When documentation conflicts with code, report the discrepancy and investigate rather than guessing.

## Working constraints
- Preserve the existing working Flask, MAFFT, GeoNames, sequence-comparison, cluster-screening, and deployed synthetic-example functionality.
- Keep the code simple, modular, readable, and explainable in an undergraduate code review. Reuse existing functions; do not replace the project with a new framework.
- Do not change scientific calculations, thresholds, metadata requirements, or validation semantics silently. Document and test any necessary changes.
- Treat individual/consensus sequence calls and pooled-sample allele frequencies as distinct evidence, even when comparing them at shared reference positions.
- Match reference genome identity, coordinate system, genomic region, reference allele, and sample provenance before attempting cross-source variant overlap. A missing variant row does NOT imply a verified zero frequency without coverage information.
- MAFFT is for homologous individual/consensus sequences; it is not a substitute for reference-based short-read alignment or pooled-sample variant calling.
- Avoid treating state-only or similarly coarse geographic labels or their centroids as precise collection coordinates for strict kilometer screening. Maintain geographic-uncertainty explanations.
- Do not claim the analysis confirms outbreaks, transmission, individual strain identity, or lineage assignments. Label prototype outputs appropriately.
- The 2021 study's population-structure method is **PCoA** from nucleotide-frequency-based dissimilarities, not automatically interchangeable with PCA. Reproduce it only if sufficient per-position frequency/coverage information is available; otherwise label an exploratory alternative honestly.
- Use official source methods and the actual available tables when building article-specific analyses; never invent original-study data or silently substitute synthetic examples.
- Keep expensive raw sequencing read processing off the small free web instance unless resource testing supports it. The hosted app should be able to analyze uploaded processed results.
- No secrets, credentials, private specimen records, or large raw data in Git. Read `GEONAMES_USERNAME` from the environment; use HTTPS verification.
- Make code changes on a dedicated development branch. Do not push to the live deployment branch, modify remote infrastructure, or perform irreversible operations without explicit approval.
- Prefer one small independently testable increment per implementation request. Explain any new package or executable dependency.

## Existing code to locate and verify
Core files on `sprint4-integration` (verified November 1, 2026):
- `app.py` (guided `/` → `/review` → `/run` workflow, `/back` navigation, `/examples/<name>` downloads, named-field `POST /` kept for scripts/tests, legacy Sprint 3 routes), `integrated_view.py`, `dashboard_guide.py` (figure reading guides, glossary), `templates/index.html` with partials `_how_it_works.html`, `_upload.html`, `_toolbar.html`, `_status.html`, `_literature.html`, `run_pipeline.py`, `run_integrated_analysis.py`
- Sprint 2 backbone: `validate_inputs.py`, `genomic_regions.py`, `align_sequences.py`, `compare_genetic_distances.py`, `add_temporal_geographic.py`, `detect_potential_clusters.py`, `geocode_locations.py`, `system_diagnostics.py`
- Sprint 3: `wastewater/processed_snv.py`, `environmental_analysis.py` (legacy pages)
- Sprint 2 + Sprint 3 integration: `genomics/` (references, unified_schema, pooled_observations, sequence_variants, quality_control, temporal_analysis, diversity_analysis, variant_matching, figures, integrated_analysis, file_formats, column_detection, input_inspection), `offline_workflow/`, `validation/check_fontenele_2021.py`, `ARTICLE_METHOD_VALIDATION.md`
- Sprint 4 AI literature enrichment (November 2, 2026): `literature/` (config, http_client, identifiers, ncbi, pubtator, sources, discovery, passages, schema, prompts, providers, budget, geography, verification, apply, workflow, evaluation), `literature_view.py`, `templates/_literature.html`, `examples/literature/`, `validation/evaluate_literature_enrichment.py`, `ai_enrichment.example.env`, `SPRINT4_AI_NLP_INTEGRATION_REPORT.md`
- `examples/` (including the synthetic paired demo in `examples/integrated/`), `tests/` (with `tests/fixtures/` and the test-only `tests/fake_mafft.py`), `Dockerfile`, `requirements.txt`, `README.md`, `INTEGRATION_IMPLEMENTATION_LOG.md`, `SPRINT_3_CHANGELOG.md`

## How to verify the current checkout
1. Start by inventorying files, git status/branch, imports, Flask routes, and schema references.
2. Read tests and README. A previous snapshot reportedly passed 87 automated tests before the later wastewater changes; do **not** assume this remains true.
3. In a suitable Python environment run `python -m unittest discover -s tests -v` if the current test suite supports it; report any environment/dependency limitations.
4. Verify MAFFT and GeoNames through the application's existing on-demand diagnostic instead of using an unrelated script as the only check.
5. Do not require real, untracked local datasets just to audit code. The repository should ideally contain explicitly labeled synthetic demo fixtures.

## Desired evolution
1. Audit current Sprint 2 + Sprint 3 code and separate actual implementation from proposed work.
2. Establish a **shared reference-aware genomic variant schema** across individual sequences and pooled samples, with source/provenance, sample date/location, depth/coverage/quality, and reference-position identifiers.
3. Bring in study-specific quantitative analyses where source data support them: variant-calling QC, depth/coverage, variant-frequency/position charts, temporal changes, and appropriately sourced PCoA.
4. Connect the two evidence types in a unified analysis and dashboard showing qualified variant overlap, trends, geographic uncertainty, and limitations — not two isolated result systems.
5. Preserve existing public hosting while gradually adding features; raw FASTQ→variant computation can be an offline/optional ingestion stage.

## Coding and communication style
Use clear function names and small modules; produce a brief file-change summary and reproducible commands. Distinguish observed results, inference, and unverified assumptions. In plain-language explanations prefer general terms such as genomic samples, individual sequences, pooled wastewater, variant profiles, and genomic surveillance; retain precise existing internal data field names where code needs them.

## Session start instruction
The first-session audit is done and the Sprint 2 + Sprint 3 integration has been implemented on `sprint4-integration`. Sprint 4 (AI literature and metadata enrichment) was requested and implemented on November 2, 2026; read `SPRINT4_AI_NLP_INTEGRATION_REPORT.md`. Live Claude requests stay disabled until the user configures them in the server environment; never ask for the API key in chat or put it in files, and never make a paid request without the user's explicit approval. At the start of a session, read `INTEGRATION_IMPLEMENTATION_LOG.md` (verified paper facts, progress, limitations), `ARTICLE_METHOD_VALIDATION.md`, and the top section of `PROJECT_HANDOFF.md`, then run the test suite before changing code. Do not push, merge into `main`, or deploy without explicit approval.
