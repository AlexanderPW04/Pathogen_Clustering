# Pre-Sprint 4 Full Audit: Integrated Sprint 2 + Sprint 3 Genomic Surveillance System

**Date:** November 1, 2026  ·  **Branch:** `sprint4-integration` (name predates the work; everything here is Sprint 2 + Sprint 3)  ·
**Audited baseline:** commit `0c68527`, tagged `checkpoint/before-full-audit`  ·  **Scope:** code, science, inputs, interface, reliability.
Sprint 4 (NLP / literature enrichment) was **not** implemented or scaffolded.

How the audit was done:

1. Baseline: inventory, git state, full test run with real MAFFT (204 tests: 202 passed, 2 optional real-data tests skipped).
2. Four independent read-only reviews (Sprint 2 pipeline; reference mapping and importers; QC, statistics, PCoA and matching; Flask, safety and interface), each asked to reproduce suspected defects with concrete inputs and to hand-verify formulas. Every finding used below was then reproduced or checked again before fixing.
3. Fixes in priority order, each with a regression test (`tests/test_audit_regressions.py`). The regression file was run against the pre-audit checkpoint: **19 of its 20 tests fail there** (13 failures, 6 crashes) and all pass now; the remaining one is a known-answer check that was already correct.
4. An end-to-end validation matrix through the web application with **real MAFFT** (33 checks, all passing), the real NCBI H5N1 dataset, and the paper's real Supplementary Table 2.

Words used for evidence: **verified** = executed here and the result checked against an independent expectation; **mocked** = executed with a stand-in for MAFFT or GeoNames; **not verified** = not executed in this environment.

---

## A. Executive summary

**What the system does.** A local Flask application (also deployable with Docker/Render) that analyzes two kinds of genomic evidence in one workflow:

- *Individual (consensus) sequences* — the Sprint 2 pipeline: FASTA + metadata validation, genomic-region grouping, MAFFT alignment, p-distance / percent identity, days apart, Haversine distance with geographic-precision rules, and strict "potential clusters" kept separate from "genomic-temporal candidates with unresolved geography".
- *Pooled samples such as wastewater* — the Sprint 3 work: processed SNV tables (standardized, iVar, the Fontenele et al. supplementary table, SNV VCF), Fontenele-style QC applied only to values supplied, per-site frequency trends, depth-versus-SNV regression, and population structure (exact Yue–Clayton PCoA when A/C/G/T counts exist, otherwise a labeled exploratory ordination).
- *Integration* — individual sequences are aligned to a verified reference genome, so both sources share 1-based reference coordinates; each variant is classified as shared, source-specific, unresolved (coverage missing), or incomparable (reference not verified). A guided upload (upload → review → results) recognizes file roles and columns, and the dashboard presents eight sections with 21 figures, integrated sections first.

**Readiness.** After this audit's fixes, the Sprint 2 + Sprint 3 foundation is **scientifically defensible for a class prototype and technically ready for Sprint 4 to build on**, with the limits listed in Section H (open items) and Section L. No critical defect remains open. Two external checks remain unverified in this environment: a live GeoNames request and a Docker build.

**Main strengths**

- The core numerical methods were independently re-derived and are correct: p-distance and identity, Haversine, Yue–Clayton per-position dissimilarity summed over positions with zero-depth positions set to 0 (as the paper states), Gower PCoA (B = −½ J D² J), least-squares regression with a bootstrap band and permutation p-value, and the vectorized frequency trends (300 random cases matched an independent loop).
- The "missing ≠ 0%" rule holds through imports, the observation matrix, figures, downloads, and cross-source matching.
- Reference compatibility is strict: different labels are never assumed equal; Ref bases are checked against the reference.
- The synthetic paired example reproduces all 8 planted cross-source answers with **real MAFFT**; the real H5N1 data reproduce the documented Sprint 2 baseline; the real Fontenele supplementary table reproduces the paper's record, sample, state, per-sample, and D614G counts exactly.
- Upload handling is safe: user file names never reach the filesystem, Jinja autoescaping covers every user-supplied string, downloads are dictionary lookups with random tokens.

**Main weaknesses found (now fixed)** — see Section H for all 47 entries.

- Sprint 2 lost sample IDs such as `001` or `NA` when reading its own pairwise CSVs (crash or mismatch).
- Coarse or placeholder locations (`USA: California`, `USA`, `not provided`) counted as "same location" and produced strict clusters when GeoNames was off.
- A pair sharing one called base had p = 0 and could form a cluster.
- Mixed day-first/month-first dates were parsed silently in two conventions.
- iVar tables with a lowercase `sample` column were merged into one sample; reported 0% rows counted as detections.
- QC-excluded pooled calls were counted as reports in incomparable rows and as assessable non-detections in profile overlap; state-level pooled coordinates produced kilometre distances.
- Twenty small gzip uploads could expand to more than 1 GB on disk; a few malformed inputs produced HTTP 500 errors.

**Remaining issues (not critical)** — memory use of very large pooled tables on a 512 MB host, country-only labels outside the USA are not recognized as coarse without GeoNames, an unverified live GeoNames call and Docker build, and the paper methods that cannot be reproduced from processed tables (read mapping, base-quality filtering, GISAID novelty, PANGOLIN lineages, the exact PCoA for the real data).

---

## B. Architecture review

### Pipeline

```text
                         Browser                                   Command line
            GET / (upload card + examples)                 run_integrated_analysis.py
                  │ POST /review (files)                   run_pipeline.py (Sprint 2 only)
                  ▼                                                    │
 app.stage_uploads: save → genomics.file_formats.plain_copy            │
 (gzip, one 60 MB budget per request) → detect_kind → inspect_file     │
                  │                                                    │
 genomics.input_inspection: suggest_roles (evidence only) →            │
 review_dataset (columns, ambiguity, units, IDs, reference,            │
 completeness, analysis plan)  → review page (Step 2)                  │
                  │ POST /run (token, roles, column choices, settings) │
                  ▼                                                    │
 prepare_inputs: standardized copies (column renames, FASTA            │
 annotations, stated pathogen, combined VCF/iVar files)                │
                  ▼                                                    ▼
 app.run_analyses ──► run_pipeline (Sprint 2)                  same functions
   validate_inputs → geocode_locations → align_sequences (MAFFT)
   → compare_genetic_distances → add_temporal_geographic → detect_potential_clusters
                  │ validated samples, MAFFT path, pairwise table
                  ▼
 genomics.integrated_analysis.run_integrated_analysis
   references.load_reference_fasta (SHA-256)
   sequence_variants.build_individual_profiles (MAFFT to reference → 1-based calls)
   wastewater.processed_snv.load_processed_snvs (+ file_formats VCF conversion)
   pooled_observations (shared tables, coverage from depth)
   quality_control.apply_pooled_quality_control (paper / exploratory)
   temporal_analysis (observation matrix, trends, regression)
   diversity_analysis.population_structure (Yue–Clayton PCoA / exploratory)
   variant_matching.compare_evidence (+ per-sample cells, profile overlap)
   figures.* (21 PNG figures)  → provenance, analysis status, integration summary
                  ▼
 integrated_view.build_dashboard (8 sections) → templates/index.html
 app.store_result_files (in-memory CSV/PNG cache, 10 min) → /results/<id>/<file>
```

Legacy Sprint 3 pages remain separate and unchanged: `/wastewater` (processed SNV summary) and `/wastewater-data` = `/environmental` (feature-measurement tables). They do not interfere with the integrated workflow (verified by the Flask review: all render, run their examples, and serve downloads). `/system-status` runs the on-demand MAFFT and GeoNames diagnostics; `/health` returns `OK`.

### How Sprint 2 and Sprint 3 are connected

The connection is data-level, not presentational. Sprint 2's validated, geocoded sample table and MAFFT path are handed to `build_individual_profiles`, which aligns each sequence to the reference FASTA and emits variants in the same `Variant_Key = Reference_ID:Position:Ref>Alt` schema (`genomics/unified_schema.py`) as the pooled importer. `compare_evidence` then evaluates every variant in every sample of both sources using each source's own evidence (consensus base calls versus allele frequencies with depth), and the dashboard's cross-source section, overview, temporal timeline, and geography figure all draw on both sources at once.

### Inventory (tracked files)

28 Python modules (about 11,300 lines in the application, plus 24 test files), 3 templates, the synthetic examples (`examples/`, `examples/integrated/` with a deterministic generator — regenerated byte-for-byte during the audit), `offline_workflow/` (template + tested converters), `validation/`, Dockerfile, requirements, and documentation. No circular imports. Entry points: `app.py` (Flask), `run_pipeline.py`, `run_integrated_analysis.py`, `validation/check_fontenele_2021.py`, `examples/integrated/make_integrated_example.py`, `offline_workflow/convert_offline_outputs.py`.

---

## C. File-by-file code review

Line numbers refer to the audited-and-fixed working tree (commit following `0c68527`).

| File | Purpose / main functions | Inputs → outputs | Issues found | Improvements made | Remaining concerns |
|---|---|---|---|---|---|
| `validate_inputs.py` | `validate_metadata` (216), `match_fasta_to_metadata` (460), `load_fasta`, `load_and_validate_inputs` | FASTA + CSV/XLSX → validated sample table, errors, warnings | Placeholder words were errors (S2-5); mixed date conventions silent (S2-6); `Nuc_Completeness` header spelling (S2-7); version mismatch silent (S2-8); FASTA BOM (S2-12) | `MISSING_VALUE_WORDS` (67), `ambiguous_day_month` (194), completeness aliases (60), `version_mismatches` (737), `utf-8-sig` (436) | Dotted non-version IDs (`S1.1`, `S1.2`) are treated as versions of `S1` and reported as duplicates (conservative, documented) |
| `genomic_regions.py` | Region labels from metadata and FASTA descriptions | text → region info | None found; complete-genome confirmation cannot override contradictions (verified) | — | Region words outside its vocabulary need a metadata column |
| `align_sequences.py` | `find_mafft`, `run_mafft` (30), `align_samples` | grouped samples → aligned sequences | No timeout (S2-4); `.bat` with spaces in both paths failed (S2-9) | `MAFFT_TIMEOUT_SECONDS = 150` (15) → readable error; `.bat` input passed by plain name from its folder | — |
| `compare_genetic_distances.py` | `compare_sequences` (39), `read_pair_table` (29) | aligned sequences → pairwise CSV | IDs lost on re-read (S2-1) | `read_pair_table` keeps IDs as text; used by every reader | — |
| `add_temporal_geographic.py` | `haversine_km` (15), `enrich_pairs`, `run_stage` | pairs + metadata → days, km, precision | Time-of-day order dependence (S2-10); missing lat/lon columns crash (S2-11) | Calendar-day difference (93); optional columns added | — |
| `detect_potential_clusters.py` | `find_related_pairs` (40), `build_clusters` | enriched pairs → related pairs, clusters | No minimum overlap (S2-3) | `MIN_COMPARABLE_POSITIONS = 20` (18), `enough_overlap` (87) — also used for app candidates | 20 is a prototype safeguard, not a published standard |
| `geocode_locations.py` | `is_state_only_location` (30), `fill_missing_coordinates` (51), `check_geonames` | locations → coordinates + precision | Coarse US labels and placeholders not recognized (S2-2); GeoNames class L areas not coarse (S2-13) | US country / state names, placeholder words, class L and CONT/RGN/AREA codes coarse (146) | Country-only labels outside the USA are coarse only when GeoNames flags them; no cap on lookups per run |
| `run_pipeline.py` | Six-stage runner | FASTA + metadata → four CSVs | — | Prints metadata notes | — |
| `system_diagnostics.py` | MAFFT and version checks | — | `.bat` quoting (S2-9) | Same fix as `run_mafft` | — |
| `genomics/references.py` | `load_reference_fasta`, `references_compatible` (89) | FASTA → reference record (SHA-256) | None; different labels never equal without identical sequences (verified) | — | The app registry holds only the uploaded reference, so differently labeled records are always incomparable (safe) |
| `genomics/sequence_variants.py` | `map_alignment_to_reference` (61), `build_individual_profiles` | sequence + reference → calls, variants, indels | One sample's alignment failure aborted the run; overhangs labeled insertions (IMP-12) | Failed sample marked unmapped (197); "Unaligned overhang" (311) | Mapping safeguard (≥ 20 callable positions, ≥ 90% identity) is a prototype rule |
| `genomics/pooled_observations.py` | shared pooled tables, depth lookup, coverage | SNV table, metadata, depth → tables | None found; unlisted positions are "unknown" (verified) | — | — |
| `wastewater/processed_snv.py` | `load_processed_snvs`, `read_table` | CSV/TSV/XLSX/VCF → standardized SNVs + warnings | Lowercase iVar sample column (IMP-1); 0% rows as detections (IMP-2); mixed-format label (IMP-8); no frequency/count cross-check (IMP-9); conflicting duplicates silent (IMP-10) | `SAMPLE_HEADERS` (67); 0% rows set aside (269); `COMBINED_FORMAT` (73); warnings | Repeated rows keep the first value (now stated) |
| `genomics/file_formats.py` | `plain_copy`, `detect_kind`, `vcf_to_snv_table`, `allele_measurement` (233) | files → plain copies, VCF rows | Corrupt gzip crash (IMP-4); infinite values (APP-4); missing FORMAT fields silent (IMP-5); AD vs DP denominators (IMP-6); INFO/AF in sample-column VCFs (IMP-7) | `zlib.error` caught (119); non-finite → missing (225); set-aside reason; depth = allele-count total; INFO/AF only for sites-only files (338) | VCF FILTER is recorded, not applied |
| `genomics/column_detection.py` | `detect_columns`, `resolve_choices`, `apply_mapping`, `value_warnings` | headers → mappings | Auto choice labeled as user choice (IMP-11) | `resolve_choices(..., chooser)` | — |
| `genomics/input_inspection.py` | `inspect_file`, `suggest_roles`, `review_dataset`, `prepare_inputs` | staged files → review, standardized inputs | Non-FASTA reference → 500 (APP-2); corrupt XLSX → 500 (APP-3); contig check promised but absent (IMP-3); combined-file warnings lost (IMP-8) | Reference must be FASTA (456); any unreadable file refused on its own (280); contig mismatch blocks the run; per-file notes kept (795) | — |
| `genomics/quality_control.py` | `apply_pooled_quality_control` (110), `method_report` | variants + samples → QC columns, report | Masked list for another reference reported Pass (SCI-4) | "Not evaluable" for other references (131) | Paper mode keeps samples whose coverage is unknown (flagged "Not evaluable"); iVar defines PASS at p ≤ 0.05 vs our p < 0.05 (matters only at exactly 0.05) |
| `genomics/temporal_analysis.py` | `observation_matrix`, `frequency_trends` (94), `depth_regression` (190) | matrix → trends, regression | Permutation ties dropped by rounding (SCI-6) | Tolerance (223) | Observation matrix memory grows as variants × samples (Section H, PERF-1) |
| `genomics/diversity_analysis.py` | `yue_clayton_dissimilarity`, `paper_distance_matrix` (69), `principal_coordinates` (104) | A/C/G/T counts → distances, PCoA | None (verified by hand) | — | float32 storage (≤ 3.2e-7 difference); axis signs arbitrary; diversity is per-site Gini–Simpson averaged over covered positions |
| `genomics/variant_matching.py` | `compare_evidence` (122), `compare_signatures` (351), `usable_coordinates` (86) | both sources → categories, cells, overlap | QC-excluded counted as reports (SCI-1) and as non-detections (SCI-2); state-level coordinates used for km (SCI-3) | Fixed at 271, 351, 86 | — |
| `genomics/integrated_analysis.py` | orchestration, summaries, status, integration facts | paths → result dict | Individual-only dates missing (APP-6); no variation figure for individual + reference (APP-7); partial coverage "Completed" (SCI-5) | `build_summary` (745); individual-only position summary; "Partly completed" (635) | — |
| `genomics/figures.py` | 21 figure functions | tables → PNG | Very wide figures unreadable when scaled (UX-1) | Wide figures scroll at full size (template + `integrated_view.png_width`) | Matplotlib draws one patch per cell; large heatmaps are slow but capped at 40 variants |
| `integrated_view.py` | `build_dashboard`, `download_files` | result → 8 sections | Links to skipped files (APP-9); "0 shared" with incomparable references (APP-8); two "Unresolved" labels (APP-13) | Links only for cached files; incomparable count shown; `short_category` (79) | — |
| `app.py` | routes, staging, cache | requests → pages | Disk blow-up (APP-1); 500s (APP-2/3/4); uploads kept after run (APP-5); orphaned folders (APP-12); confirmation demanded without sequences (APP-14); wording (APP-15) | One decompression budget per request (884); broad, logged error handling in `/review` (949) and `/run` (980); uploads deleted after a run; orphan cleanup in `clean_staged` (848) | Settings checks are partly duplicated between `index()` and `check_settings`; `INPUT_LABELS` and `ROLE_LABELS` overlap |
| `templates/index.html` | three-step page, section macro | context → HTML | Sprint 2 counts shown twice; misleading section-8 message when Sprint 2 failed (APP-10, APP-13) | Duplicates removed; correct message | — |
| `run_integrated_analysis.py` | CLI | paths → CSV/PNG folder | Unused import (found during the audit) | Removed | `--output-folder` equal to the input folder raises `SameFileError`; `--example` overrides other inputs |
| `environmental_analysis.py`, `templates/wastewater.html`, `templates/environmental.html` | legacy Sprint 3 pages | — | None found | — | Kept for comparison; not part of the integrated analysis |

Duplication noted but intentionally not refactored (it is stable, tested, and a refactor would add risk before Sprint 4): table readers exist in `validate_inputs.load_metadata`, `pooled_observations.read_any_table`, `processed_snv.read_table`, and `file_formats.read_table`; header normalization exists in `pooled_observations.standard_name` and `column_detection.squash`.

---

## D. Scientific validation

All expected values below were computed by hand or by an independent loop, not by the function under test.

**Genetic distance (Sprint 2).** p = differences / comparable positions, identity = (1 − p) × 100; only A/C/G/T versus A/C/G/T positions count, so gaps, N, and IUPAC codes are excluded. Each unordered pair appears once (i < j). Known answer: `ACGTN-` vs `ACGAAA` → 4 comparable, 1 difference, p = 0.25. *Change:* a pair now needs at least 20 comparable positions to count toward relatedness (before, `A` + 27 N versus a full sequence gave p = 0 and a strict cluster). Unequal original lengths are handled by alignment; MAFFT output IDs, lengths, and characters are validated, and output order does not matter.

**Reference coordinates.** Verified with hand-made alignments and real MAFFT: insertions sit between reference positions and never shift coordinates; internal gaps are deletions (D); terminal gaps are "outside the sequenced span" (.); N/IUPAC are uncallable; a substitution at sample index 120 of a sequence starting at reference base 21 maps to position 141. Overhangs beyond the reference ends are now labeled unaligned, not insertions. A sample that fails alignment is now reported as unmapped instead of stopping the run.

**Variant import and frequencies.** Positions must be positive integers; Ref/Alt one of A/C/G/T and different; frequencies in 0–1 (percentages block the run on review); depth positive integers; alternate reads ≤ depth; sample IDs must exist in pooled metadata. *Changes:* 0% rows are set aside (they are not variant calls and are not coverage evidence); frequency/count disagreements > 0.05 are reported; conflicting duplicate rows are reported; VCF depth equals the allele-count total when frequency is computed from allele counts; INFO/AF is used only for sites-only VCFs; VCF lines without sample fields are counted as set aside.

**Missing ≠ zero.** Traced through: importer (blank depth stays blank), `prepare_inputs`/`combine_variant_tables` (blank → blank, never 0), depth lookup (−1 = unknown), observation matrix (`Coverage unknown`, NaN frequency), trajectory and heatmap figures (unknown drawn in a separate strip/hatch), downloads (blank), and cross-source matching (unknown coverage → Unresolved). Verified by the import review and by unit tests.

**Quality control.** Hand checks: coverage 90.0% fails (> 90% required); depth 20 passes, 19 fails; p = 0.05 fails, 0.0499 passes; base quality is "Not verifiable"; unverifiable checks are never reported as passed. Paper mode excludes failures from analyses but keeps them in downloads; exploratory mode flags them. *Changes:* a masked-site list for another reference now gives "Not evaluable"; QC-excluded calls are never counted as reports or non-detections anywhere; the analysis-status table says "Partly completed" when coverage is known for only some samples.

**Temporal and geography.** Days apart are calendar days (time of day no longer changes the result); year-only dates are not given a day; mixed day/month conventions are rejected with an explanation; placeholder words are missing values. Haversine: one degree of latitude = 6371 × π / 180 km (verified to 6 decimals). Coordinates are validated to ±90/±180. State/country-level labels, GeoNames administrative/area results, and missing locations never satisfy the kilometre cutoff and never produce cross-source kilometre distances. Slopes are least squares on actual day offsets (irregular intervals respected); hand slope 15:T>C North = 0.200 per week (exact).

**Cluster detection.** Thresholds are inclusive (≤). Missing dates exclude a pair; kilometres take priority over the same-name fallback; coarse locations give neither a distance nor a same-location match. Clusters are connected components within pathogen + region groups (single linkage, documented). Strict clusters and genomic-temporal candidates with unresolved geography stay separate; neither is called transmission.

**Cross-source matching** (`compare_evidence`). For each variant on the verified reference: individual carriers (exact Variant_Key), reference calls, other alleles, or not callable; pooled reports (QC-included only), adequate depth without a report, low/no coverage, unknown coverage, or QC-excluded. Shared = carriers and reports; Individual only requires adequate pooled depth; Pooled only requires callable individual reference calls; otherwise Unresolved; records whose reference label differs or whose Ref base disagrees are Incomparable. All 8 planted answers of `examples/integrated/make_integrated_example.py` match in both QC modes, with mocked and with real MAFFT; every per-variant count matched a recount from the per-sample cells. Date-window pairs (10 for 15:T>C; 2 for 42:T>A) and the date offsets (−1 and −17 days) were checked by hand.

**PCoA.** Yue–Clayton: [1,0,0,0] vs [.5,.5,0,0] = 0.5; identical = 0; [.6,.4,0,0] vs [.2,.2,.6,0] = 0.7368; a 3-sample matrix with a zero-depth position gives 0.5 / 1.5 / 1.5 with that position contributing 0, as in the paper. Matrix symmetric with zero diagonal. Gower PCoA double-centres D²; on the 3-4-5 triangle the eigenvalues are 12.964 and 3.703 (trace 16.667) and the coordinates reproduce the distances exactly; a non-Euclidean example gives eigenvalues 12.5, 0.5, −5.5 with one negative eigenvalue reported (ratio 0.423). Explained variation uses the positive-eigenvalue sum and is labeled that way. Sparse SNV-only tables never use this method; they get a Jaccard ordination titled "EXPLORATORY – not the paper's Yue–Clayton PCoA".

**Depth-versus-SNV regression.** Hand least squares for x = 10..50, y = [1,3,2,5,4]: slope 0.08, intercept 0.6, R² 0.64 (matched); the bootstrap 95% band contains the fit; p = (extreme + 1)/(n + 1); samples without mean depth or excluded by QC are dropped; fewer than 10 samples are labeled descriptive only. *Change:* tied permutation slopes are now counted (x = [0.1,0.2,0.3], y = [0.1,0.7,0.3] has exact p = 1; before the fix it returned 0.834).

---

## E. Article comparison (Fontenele et al. 2021, Water Research 205:117710)

The detailed per-method table with paper evidence is in `ARTICLE_METHOD_VALIDATION.md`; this summary reflects the audit's re-check.

| Published method | Current implementation | Fidelity | Required input | Verification | Remaining limitation |
|---|---|---|---|---|---|
| Fig 1 / §2.4 coverage: > 90% of the genome at ≥ 10x (Q20) | Coverage filter on supplied or computed coverage; coverage bars with the 90% line | Partial (threshold faithful; Q20 not verifiable) | Coverage in pooled metadata or a complete per-position depth table | Boundary tests; synthetic example excludes WW_S4 (78%) | The real supplementary table has no coverage column |
| §2.4 iVar calling: depth ≥ 20, Q20, no frequency cut-off | Depth ≥ 20 applied; no frequency cut-off; Q20 reported as not verifiable | Partial | Read depth per record | Unit tests; real table 7,961 of 7,973 retained | 12 published records are below 20x |
| §2.4 Fisher's exact p < 0.05 | Applied to supplied p-values or iVar PASS flags; never recalculated | Partial | iVar PVAL / PASS | Unit tests | Not evaluable for the published table |
| §2.4 De Maio problematic sites | Applied to an uploaded list (reference-aware) | Partial | Masked-site table | Unit + regression tests | The historical list is not bundled |
| §2.4 BWA-MEM mapping, iVar trim | Offline template only | Not reproduced in the app | Raw reads | Converters unit-tested; template not run | Processed uploads are not equivalent to raw-read processing |
| Fig 2A SNVs per sample | Bar chart by site, QC exclusions hatched | Faithful | Pooled SNV table | Real table: 24–793 per sample (matches the paper) | — |
| Fig 2B SNVs vs mean depth with 95% CI | Least squares, bootstrap band, permutation p | Partial (model not stated in the paper) | Mean depth per sample | Hand-verified statistics | Needs the paper's Table 1 depths as metadata |
| Fig 3 SNVs by genome position | Position plot: samples per variant, Alt base, source, multi-allele rings, published novelty stars | Partial (all SNVs; clinical-follow-up flag not drawn) | SNV table and/or individual profiles | Real table drawn (5,669 SNVs) | Annotation bands span annotated records only |
| Fig 4 lineage support / clinical comparison | Documented-profile overlap (never a lineage call); cross-source timeline with individual sequences | Exploratory | Profile table; individual sequences | Synthetic example | No GISAID / PANGOLIN |
| Fig 5 Yue–Clayton PCoA | Exact per-position Yue–Clayton summed over positions, zero-depth = 0, Gower PCoA, eigenvalue panel | Faithful when A/C/G/T counts exist | Per-position A/C/G/T counts | Hand-computed matrices and eigenvalues | Real data need counts regenerated from SRA reads |
| Fig 5 without counts | Jaccard ordination of reported SNV sets | Exploratory (labeled) | SNV table | Real table: 52 samples | Axis 1 explains 4.8%; depth can drive distances |
| Novelty vs GISAID | Published flag shown only from the source table | Partial (display) | Published table | 589 multi-sample novel SNVs vs the paper's 548 | GISAID comparison not possible |

**Published supplementary data (real).** `mmc2.xlsx` was obtained through the Europe PMC supplementary-files service (PMC's own links use a bot challenge, which was not bypassed) and kept outside Git. The parser recognized the real layout (descriptive first row, published headers, lower-case `REF`), preserved all published annotations as `Published_*` columns, and kept sample IDs as text. Exact matches with the paper: 7,973 records, 52 samples, 10 states, 24–793 SNVs per sample, D614G in all 52, 4 novel SNVs in more than 8 samples. Reported differences: 5,677 vs 5,680 distinct SNVs; 4,368/1,309 vs 4,372/1,308 non-synonymous/synonymous; 589 vs 548 multi-sample novel SNVs; 12 records below 20x (excluded by paper-mode QC, kept in downloads). Re-run after all fixes with identical results. No missing measurement was invented to make a figure reproducible.

---

## F. Input format and metadata compatibility

| Format | Parser | Status | Tested |
|---|---|---|---|
| FASTA `.fasta/.fa/.fna` | Biopython via `validate_inputs.load_fasta`; content detection | Full | Yes (unit, Flask, real NCBI data, real MAFFT) |
| Gzip FASTA | `file_formats.plain_copy` (one 60 MB budget per upload) | Full | Yes (unit, Flask; corrupt gzip and multi-file budget) |
| CSV / TSV / TXT | `file_formats.read_table` and the role loaders | Full | Yes |
| XLSX (first sheet) | pandas/openpyxl | Full (corrupt files refused per file) | Yes |
| Gzip CSV/TSV | as above | Full | Yes |
| Standardized processed SNV table | `processed_snv.load_processed_snvs` | Full | Yes |
| iVar variants TSV (combined or one per sample) | same, iVar branch; any capitalization of the sample column | Full for SNVs (indels set aside) | Yes (incl. lowercase `sample` regression) |
| Fontenele et al. Supplementary Table 2 | same, published branch | Full | Yes (real workbook) |
| SNV VCF / `.vcf.gz` | `file_formats.vcf_to_snv_table` | Partial by design (allele-specific AD, AO/RO, FORMAT/AF, iVar ALT_FREQ; INFO/AF only in sites-only files) | Yes (multi-allelic, GT-only, indel, sites-only, cohort INFO, missing FORMAT, AD/DP, infinite values) |
| Reference FASTA (one record, A/C/G/T/N) | `references.load_reference_fasta` | Full | Yes |
| Per-position depth or A/C/G/T counts | `pooled_observations.load_position_table` | Full | Yes |
| Masked sites, documented profiles | QC and matching loaders | Full | Yes |
| FASTQ, SAM/BAM/CRAM | recognized and refused with a pointer to `offline_workflow/` | Not accepted (by design) | Yes |

**Role detection** suggests a role only with stated evidence (shared IDs with another file, the reference named in a variant table, the word "reference" in a FASTA description, required headers, published or iVar layouts); a one-sequence FASTA without evidence and metadata matching no other file are left for the user to choose. **Column detection** is role-specific, ignores case/spaces/underscores/hyphens/dots, blocks the run when two columns could be the same field, refuses fraction-coverage and percentage-frequency units, and labels automatic choices as automatic. **Analysis-specific requirements:** genetic distance needs IDs, valid sequences, a pathogen (column, FASTA annotation, or stated once), and a region (metadata, FASTA description, or confirmed complete genomes); temporal screening needs full dates; geographic screening needs adequately precise locations; cross-source comparison needs a verified common reference; exact PCoA needs A/C/G/T counts. Missing optional data disables only the affected analysis (verified: FASTA + `Accession,Pathogen` metadata alone runs 10 comparisons with a "No Collection_Date column" note). Essential identifiers and reference checks stay strict.

**Upload safety.** 10 MB per request (413 with a readable message); at most 20 files; one 60 MB decompression budget per request (a 3 × 9 MB-expansion bomb is refused); user file names never used as paths; temporary folders deleted after a run, after 30 minutes, or (for folders left by an earlier process) at the next upload; random 128-bit tokens; download names are dictionary keys (no traversal); every user string autoescaped. Malformed CSV/TSV/XLSX/VCF/gzip, empty files, and protein or FASTQ content produce messages, not 500 errors.

---

## G. Visualization review

Every section now has a purpose statement, summary statistics ("Not available" never shown as 0), figures with captions and PNG downloads, tables with CSV downloads (links only for files actually cached), and an "Uncertainty and limitations" list. Integrated sections come first.

| Section | Figures (data requirement) | Audit result |
|---|---|---|
| 1 Overview | Metadata completeness (any source); integration facts; analyses run/skipped | Fixed: individual-only date range; "Reference FASTA" label; incomparable count shown next to "0 shared" |
| 2 Cross-source evidence | Category counts; variant-overlap matrix; positions by category; timeline (both sources + reference) | Verified against the generator's answers; distinct "Unresolved (pooled/individual)" labels |
| 3 Genomic variation | Variants along the reference (pooled and/or individual) | Fixed: now drawn for individual sequences + reference without pooled data; marks scale for thousands of SNVs |
| 4 Temporal patterns | Sampling timeline (dates); trajectories; first-vs-last changes; missing observations; heatmap (pooled dates) | Irregular dates plotted at true positions; unknown never 0 |
| 5 Geographic context | Precision classes + plain coordinate plot (locations) | Coarse labels never plotted as points; not a map |
| 6 Sequencing quality | Coverage by sample; QC outcomes (pass/fail/not verifiable); depth distributions; callable positions; SNVs per sample; depth vs SNVs | Hand-verified regression; unverifiable checks drawn separately |
| 7 Population structure | PCoA (date colour, site/state markers) + eigenvalue panel | Exploratory method labeled in title and badge |
| 8 Individual sequence relationships | Distance distribution + genetic-vs-temporal plot; relationship network (≤ 40 sequences) | Duplicate counts removed; correct message when the comparison itself fails |

The nine originally reported figures were re-checked for denominators, labels, and missing-value handling. Figures wider than 1,400 px (for example SNVs per sample, missing observations, and the heatmap with 52 real samples) are now shown at full size in a scrolling box instead of being shrunk until unreadable (verified in the browser DOM). Colours follow a colour-blind-checked palette with legends and marker shapes as secondary encoding.

---

## H. Issue register

Severity: C critical, H high, M medium, L low. Status: **Fixed** (with regression evidence), **Documented** (accepted limitation), **Open**.

| ID | Sev | File / function | Description | Evidence | Fix | Verification | Status |
|---|---|---|---|---|---|---|---|
| S2-1 | H | `add_temporal_geographic.run_stage`, `detect_potential_clusters.run_stage`, `app.run_relationship_analysis` | IDs like `001`/`NA` changed on CSV re-read → crash or mismatched IDs | Reproduced ("pair 1 and 2") | `compare_genetic_distances.read_pair_table` | Regression test fails before, passes after | Fixed |
| S2-2 | H | `geocode_locations.is_state_only_location` | `USA: California`, `USA`, `not provided` matched as same location → strict clusters | Reproduced (cluster formed) | US country/state names, placeholder words | Regression | Fixed |
| S2-3 | M | `detect_potential_clusters.find_related_pairs`, `app.find_geography_unresolved_candidates` | One comparable base gave p = 0 | Reproduced | Minimum 20 comparable positions | Regression; four 4-base test fixtures lengthened to 20 bases with identical p-distances | Fixed (documented rule change) |
| S2-4 | M | `align_sequences.run_mafft` | No timeout | Code review | 150 s timeout, readable error | Suite | Fixed |
| S2-5 | M | `validate_inputs.validate_metadata` | `missing`/`NA` dates or coordinates were errors | Reproduced | Placeholder words → missing with a note | Regression | Fixed |
| S2-6 | M | `validate_inputs.validate_metadata` | `03/04` and `13/04` parsed in different conventions | Reproduced | `ambiguous_day_month` error | Regression (US-only dates still accepted) | Fixed |
| S2-7 | M | `validate_inputs` aliases | `Nuc Completeness` spelling bypassed the complete-genome safeguard | Reviewer reproduction | Aliases | Suite | Fixed |
| S2-8 | L | `validate_inputs.load_and_validate_inputs` | Different accession versions matched silently | Reviewer reproduction | Visible note | Suite | Fixed |
| S2-9 | L | `align_sequences.run_mafft`, `system_diagnostics.check_mafft` | `.bat` with spaces failed | Reproduced | Plain input name from its folder | Regression + real MAFFT in a spaced folder | Fixed |
| S2-10 | L | `add_temporal_geographic.enrich_pairs` | Days apart depended on time of day and order | Reviewer reproduction | Calendar days | Regression | Fixed |
| S2-11 | L | `add_temporal_geographic.enrich_pairs` | Standalone stage crashed without coordinates | Reviewer reproduction | Optional columns added | Suite | Fixed |
| S2-12 | L | `validate_inputs.load_fasta` | UTF-8 BOM rejected confusingly | Reviewer reproduction | `utf-8-sig` | Suite | Fixed |
| S2-13 | L | `geocode_locations.fill_missing_coordinates` | GeoNames area/continent results not coarse | Code review (not live) | Class L, CONT/RGN/AREA coarse | Mocked only | Fixed (not live-verified) |
| IMP-1 | H | `processed_snv.load_processed_snvs`, `input_inspection.inspect_table` | iVar lowercase `sample` column → all samples merged | Reproduced | Case-insensitive sample column | Regression | Fixed |
| IMP-2 | H | `processed_snv.load_processed_snvs` | Reported 0% rows counted as detections | Reproduced (became Shared) | Set aside with a warning | Regression | Fixed |
| IMP-3 | M | `input_inspection.compatibility_review` | Promised contig-length check absent | Code review | Mismatch blocks the run | Regression | Fixed |
| IMP-4 / APP-3 | M | `file_formats.plain_copy`, `input_inspection.inspect_file` | Corrupt gzip/XLSX → HTTP 500 | Reproduced | Caught per file | Regression + validation matrix | Fixed |
| IMP-5 | M | `file_formats.vcf_to_snv_table` | Lines without sample fields dropped silently | Reproduced | Counted as set aside | Regression | Fixed |
| IMP-6 | M | `file_formats.allele_measurement` | AD/DP mismatch rejected whole file; denominators differed | Reproduced | Depth = allele-count total | Regression | Fixed |
| IMP-7 | M | `file_formats.vcf_to_snv_table` | Genotype-based INFO/AF accepted as read frequency | Reproduced | INFO/AF only for sites-only VCFs | Regression | Fixed |
| IMP-8 | L | `input_inspection.combine_variant_tables`, importer | Per-file warnings lost; mixed format labeled VCF | Reproduced | Notes kept per file; combined label | Regression | Fixed |
| IMP-9 | L | importer | Frequency vs counts not cross-checked | Reproduced | Warning | Suite | Fixed |
| IMP-10 | L | importer | Conflicting duplicates silent | Reviewer reproduction | Warning names the conflict | Suite | Fixed |
| IMP-11 | L | `column_detection.resolve_choices` | Automatic choice labeled as user's | Code review | `chooser` wording | Suite | Fixed |
| IMP-12 | L | `sequence_variants.build_individual_profiles` | Overhang called insertion; one failure aborted all | Reviewer reproduction | Unaligned overhang; per-sample unmapped | Suite | Fixed |
| SCI-1 | M | `variant_matching.compare_evidence` | QC-excluded records counted as Incomparable reports | Reproduced | Included records only; excluded counted separately | Suite | Fixed |
| SCI-2 | M | `variant_matching.compare_signatures` | QC-excluded call counted as assessable non-detection | Reproduced | Neither present nor assessable | Suite | Fixed |
| SCI-3 | M | `variant_matching.usable_coordinates`, `geographic_summary` | State-level pooled coordinates gave km distances | Reproduced (4.99 km) | Coarse labels excluded | Regression | Fixed |
| SCI-4 | L | `quality_control.apply_pooled_quality_control` | Masked list for another reference → Pass | Reproduced | "Not evaluable" | Regression | Fixed |
| SCI-5 | L | `integrated_analysis.analysis_status` | Partial coverage shown as Completed | Reproduced | "Partly completed" | Suite | Fixed |
| SCI-6 | L | `temporal_analysis.depth_regression` | Tied permutation slopes dropped (p 0.834 vs exact 1) | Reproduced | Relative tolerance | Regression | Fixed |
| APP-1 | H | `app.stage_uploads` | 20 small .gz files expanded to 1.15 GB | Reproduced | One 60 MB budget per request | Validation matrix | Fixed |
| APP-2 | H | `input_inspection.review_dataset` | Non-FASTA as reference → HTTP 500 | Reproduced | Problem message | Regression | Fixed |
| APP-4 | M | `file_formats.numbers` | `DP=inf` → HTTP 500 | Reproduced | Non-finite → missing | Regression | Fixed |
| APP-5 | M | `app.run_reviewed` | Page said uploads deleted; they persisted | Reproduced | Deleted after the run | Validation matrix | Fixed |
| APP-6 | M | `integrated_analysis.build_summary` | Individual-only runs: dates "Not available" | Reproduced | Dates from any source | Validation matrix | Fixed |
| APP-7 | M | `integrated_analysis.run_integrated_analysis` | Variation section empty for individual + reference | Reproduced | Position summary and figure | Validation matrix | Fixed |
| APP-8 | M | `integrated_view.build_dashboard` | Incompatible reference read as "0 shared" | Reproduced | Incomparable count shown | Validation matrix | Fixed |
| APP-9 | M | `integrated_view.build_dashboard` | Links to files too large to cache | Code review | Only cached files linked | Suite | Fixed |
| APP-10 | L | template | Section 8 told users to upload FASTA after a failure | Reproduced | Correct message | Manual | Fixed |
| APP-11 | L | `file_formats.plain_copy` | Messages named internal files | Reproduced | User file names | Regression | Fixed |
| APP-12 | L | `app.clean_staged` | Folders from earlier processes never deleted | Observed 91 folders | Orphan cleanup | Manual | Fixed |
| APP-13 | L | template, view | Sprint 2 counts twice; two "Unresolved" cards | Observed | Removed / distinct labels | Suite | Fixed |
| APP-14 | L | `app.run_reviewed` | Complete-genome confirmation demanded without sequences | Reproduced | Not required without FASTA | Manual | Fixed |
| APP-15 | L | template | "10 minutes" ignored early eviction | Observed | Wording | — | Fixed |
| UX-1 | M | figures / template | 2,100–2,400 px figures unreadable when scaled | Browser check (real table) | Scrolling full-size box | Browser DOM | Fixed |
| PERF-1 | M | `temporal_analysis.observation_matrix` | Real 7,973-record table: ~330 MB peak Python allocations; the observation-matrix CSV exceeds the 8 MB download limit | Measured | Large CSVs are listed as "too large" (not linked) | Measured | Documented (use the CLI for large tables on a 512 MB host) |
| DOC-1 | L | README | VCF INFO/AF statement outdated after IMP-7 | Review | One-line correction | — | Fixed |
| OPEN-1 | M | `geocode_locations.is_state_only_location` | Country-only labels outside the USA are not coarse without GeoNames | Design limit | — | — | Documented |
| OPEN-2 | L | `geocode_locations.fill_missing_coordinates` | No cap on distinct GeoNames lookups per run (8 s timeout each) | Code review | — | — | Open |
| OPEN-3 | L | `quality_control.fisher_check` | iVar PASS is p ≤ 0.05 per iVar docs (not checked); our p < 0.05 differs only at exactly 0.05 | Reviewer note | — | — | Open |
| OPEN-4 | L | `diversity_analysis` | PCoA axis signs arbitrary; diversity is mean per-site Gini–Simpson without depth correction | Reviewer note | — | — | Documented |
| OPEN-5 | L | `app.index` / `check_settings`; label dictionaries | Duplicated settings checks and labels | Code review | — | — | Open (refactor deferred) |
| OPEN-6 | L | `run_integrated_analysis.py` | `--output-folder` equal to the input folder raises `SameFileError` | Reviewer note | — | — | Open |

---

## I. Test results

| Run | Total | Passed | Failed | Skipped |
|---|---|---|---|---|
| Baseline before the audit (real MAFFT on PATH) | 204 | 202 | 0 | 2 (optional real-data tests) |
| After the audit (real MAFFT on PATH) | 224 | 222 | 0 | 2 (optional real-data tests) |
| After the audit, MAFFT not on PATH | 224 | 221 | 0 | 3 (adds the real-MAFFT class) |

- **New tests:** `tests/test_audit_regressions.py` (20 tests). Against the pre-audit checkpoint 19 of them fail (13 failures, 6 errors), demonstrating that they detect the defects; all 20 pass now.
- **Changed tests:** four Sprint 2 tests in `test_app.py` / `test_pipeline.py` used 4-base toy sequences; they now use the same patterns repeated to 20 bases so they satisfy the new overlap rule, with identical p-distances and unchanged assertions.
- **End-to-end validation matrix (real MAFFT, Flask test client, 33 checks, all pass):** Test A Sprint 2 example (8 aligned, 28 comparisons, 2 clusters — the prior baseline); Test B pooled-only (9 samples, exploratory ordination labeled, readable CSVs); Test C integrated example (6 roles detected, all 8 generator answers, 4 of 5 mapped, paper PCoA, 21 figures, provenance, uploads deleted); Test D invalid inputs (empty, invalid characters, corrupt gzip, corrupt XLSX, FASTQ, ambiguous mapping, percentage frequencies, 11 MB → 413, multi-file decompression budget, missing dates/locations); Test E single-source and incompatible references; Test G pages and 404s.
- **Real data:** NCBI H5N1 segment 4 (10 sequences) through the guided workflow with real MAFFT: 45 comparisons, 22 genomic-temporal candidate pairs among 9 samples, 0 strict clusters (state-level locations) — unchanged from the documented baseline. Fontenele Supplementary Table 2: counts in Section E, unchanged after fixes. Optional tests `test_optional_real_ncbi_data.py` and `test_optional_fontenele_table.py` pass when their untracked files are supplied.
- **External dependencies:**
  - *MAFFT:* real executions on Windows (`C:\Tools\mafft-win\mafft.bat`, v7) — the real-MAFFT test class, the examples, the H5N1 data, and an input folder containing spaces. Linux MAFFT (Docker) not executed here.
  - *GeoNames:* **not executed live** — `GEONAMES_USERNAME` is not set in this environment. All GeoNames behaviour is covered by mocked tests (success, failure, coarse results, caching). Unverified live: the actual service response format today and the new class-L/area coarseness rule.
  - *Docker:* build **not executed** (Docker daemon not running). The Dockerfile copies every module the app imports.
- **Performance (Windows, wall clock):** Sprint 2 example 2.0 s; integrated example 2.1 s (Sprint 2 stage) + 14.2 s (reference alignment of 5 sequences through the MAFFT batch file, analyses, 21 figures); pooled-only synthetic 2.4 s; real 7,973-record table 14.3 s, ~330 MB peak Python allocations, 37 MB of CSV output (the largest CSVs are not cached for download). All inputs used are under the 10 MB upload limit (`mmc2.xlsx` is 0.95 MB). A genome-wide A/C/G/T table for many samples would exceed 10 MB and belongs in the CLI.

---

## J. Local execution instructions (Windows)

Command Prompt, from the project folder (replace the path):

```bat
cd /d "C:\path\to\sprint4_claude"
py -3.11 -m venv .venv
.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\python.exe -m pip install -r requirements.txt
```

Verify MAFFT (Windows package with `mafft.bat`, here installed in `C:\Tools\mafft-win`):

```bat
set "PATH=C:\Tools\mafft-win;%PATH%"
where mafft
```

Configure GeoNames (optional; use your own free GeoNames username, never commit it):

```bat
set "GEONAMES_USERNAME=your_geonames_username"
```

Run the test suite:

```bat
.venv\Scripts\python.exe -m unittest discover -s tests -v
```

Start Flask, then open http://127.0.0.1:5000 (use **Check MAFFT & GeoNames** for the live diagnostics):

```bat
.venv\Scripts\python.exe app.py
```

Integrated synthetic example: click **Run Integrated Example** on the first page, or from the command line:

```bat
.venv\Scripts\python.exe run_integrated_analysis.py --example --output-folder results\integrated_example
```

Command-line analysis of your own files (gzip inputs are accepted):

```bat
.venv\Scripts\python.exe run_integrated_analysis.py --fasta sequences.fasta --metadata metadata.csv --reference reference.fasta --snv pooled.vcf.gz --output-folder results\my_run
```

Outputs: every table (CSV) and figure (PNG) is written to the `--output-folder`; the Sprint 2 pairwise, related-pair, cluster, and region-assignment CSVs are written next to the copied metadata in that folder. In the web app, use the **Downloads** links in each section or the list under **Provenance and downloads** (kept in memory for up to 10 minutes).

Optional real-data checks (files are not in Git):

```bat
set FONTENELE_MMC2=C:\path\to\mmc2.xlsx
.venv\Scripts\python.exe validation\check_fontenele_2021.py %FONTENELE_MMC2%
```

---

## K. Future manual validation checklist (after the Sprint 4 integration)

Uploads (one at a time, then mixed in one upload):
- [ ] FASTA, FA, FNA, and `.fasta.gz` individual sequences; a one-sequence FASTA as reference (role suggested only with evidence)
- [ ] Metadata as CSV, TSV, XLSX, `.csv.gz`; NCBI Virus export; headers with odd capitalization and spacing; a file with both `Date` and `Collection_Date` (must ask)
- [ ] Standardized SNV table, combined iVar TSV (uppercase and lowercase sample column), per-sample iVar files, the Fontenele `mmc2.xlsx`, multi-sample VCF with AD, sites-only LoFreq VCF, `.vcf.gz`, GT-only VCF (must be set aside)
- [ ] Pooled metadata, per-position depth, per-position A/C/G/T counts, masked sites, documented profiles
- [ ] Refusals: FASTQ, BAM, empty file, corrupt gzip, corrupt XLSX, over 10 MB, many large `.gz` files

Analysis modes:
- [ ] Individual only (with and without metadata; stated pathogen; complete-genome confirmation)
- [ ] Individual + reference; pooled only; both sources; both with a different reference label; both with a same-label but different-length reference (must block)
- [ ] Paper versus exploratory QC; changed thresholds and date window

Dashboard (for each mode, all eight sections):
- [ ] Purpose, statistics (no misleading 0), figures with captions, limitations, downloads
- [ ] Cross-source categories against a known example; overlap matrix colours; incomparable explanation
- [ ] Temporal sampling timeline and trajectories; PCoA method label (paper vs exploratory)
- [ ] Geography precision classes; coarse labels never plotted as points
- [ ] Sprint 2 clusters vs candidates; network wording (never transmission)
- [ ] Wide figures scroll; phone-width layout has no horizontal page scroll

Sample results to compare: integrated example (2 shared, 1 individual-only, 1 pooled-only, 1 + 1 unresolved, 2 incomparable); Sprint 2 example (28 comparisons, 2 clusters); real H5N1 (45 comparisons, 22 candidates, 0 clusters); real `mmc2.xlsx` (Section E counts).

Charts and downloads: every PNG opens; every listed CSV downloads and opens; provenance lists inputs, versions, reference checksum, parameters, and exclusions.

Error handling: each refusal above shows a readable message, never a server error page; expired upload links ask for a new upload.

Geographic enrichment (including any Sprint 4 literature-derived locations): provenance of each location visible; coarse or inferred locations never used for kilometre screening; GeoNames live lookup with a real username (Check MAFFT & GeoNames shows Operational).

External tools: MAFFT on Windows and in Docker; GeoNames live; Docker build and container start; Render deployment only after approval.

---

## L. Readiness for Sprint 4 (NLP / literature-based metadata enrichment)

**Assessment: ready to build on, with clear integration points.** No Sprint 4 code was added.

Existing interfaces a future evidence-based enrichment component can use:

1. **Metadata input boundary.** All individual metadata passes through `validate_inputs.validate_metadata` (and, for the guided upload, `input_inspection.prepared_metadata`, which already adds values from explicit FASTA annotations with a `<Field>_Source` column and never overwrites existing values). A literature-derived location or date could enter the same way: as a new value with its own source column, applied only where the field is missing.
2. **Geographic precision model.** `Coordinate_Source` and `Geographic_Precision` (set in `geocode_locations.fill_missing_coordinates`) already distinguish provided, GeoNames, coarse, and missing locations; `integrated_analysis.precision_category` maps them to five display classes. An enrichment step needs its own precision class (for example "literature-derived, unverified") that `usable_coordinates` and the Sprint 2 kilometre rule treat as **not** precise unless it is verified.
3. **Provenance.** `integrated_analysis.provenance` records data sources, versions, parameters, and exclusions per run; `analysis_status` records what ran. Enrichment should add rows (source document, extraction method, confidence) rather than a separate report.
4. **Review step.** The Step 2 review page already shows "what was recognized, inferred, and missing". Enriched values should appear there for user confirmation before a run, consistent with the no-silent-guessing rule.

Compatibility requirements for Sprint 4:

- Never overwrite a value supplied in the uploaded metadata; record the source of every added value.
- Keep inferred or literature-derived locations out of strict kilometre screening unless verified at site level; keep GeoNames HTTPS, timeouts, and caching, and add a cap on lookups per run (OPEN-2).
- Keep sample identifiers as text end to end (see S2-1) and match by the same normalized accession rules.
- Keep the new component off the request path or bounded in time (the hosted worker stops after 180 s).
- Add tests in the same style (known answers, provenance assertions), and keep mocked versus live external calls clearly separated.
