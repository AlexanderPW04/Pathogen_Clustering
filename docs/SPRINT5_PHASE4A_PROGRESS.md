# Sprint 5 Phase 4A — Progress (generalized import and multi-dataset validation)

The goal loop resumes from the last milestone marked **Done**. Nothing is committed, pushed,
merged, or deployed; `AI_ENRICHMENT_ENABLED=false` in every run; published files stay unchanged.

Test command used throughout (project `.venv`, MAFFT on PATH, real workbook):

```bash
FONTENELE_MMC2=<scratchpad>/data/mmc2.xlsx AI_ENRICHMENT_ENABLED=false PATH=/c/Tools/mafft-win:$PATH .venv/Scripts/python.exe -m unittest discover -s tests
```

## Milestone 0 — Current application (attempt 1) — Done

**Verified baseline (this session, before any change):**

- Branch `sprint4-integration`, HEAD `b641533`, 67 uncommitted paths (Sprint 5 Phases 1–4).
- `Ran 636 tests in 448.2 s — OK (skipped=5)`.
  - The goal text quotes 585 tests and 4 skips. That was the Phase 3 baseline. Phase 4 (mutation
    annotation) added 51 tests and one optional skip; its report gives 636 / 5. The 5 skips are
    3 SciPy cross-checks, the NCBI-files test, and the mpox test.
- `mmc2.xlsx` md5 `c0a75b7fa5bde8db09746613354a9f64` (unchanged).
- Reference outputs from the current code were saved before any change, for the stability
  comparison in Milestone 9:
  - `p4a_base_mmc2`: `mmc2.xlsx` alone;
  - `p4a_base_mmc2_ann`: with the NCBI GenBank annotation;
  - `p4a_base_example`: the synthetic paired example.

**Reports read.** Sprint 3 audit, Phase 1–4 reports (Phase 3 in full, the others by section).
Phase 3 conclusions that bind this phase:
- depth-aware grouping failed criteria P1/P2 → opt-in only, labelled "not validated";
- replicate pairs need a shared `Specimen_ID` and distinct documented `Library_ID`s;
- unknown is never 0.

**Architecture (verified by reading code).**

| Stage | Where | What it does today |
|---|---|---|
| Upload | `app.review_upload` → `stage_uploads` | 10 MB request limit, 30-min staging, ≤ 20 files |
| File kind | `file_formats.detect_kind` | content-based: FASTA, VCF, GenBank/GFF3, FASTQ/BAM refused, CSV/TSV/XLSX |
| Worksheet and header row | `excel_tables.describe_sheets`, `read_sheet`, `find_header_row` | header = row with most recognized names in the first 15; several data sheets → user chooses |
| Format | `input_inspection.inspect_table` | Fontenele Table 2 (title row + exact headers), iVar, VCF, or "table with recognizable columns" |
| Column mapping | `column_detection.detect_columns`, `resolve_choices`, `value_warnings` | alias sets per role; 2+ candidates → user must choose; missing required field → problem |
| Roles | `suggest_roles`, `review_dataset` | FASTA, metadata, reference, annotation, SNV, pooled metadata, positions, masked sites, signatures |
| Review | `individual_review`, `pooled_review`, `compatibility_review`, `analysis_plan`, `metadata_readiness` | IDs, dates, references, size limit, per-analysis plan (3 states) |
| Standardized copies | `prepare_inputs`, `prepare_variant_file`, `combine_variant_tables` | mapped columns renamed; multiple SNV files combined |
| Load and validate | `processed_snv.load_processed_snvs` | bases, positions, frequency in [0,1] (no rescaling), depth, Alt/Ref counts, metadata join (raises on conflicts or unknown samples), duplicates flagged |
| Shared schema | `unified_schema`, `pooled_observations` | variant rows, sample rows, per-position depth / A/C/G/T lookup |
| Analyses | `integrated_analysis.run_integrated_analysis` | QC, Sprint 3 summaries, observation matrix, trends, PCoA (Yue–Clayton or Jaccard), Sprint 5 VAF/Wilson, beta-binomial, overdispersion, profiles/grouping, annotation, cross-source comparison |
| Results | `integrated_view.build_dashboard`, `download_files`, `app.store_result_files` | CSV + PNG only, zlib in memory 10 min, files > 8 MB skipped and listed |

**Gaps found (what Phase 4A must add; nothing below exists yet):**
1. No frequency-scale choice: percentages are refused, never rescaled.
2. No 0-based position option; no stated reference when a table has no reference column.
3. `Read_Depth` is a required column even for frequency-only tables.
4. No wide (sample × mutation) matrix reader; no base-count (A/C/G/T + Ref) variant derivation.
5. No header-row override; no reusable mapping profile.
6. Metadata join: conflicts, duplicate IDs, and unmatched samples stop the run with an error;
   there is no matching report and no way to keep a conflict explicitly unresolved.
7. Analysis plan has 3 states and no "Experimental" state; it is not per Sprint 3 + 5 method.
8. No XLSX or PDF export; no documented canonical tables or import-provenance download.

**Published datasets obtained (originals kept unchanged in `<scratchpad>/datasets/`).**
- **Crits-Christoph et al. 2021** (mBio, doi 10.1128/mBio.02703-20).
  - Supplementary Tables S1–S4, via the Europe PMC supplementary-files API. PMC itself served
    a proof-of-work bot check, which was not bypassed.
  - Also downloaded from GitHub `alexcritschristoph/wastewater_sarscov2`: `reference.fasta`,
    `tables/`, `final_tables/`, the inStrain SNV files, and the metadata. GPL-3 licence.
- **Schumann et al. 2022** (Sci Total Environ, doi 10.1016/j.scitotenv.2022.158931).
  - mmc1–mmc4 via the Europe PMC API.
  - PiGx SARS-CoV-2 source at the paper's commit `524ed48`, used to define the table columns:
    the reference is `NC_045512.2`, frequencies are LoFreq AF, and NA means the mutation was not
    called in that sample.

## Milestone 1 — Canonical sequencing-data model (attempt 1) — Done

- **New `genomics/canonical_model.py`.** Builds three linked tables plus an import record from
  values the analyses already computed. It selects and renames only; it calculates nothing.
  - `canonical_samples.csv`: one row per sequenced sample (library). Columns are Specimen,
    Library, date (parsed and original text), site, coordinates, organism, references (one per
    segment), platform, coverage, QC, source dataset, and the source of the metadata.
  - `canonical_variant_observations.csv`: one row per reported observation.
    - `Reported_VAF` and `Calculated_VAF` are separate columns.
    - When the importer calculated a frequency from base counts, `Reported_VAF` is left blank
      and `Frequency_Provenance` says so.
    - Also carries count, depth, and frequency provenance; `Measurement_Class` (counts + depth,
      frequency + depth, or frequency only); QC; reference status; and every
      `Source_*_Original` value an import option changed.
  - `canonical_coverage.csv`: per-position A/C/G/T, callable depth, and coverage status from a
    per-position table. It is empty without one: never zero-filled, and unreported stays
    "coverage unknown".
  - `import_provenance.csv`: per file, one row for each layout, worksheet, header row, field
    column, and transformation, and who decided it (detected, user, documented format).
- **Wiring.**
  - `run_integrated_analysis(import_details=...)` adds these tables plus
    `analysis_readiness.csv` and `metadata_matching.csv`.
  - `pooled["metadata"]` is kept so the sample table can show the platform.
- **Column detection.** New `Sequencing_Platform` field. A column named only "sequencing" is
  deliberately not an alias: in Crits-Christoph Table S1 it means enriched or unenriched.
- **Check.** The synthetic example gives 9 samples, 22 observations, and 960 coverage records.

## Milestone 2 — Flexible detection and mapping (attempt 1) — Done

- **New `genomics/table_layouts.py`.**
  - **Wide matrices.** A header is a mutation only if it states one single-nucleotide change:
    `C10029T`, `10029C>T`, `g.10029C>T`, or PiGx `GENE:AA::C10029T`.
    - Non-SNV change headers are listed and not analyzed.
    - Blank or "NA" cells become no record (never 0%); sample rows with no records are listed.
  - **Base-count tables** (Ref + A/C/G/T, no Alt). Each non-reference base with at least the
    user-confirmed minimum number of reads becomes one record.
    - Counts are measured; depth = A+C+G+T.
    - The frequency is labelled as calculated by the importer.
    - Columns that could be mistaken for derived values are renamed `Source_*`. A source depth
      column is compared with A+C+G+T.
    - inStrain's `var_base` and `var_freq` are never interpreted: `var_base` can be the
      reference base (Table S4 row 1738).
- **New `genomics/import_options.py`.** Per-file choices on the review page:
  - frequency scale (percent ÷ 100 with exact `Decimal` arithmetic; original kept);
  - 0-based or 1-based positions (+1; original kept);
  - stated reference (recorded as user-stated);
  - depth definition (provenance only);
  - base-count minimum;
  - wide-layout confirmation;
  - header row;
  - sample ID from the file name.
- **Numbering check.** `reference_check` / `numbering_warning` warn, but never shift, when Ref
  bases fit the other numbering. inStrain files: 16–30% match as 1-based, 100% as 0-based.
- **`column_detection.py`.**
  - `Reference_Genome` is no longer required: a stated reference replaces it, and the review
    asks for one.
  - `Read_Depth` is now optional: a frequency-only table is accepted with depth left blank.
  - Added `Ref_Depth`, coordinate fields in the SNV role, and aliases (`genome_position`,
    `coverage`, `position_coverage`, `scaffold`, `dates`, `coordinates_lat/long`, alternative
    read counts).
  - "site" is still not a position alias.
- **Header row.** `excel_tables.read_sheet(header_row=...)` and `file_formats.read_table(header_row=...)`.
- **`input_inspection.py`.** Layouts at inspection; `table_for(path, info, role, options)`;
  `interpretation_review`; `numbering_review`; `interpreted_table`; `import_provenance`;
  `Source_File` per record when files are combined.
- **Not done here.** Saved mapping profiles (JSON) are done in Milestone 6 with the review form.

## Milestone 3 — Multiple files and metadata matching (attempt 1) — Done

- **New `genomics/metadata_matching.py`.** Joins only identical `Sample_ID` text. It reports:
  - matched, unmatched, and unused samples;
  - rows without an ID (for example Table S1's "average:" row);
  - identical and conflicting duplicate IDs;
  - conflicts between the variant table and the metadata, and within one sample;
  - possible matches that differ only in case, punctuation, or leading zeros, listed and
    **not** joined;
  - specimens with several samples, and whether each has distinct documented `Library_ID`s;
  - segments;
  - unusable dates.
- **Conflicts.** The user picks one choice for all conflicts: "metadata", "variant table", or
  "unresolved" (both blanked and recorded). Duplicate-row and within-sample conflicts accept
  only "unresolved". No choice → blocking problem.
- **Unmatched samples** stay a blocking problem unless the user ticks "Continue without
  metadata"; their values stay unknown.
- **Coordinates in a variant table** are moved into the sample metadata.
- **Fix.** The review counted every metadata row with a date, including rows of samples with no
  variant records: Crits-Christoph showed "27 of 7". It now counts only matched samples.

## Milestone 4 — Readiness assessment (attempt 1) — Done

- **New `genomics/analysis_readiness.py`.** Gives one of four statuses (Available / Partially
  available / Unavailable / Experimental), with a reason and needs, for:
  - the 14 Sprint 3 + 5 methods in the goal;
  - the Phase 4 annotation;
  - the Sprint 2 rows when FASTA files are present.
- **Basis.** Statuses come from mapped and checked inputs: depth fill, counts (measured,
  derived, reconstructable), reported versus calculated frequency, dates, A/C/G/T tables,
  replicate candidates, reference FASTA, annotation match, and the grouping setting. A paper or
  file name is never used.
- **Example reasons** (goal wording): "Depth is available, but alternative read counts were not
  directly provided"; "Technical replicates were not verified…"; "Full nucleotide coverage is
  unavailable…"; "Mutation grouping is experimental and disabled by default".

## Milestone 5 — Published datasets (attempt 1) — Done

- **New `validation/check_published_datasets.py`.** Uses the website's import functions with
  the user's choices, then independent openpyxl/pandas checks. md5 checksums are taken before
  and after.
- **Results** (`Ran … FAILED checks: 0`):

| Dataset | Structure | Checks |
|---|---|---|
| Fontenele 2021 `mmc2.xlsx` | Published long table, title row | 7/7: 7,973 records, 52 samples, frequencies identical, depth = TOTAL_DP, 7,973/7,973 Ref bases match NCBI NC_045512.2, counts labelled reconstructed, md5 unchanged |
| Crits-Christoph 2021 Table S4 + S1 | Long base-count table (inStrain), no reference column, no Alt column; separate metadata | 13/13: join refused (0 of 7 IDs match S1), "average:" row reported, 294 records = independent derivation at ≥ 2 reads (paper's minimum), calculated VAF = the study's con/var frequency on 293/293 rows, 1-based (293/293 Ref match WIV04), WIV04 ≠ NC_045512.2 (29,840 vs 29,903 bases) |
| Crits-Christoph GitHub | 8 per-sample inStrain TSVs (no sample column, 0-based) and `final_tables/Wastewater_SNVs_v3.csv` (0-based long) | 0-based warning on 8/8 files; after +1, v3 reproduces S4 exactly (294/294); the per-sample files are an earlier run (250 of 263 sites shared; one-base error would give 25) |
| Schumann 2022 mmc1 S2.1 | Wide 38 × 1,907 matrix, title row, NA cells, no depth, no reference column | 12/12: 5,145 records = non-NA SNV cells, 5,145/5,145 frequencies exact, 67,093 NA → no record, depth never invented, per-sample counts = authors' hidden `total_muts` 38/38, 2 samples with no SNV listed, dates 36/36, Ref bases match NC_045512.2 (PiGx default), **Phase 4 amino-acid predictions = authors' VEP labels 1,830/1,830**, mmc3 (14 MB, 4.7 M cells) refused with an actionable message |

- **Real-data findings.**
  - Crits-Christoph S4 and S1 sample names do not correspond in any documented way: `5_28_S`
    against `5_28_I`, `5_19_S` / `5_19_S_2` against `5_19_F_S3` / `5_19_A_S4`. No crosswalk was
    invented. Dates are therefore unavailable, and temporal analyses are honestly Unavailable.
  - Schumann's NA is PiGx "not called" (`full_join` of per-sample tables), so it means not
    reported, never 0.

## Milestone 6 — Understandable upload workflow (attempt 1) — Done

- Five step labels: Upload → Review detected data → Review available analyses → Run → Results.
- **Step 2.**
  - "How to read each table" (`templates/_review_import.html`) shows only the options that apply
    to each table.
  - "Sample metadata matching" holds the report, the conflicts, the resolution choice, and the
    "continue without metadata" box.
  - Mapping profile save/load: the `/review/profile` route and the `use_profile` helper.
  - Technical tables are collapsed.
- **Step 3** (`templates/_readiness_step.html`): status counts and one row per analysis.
- **Tested headings kept:** "Step 2 · Review your dataset" and "Analyses available with these
  files". The run step is now "Step 4 · Run analysis".
- **Fix during testing.** The new "use the file name as the sample ID" hint also fired on tables
  missing frequencies. It now requires a table that is complete except for the sample column.
  `test_unsupported_workbook_says_which_values_are_missing` passes again.
- **Fix during testing.** Profiles now also record FASTA and annotation roles (matched by kind; a
  reference needs exactly one sequence).

## Milestone 7 — Results and downloads (attempt 1) — Done

- **New `genomics/report_exports.py`.**
  - `analysis_results.xlsx`: a Summary sheet plus up to 19 key tables. It uses openpyxl write-only
    mode, and the Excel row limit is reported, not silent.
  - `analysis_report.pdf`: dataset, the dashboard's findings, Step 3 readiness, methods, import
    choices, limitations, and up to 8 figures (wide figures on landscape pages). Built with
    matplotlib PdfPages; no new dependency.
- **`integrated_view.download_files`** adds both. A failure becomes a page message and never
  fails the run.
- **The dashboard** gets a download bar and collapsible readiness, provenance, and matching
  sections.
- **`run_integrated_analysis.py --reports`** writes both files locally.
- **Check.** For the Schumann run, the XLSX sheets Observations (5,145 × 23), Samples, and
  Readiness equal the CSV downloads exactly. The PDF was rasterized with PyMuPDF and inspected:
  10 pages.

## Milestone 8 — MAGOS-inspired foundation (attempt 1) — Done

- **Grouping code, simulations, and benchmarks unchanged.** `DEFAULT_SETTINGS["grouping_mode"]`
  is still `"profiles"`. Readiness always labels grouping "Experimental", with the depth-aware or
  ordinary variant named by the counts available.
- **Inputs it can now use.** Base-count layouts feed measured counts to the depth-aware variant;
  frequency-only matrices feed only ordinary clustering. Grouping is never promoted.
- **Future-work assessment** (annotation, consequences, maps, real-data group validation, FASTQ
  outputs, other pathogens): report §15–16.

## Milestone 9 — Tests, performance, integration (attempt 1) — Done

- **New tests:** `tests/test_phase4a_import.py` (29), `tests/test_phase4a_web.py` (6), and
  `tests/test_phase4a_published.py` (3 optional; they ran with the real files).
- **Final suite:** `Ran 674 tests in 651.6 s — OK (skipped=5)`.
- **Fixes made during testing:**
  - The sample-ID hint is now limited to tables that are otherwise complete.
  - Two Sprint 4 tests checked the old step label. The label now reads "4 · Run · Metadata
    enrichment (Sprint 4, optional)" and the assertions were updated.
  - The pooled metadata readiness listed depth as Required. It is now Recommended, and a stated
    reference counts as present. One test expectation was updated.
  - The note "no reference FASTA uploaded" appeared even when one was uploaded. It now checks
    whether any file has the reference role.
- **Stability.** All 42 / 48 / 64 earlier `mmc2`, `mmc2` + GenBank, and example output files are
  identical. The only additions are 4 canonical tables, plus XLSX/PDF with `--reports`.
- **Performance regression, found and fixed.**
  - Building the XLSX and PDF inside `/run` made the 10%-CPU run take 198 s, above the 180 s
    gunicorn timeout.
  - They are now built on demand from the stored files. Large detail tables stay CSV-only, and PDF
    figures are 8-bit and at most 1,200 px.
  - Final measurements: `mmc2` 308 MB, 12.5 s; Render-like 300 MB, 126 s. Details in report §13.
- **Visual check.** A headless Edge screenshot of the saved review page found the two readiness
  contradictions above, and both are fixed. The in-app browser cannot open localhost.

## Milestone 10 — Audit preparation and final report (attempt 1) — Done

- `docs/SPRINT5_PHASE4A_AUDIT_CHECKLIST.md`: 14 areas. Every check is tagged with the kind of
  evidence it needs: code execution, mathematics, scientific validity, real data, or hosting.
- `docs/SPRINT5_PHASE4A_IMPLEMENTATION_REPORT.md`: 19 sections plus the completion table.
- **Attempts used:** 1 per milestone (18 allowed). Fixes found during testing are listed above.

**SPRINT 5 PHASE 4A COMPLETE**
