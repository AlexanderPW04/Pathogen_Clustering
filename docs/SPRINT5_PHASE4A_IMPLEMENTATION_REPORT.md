# Sprint 5 Phase 4A — Implementation Report

**Generalized sequencing-data import and multi-dataset validation**

- Branch `sprint4-integration`, HEAD `b641533`. All Sprint 5 work (Phases 1–4A) is local and
  uncommitted.
- Nothing was committed, pushed, merged, or deployed.
- No paid API call was made: `AI_ENRICHMENT_ENABLED=false` in every run.
- Published files were only read. Their md5 checksums were identical before and after every
  validation run.

Labels used below:
- **verified**: run or read in this session;
- **synthetic**: data invented by the tests;
- **real**: a published file.

The step-by-step record is in `docs/SPRINT5_PHASE4A_PROGRESS.md`, and the later audit plan is in
`docs/SPRINT5_PHASE4A_AUDIT_CHECKLIST.md`.

## 1. Executive summary

The Sprint 3 + 5 pipeline now reads processed pathogen-sequencing results from other laboratories.
Their tables, worksheets, column names, layouts, metadata, and available measurements can differ.
None of the existing scientific calculations was changed. What was added:

- **Two new readers for common published layouts.**
  - *Wide sample-by-mutation matrices*: PiGx-style tables (Schumann et al. 2022).
  - *Base-count tables*: Sample, position, Ref, and A/C/G/T counts with no Alt column
    (inStrain-style tables, Crits-Christoph et al. 2021).
- **Choices only the source documents can settle.** Where headers cannot settle a question, the
  user confirms the answer from the source's documentation:
  - frequency scale (fraction or percent);
  - 0-based or 1-based positions, checked base by base against a reference FASTA;
  - the reference when a table has no reference column;
  - the depth definition;
  - the minimum reads for a base-count variant;
  - the wide layout;
  - the header row;
  - the file name as the sample ID.
- **Provenance.** Every choice is recorded, and every changed value keeps its original.
- **Mapping profiles.** A confirmed mapping can be saved and reused as a small JSON profile (no
  records). It is applied only to files with the same layout and every mapped column.
- **A metadata-matching report.** It joins only on identical IDs; near matches are listed, never
  joined. It reports duplicates, rows without IDs, and conflicts. Conflicts are settled only by an
  explicit user choice, or kept explicitly unresolved.
- **A readiness assessment** with Available / Partially available / Unavailable / Experimental
  states. It covers 14 Sprint 3 + 5 methods, the Phase 4 annotation, and the Sprint 2 rows. It is
  decided from the mapped, checked values, never from a file or paper name.
- **Canonical tables.**
  - Sample metadata.
  - Variant observations, with separate reported and calculated VAF and count/frequency
    provenance.
  - Coverage records.
  - Import provenance.
- **Two new downloads:** one consolidated XLSX workbook and one concise PDF report. They are built on
  demand from the stored CSV/PNG downloads, so the workbook equals the CSVs value for value and the
  analysis request is not slowed down.
- **A five-step workflow:** upload → review detected data → review available analyses → run →
  results. Technical details sit in collapsed sections.

**Real-data validation (verified).** Three published wastewater studies, with five structurally
different files, were imported end to end. Each was compared with values computed independently
from the original files:

| Dataset | Structure | Independent checks |
|---|---|---|
| Fontenele et al. 2021 | Long table with title row, published annotations | 7/7 pass |
| Crits-Christoph et al. 2021 | Long base-count table (S4) with separate metadata (S1); 8 per-sample inStrain TSVs (0-based, no sample column); 0-based long `final_tables` CSV | 13/13 pass |
| Schumann et al. 2022 | Wide 38 × 1,907 matrix with title row and NA cells; 109 × 8,884 matrix; 384 × 12,359 matrix (too large, refused with an explanation) | 12/12 pass |

**Regression and safeguards.**
- Every one of the 42 / 48 / 64 earlier output files is identical to the saved baseline: `mmc2`,
  `mmc2` with annotation, and the synthetic example, including all figures. The only additions
  are four canonical tables.
- Experimental mutation grouping is still off by default.

**Tests (verified).** `Ran 674 tests — OK (skipped=5)`: the 636 pre-existing tests plus 38 new.
- The 5 skips are the same optional checks as at baseline.
- All real-data tests ran: Fontenele, Crits-Christoph, Schumann, and the real-workbook web export.

**Performance (verified).** One `mmc2` run takes 308 MB and 12.5 s, against 296 MB and 12.0 s
before. Under 512 MB and 10% CPU it takes 300 MB and 126 s, against 288 MB and 125 s. The XLSX and
PDF are built only when downloaded.

## 2. Input formats supported

| Input | Before Phase 4A | Now |
|---|---|---|
| Fontenele Table 2, iVar TSV, VCF (AD/AO/AF), standard long CSV/TSV/XLSX | Supported | Unchanged (format paths untouched) |
| Long table **without a depth column** | Refused (Read_Depth required) | Accepted; depth stays blank (unknown); depth-based analyses are Unavailable |
| Long table **without a reference column** | Refused | Accepted after the user states the reference; it is checked against a reference FASTA, or labelled unverified |
| **Percentages** (0–100) | Refused, never rescaled | Converted only after the user chooses "Percentages"; exact decimal division; original kept in `Source_Alt_Frequency_Original` |
| **0-based positions** | Not supported; Ref mismatches failed QC | Converted only after the user's choice (+1; original kept); a warning appears when Ref bases fit the other numbering |
| **Wide sample-by-mutation matrix** | Not recognized | Recognized when 3 or more headers each state one SNV; the user must confirm; NA → no record |
| **Base counts with Ref, no Alt** | Recognized as a positions table only (and only with a reference column) | Also usable as variant records with measured counts after the user enters the minimum reads; a source depth column is compared with A+C+G+T |
| **Per-sample files without a sample column** | Refused (except iVar) | "Use the file name as the sample ID" option |
| **Header row below notes** | Detected for Excel only | Detected, or set by the user for CSV, TSV, and Excel |
| **Coordinates in the variant table** | Dropped | Moved into the sample metadata |
| Raw FASTQ / BAM / CRAM | Refused | Unchanged: refused (offline workflow only) |

## 3. Canonical data model

The import still writes standardized copies that the validated loader reads, so no calculation
changed. `genomics/canonical_model.py` then exposes the result as three linked tables, plus a
record of the import. It selects and renames values only.

**A. `canonical_samples.csv`.** One row per sequenced sample (library).
- Identity: Sample_ID, Specimen_ID, Library_ID, Sample_Type.
- Date and place: Collection_Date (parsed) and Collection_Date_Original, Collection_Site,
  Latitude, Longitude, Geographic_Precision.
- Genome: Organism, Reference_IDs, Genome_Segments (one reference per segment).
- Sequencing: Sequencing_Platform, Genome_Coverage_Percent, Mean_Depth, Coverage_Information.
- Counts and QC: Reported_Observations, Sample_QC.
- Provenance: Source_Dataset, Metadata_Source.

**B. `canonical_variant_observations.csv`.** One row per reported observation.
- Position: Reference_ID, Genome_Position (1-based), Ref, Alt.
- Counts and depth: Ref/Alt read counts, Total_Depth.
- Frequencies:
  - **Reported_VAF**: blank when the source reported none;
  - **Calculated_VAF**: reads / depth;
  - VAF_Used.
- Provenance: Count_Provenance (measured, reconstructed (approximate), or frequency only), plus
  depth and frequency provenance.
- Measurement_Class.
- QC_Status, Duplicate_Status, Reference_Status.
- Every `Source_*_Original` value that an import option changed.

**C. `canonical_coverage.csv`.** Per-position A/C/G/T, total callable depth, and coverage status
from a per-position table. It is empty without one. It is never zero-filled, and unreported
variants stay "coverage unknown".

**`import_provenance.csv`.** One row for each file, layout, worksheet, header row, field column,
and transformation, with who decided it: detected, user, documented format, or importer with a
user-confirmed minimum.

Linking keys:
- samples ↔ observations: `Sample_ID`;
- observations ↔ coverage: (`Sample_ID`, `Reference_ID`, `Genome_Position`).

The original uploads are never changed. Transformed copies exist only in the temporary run
folder.

## 4. Column detection and mapping rules

- **Alias matching** (`genomics/column_detection.py`, unchanged algorithm):
  - case, spaces, and punctuation are ignored;
  - exactly one candidate → used, with the evidence recorded;
  - two or more candidates → the user must choose;
  - a missing required field → a stated problem.
- **New aliases.**
  - `genome_position`, `genomic_position`, `nucleotide_position`.
  - `coverage`, `position_coverage`, `base_coverage` as per-record depth (genome-coverage
    percentages keep their own field).
  - `scaffold` as the reference.
  - `dates`.
  - `coordinates_lat` / `coordinates_long`.
  - Alternative and reference read-count names.
  - `Sequencing_Platform`.
- **Similar names that are deliberately not equated.**
  - "Site" is never a position.
  - A column named only "sequencing" is not a platform: in Crits-Christoph S1 it means
    enriched/unenriched.
  - inStrain `var_base` and `var_freq` are never Alt or frequency: in Table S4 row 1738,
    `var_base` is the reference base.
- **Values are checked, not only headers.**
  - Frequencies above 1 → "may be a percentage" (blocking until the user chooses a scale).
  - Percent values outside 0–100 → problem.
  - Non-numeric positions, frequencies, or depths → problem.
  - Fractional or negative base counts → refused.
  - 0- versus 1-based fit against the reference FASTA (warning above 95% / below 50%).
- **Profiles.** Profiles store `{role, kind, format, layout, sheet, fields, options}`. A profile
  entry fits a file only with the same layout and format and every mapped column. Incompatible
  files are listed with the missing columns.

## 5. Metadata matching (`genomics/metadata_matching.py`)

- **Joins use identical `Sample_ID` text only** (surrounding spaces removed). IDs that differ only
  in case, punctuation, or leading zeros are listed as "Possible matches that were NOT joined".
- **The report covers:**
  - matched, unmatched, and unused samples;
  - rows without an ID;
  - duplicate IDs (identical rows kept once; differing rows are conflicts);
  - variant-table versus metadata conflicts and within-sample conflicts;
  - specimens with several samples, and whether each has distinct documented Library_IDs;
  - reference segments;
  - unusable dates.
- **Conflicts.** One choice is applied to all listed conflicts:
  - "Use the metadata file's value";
  - "Use the variant table's value";
  - "Leave unresolved": both are blanked and the conflict is recorded in
    `Source_Metadata_Resolution` and `metadata_matching.csv`.

  Duplicate-row and within-sample conflicts accept only "unresolved". With no choice, the run is
  blocked.
- **Unmatched variant samples** block the run until the user ticks "Continue without metadata for
  unmatched samples". Their dates and sites then stay unknown.
- **Fix to an existing count.** The review counted every metadata row with a date, including rows
  of samples with no variant records (Crits-Christoph showed "27 of 7 samples"). Only matched
  samples are counted now.
- **Equivalence (test).** Separate variant and metadata files give results identical to a combined
  table: `sprint5_vaf_observations.csv`, `temporal_trends.csv`, and
  `pooled_samples_with_qc.csv` are equal apart from the data-source label.

## 6. Data-readiness criteria (`genomics/analysis_readiness.py`)

| Analysis | Available when | Partially available when | Unavailable when |
|---|---|---|---|
| Basic SNV summaries | Any pooled variant file | — | No pooled file |
| Sequencing-depth QC | Every record has depth | Some records have depth (others "not evaluable", never passed); VAF depth checked per record | No depth |
| Reported VAF summaries | Frequencies reported | Only importer-calculated frequencies (base counts) | — |
| Count-based VAF | Alternative counts measured (supplied, or A/C/G/T) | Depth + reported frequency (counts reconstructed where unique, labelled approximate); counts for some records | No counts and no depth |
| Wilson intervals | Measured counts | Reconstructed or partial counts | No counts or depth |
| Beta-binomial sensitivity | As Wilson (labelled scenario values) | As Wilson | As Wilson |
| Empirical calibration | Specimens with distinct documented Library_IDs and counts | Same, with reconstructed counts | No verified replicates |
| Temporal mutation-frequency analysis | Every sample dated (2 or more) | Some dated | Fewer than 2 dated |
| Full-count Yue–Clayton | Per-position A/C/G/T table | — | No per-position counts (sparse base counts do not count) |
| PCoA | Yue–Clayton available, 3+ samples | Exploratory Jaccard only | Fewer than 3 samples |
| Coverage-aware non-detection | Per-position table | — | None |
| Reference-based comparisons | Reference FASTA | Labels only (Ref bases unverified) | — |
| Gene and protein annotation (Phase 4) | GenBank/GFF3 of the same accession.version | — | None or a different reference |
| Experimental grouping | — | — | — (always **Experimental**; off by default) |
| Sprint 2 + 3 cross-source | FASTA + reference + pooled | — | Otherwise |

Results on the real data are in §7.

## 7. Results from each published dataset (real, verified)

Full provenance is in §7.4. Each table lists the readiness actually produced.

### 7.1 Fontenele et al. 2021 (`mmc2.xlsx`, regression dataset)

- **Checks.**
  - 7,973 records and 52 samples imported.
  - Reported frequencies are identical to the workbook (max difference 0).
  - Depth equals TOTAL_DP.
  - 7,973/7,973 Ref bases match NCBI NC_045512.2.
  - Counts are labelled "Reconstructed" (7,973), never measured.
  - md5 unchanged.
- **Readiness.**
  - Available: basic summaries, depth QC, reported VAF, reference comparisons, temporal analysis.
  - Partially available: count-based VAF and Wilson (reconstructed counts), exploratory PCoA.
  - Unavailable: Yue–Clayton, coverage-aware non-detection, empirical calibration.
  - Experimental: grouping.

### 7.2 Crits-Christoph et al. 2021 (mBio, doi 10.1128/mBio.02703-20)

- **Table S4 (`mBio.02703-20-st004.xlsx`, sheet "S4- variants", 293 rows).**
  - Columns: `sample, genome_position, position_coverage, allele_count, ref_base, con_base,
    var_base, ref_freq, con_freq, var_freq, A, C, T, G`.
  - No reference column. The paper states the reference was WIV04 (EPI_ISL_402124).
  - Positions are 1-based: 293/293 Ref bases match WIV04 at position − 1.
  - Base-count layout, derived with the paper's stated minimum of 2 reads:
    - 294 records, equal to an independent derivation;
    - calculated VAF equals the study's own reported allele frequency on 293/293 rows;
    - one extra allele with ≥ 2 reads that the study did not list: `6_30_S_MOS:6257A>G`;
    - `position_coverage` equals A+C+G+T on every row.
- **Table S1** (28 rows, sheet "S1-samples").
  - Sample, run, location, method, date, concentrations, Ct, enrichment, breadth, and coverage.
  - It ends with an "average:" row that has no sample ID; this is reported.
  - **None of the 7 S4 sample names appears in S1** (`5_28_S` against `5_28_I`; `5_19_S` and
    `5_19_S_2` against `5_19_F_S3` and `5_19_A_S4`). The join is refused. No crosswalk was
    invented, so dates stay unknown and temporal analyses are Unavailable.
- **GitHub repository** (`alexcritschristoph/wastewater_sarscov2`, GPL-3.0).
  - `final_tables/Wastewater_SNVs_v3.csv` uses 0-based positions (the README says so). After the
    user's 0-based choice it reproduces Table S4 exactly (294/294).
  - The 8 per-sample `*_SNVs.tsv` files have no sample column and 0-based positions:
    - the importer warned on 8/8 files (16–30% Ref fit as 1-based, 100% as 0-based);
    - after conversion, 250 of 263 sites coincide with S4 (a one-base error would give 25);
    - these files are an earlier inStrain run: 13 sites appear only in them and 43 only in S4.
  - `reference.fasta` (WIV04, 29,840 bases) is a prefix-identical but different sequence from
    NC_045512.2 (29,903 bases). It is never treated as NC_045512.2, so the NC_045512.2
    annotation is not applied.
- **Readiness.**
  - Available: basic summaries, depth QC, count-based VAF, Wilson, sensitivity, reference.
  - Partially available: reported VAF (calculated only), PCoA (exploratory).
  - Unavailable: temporal analysis (no joinable dates), Yue–Clayton (counts only at listed
    sites), coverage-aware non-detection, calibration.

### 7.3 Schumann et al. 2022 (Sci Total Environ, doi 10.1016/j.scitotenv.2022.158931)

- **`mmc1.xlsx` (16 sheets, 4 hidden).** "S2.1 - mutations_found_dataset-" holds the Berlin
  dataset:
  - a title row, then `samplename, dates, location_name, coordinates_lat, coordinates_long`;
  - then 1,907 mutation columns named `GENE:AA::RefPosAlt` (1,901 are single SNVs, 6 are
    multi-allele or invalid headers).
- **Definitions** (PiGx source at the paper's commit `524ed48`):
  - reference `NC_045512.2.fasta`;
  - frequencies are LoFreq AF fractions;
  - `NA` is a mutation not called in that sample (`full_join` of per-sample tables);
  - no depth.
- **Checks.**
  - 5,145 records = non-NA SNV cells; 5,145/5,145 frequencies exact.
  - 67,093 NA cells → no record.
  - No depth invented.
  - Per-sample counts equal the authors' hidden "S4 - mutation_counts" `total_muts` for 38/38
    samples (SNV records + non-SNV columns).
  - Two samples have no SNV and are listed.
  - Dates 36/36.
  - Ref bases match NC_045512.2.
  - **The Phase 4 amino-acid predictions agree with the authors' VEP protein labels for
    1,830/1,830 labelled SNVs.**
- **`mmc2.xlsx`** (Berlin35, 109 samples × 8,862 SNV columns): imported locally. See §13 for
  memory; it is refused under the Render cell limit with an explanation.
- **`mmc3.xlsx`** (UCSD spike-in, 14 MB, 384 × 12,359): refused with "larger than this server can
  process safely (more than 4,000,000 cells)…". It is also above the 10 MB upload limit.
- **Readiness.**
  - Available: basic summaries, reported VAF, temporal analysis, reference, annotation.
  - Partially available: PCoA (exploratory).
  - Unavailable: depth QC, count-based VAF, Wilson, beta-binomial, Yue–Clayton, coverage-aware
    non-detection.

### 7.4 Dataset provenance

| | Fontenele 2021 | Crits-Christoph 2021 | Schumann 2022 |
|---|---|---|---|
| Title | High-throughput sequencing of SARS-CoV-2 in wastewater provides insights into circulating variants | Genome Sequencing of Sewage Detects Regionally Prevalent SARS-CoV-2 Variants | SARS-CoV-2 infection dynamics revealed by wastewater sequencing analysis and deconvolution |
| DOI | 10.1016/j.watres.2021.117710 (PMC8464352) | 10.1128/mBio.02703-20 (PMC7845645) | 10.1016/j.scitotenv.2022.158931 (PMC9549760) |
| Download | Article supplement | Europe PMC supplementaryFiles API; GitHub repository | Europe PMC supplementaryFiles API; PiGx GitHub at commit 524ed48 |
| Files (md5) | `mmc2.xlsx` (c0a75b7f…) | `st001` (6041fcdf…), `st004` (6d3fe78d…), `reference.fasta` (c04b26e5…), 8 `*_SNVs.tsv`, `Wastewater_SNVs_v3.csv` | `mmc1` (d01b7222…), `mmc2` (a01f839d…), `mmc3` (e25a95ed…) |
| Format / sheets | XLSX, one sheet with a title row | XLSX "S1-samples", "S4- variants"; TSV; CSV | XLSX, 16 sheets (S2.1 used; hidden counts sheet for checking); `data_mutation_plot.csv` sheets |
| Samples | 52 | S4: 7; S1: 27 (+ 1 average row); GitHub: 8 | 38 (36 with an SNV); 109; 384 |
| Variant records | 7,973 | 294 (S4, ≥ 2 reads) | 5,145; 24,575 non-NA cells in mmc2 |
| Reference | NC_045512.2 (paper methods) | WIV04 EPI_ISL_402124 (paper); not NC_045512.2 | NC_045512.2 (PiGx default) |
| Depth | TOTAL_DP | A+C+G+T at listed sites | None |
| Allele counts | Reconstructable from frequency × depth | Measured (A/C/G/T) | None |
| Dates / locations | Yes (state: location) | In S1 only; not joinable | Yes, with coordinates |
| Missing | Ref counts, coverage, replicates | Joinable dates, genome-wide coverage, replicates | Depth, counts, coverage |

**Retrieval.** PMC's own download pages answered automated requests with a proof-of-work bot
check, which was not bypassed. The Europe PMC REST endpoint
`https://www.ebi.ac.uk/europepmc/webservices/rest/<PMCID>/supplementaryFiles` returned the same
files. The commands are in `validation/check_published_datasets.py`.

## 8. Numerical validation

- **Independent of application code** (`validation/check_published_datasets.py`; openpyxl/pandas
  reads of the originals): 32 checks across the three studies. Result `FAILED checks: 0`.
- **Synthetic equivalence tests** (`tests/test_phase4a_import.py`). The same six observations, in
  several layouts, give identical canonical values:
  - renamed CSV, TSV, and XLSX;
  - a percent table after conversion: reported VAF equal to 1e-15, and original text kept;
  - a 0-based table after conversion;
  - a header row below notes;
  - a wide matrix: reported VAF equal; depth correctly absent;
  - a base-count table: identical counts and depth; frequency marked calculated.
- **Regression:**
  - 42/42 (`mmc2`), 48/48 (`mmc2` + GenBank), and 64/64 (synthetic example) earlier output files
    are identical, as DataFrames and as PNG bytes;
  - the only additions are `canonical_*` (3) and `import_provenance.csv`;
  - on the web path, `analysis_readiness.csv` and `metadata_matching.csv` are also added.
- **Downloads.** For a Schumann run, the XLSX sheets Observations (5,145 × 23), Samples (36 × 22),
  and Readiness (15 × 5) equal the CSV downloads value for value. The web test repeats this on
  synthetic data.

## 9. Website interface changes

- **Steps in the page header:** 1 Upload files · 2 Review detected data · 3 Review available
  analyses · 4 Run (optional Sprint 4 enrichment) · 5 View and download results.
- **Step 2 (`templates/_review_import.html`).**
  - **"How to read each table":** one compact block per table, showing only the options that apply
    to it.
  - **"Sample metadata matching":** the report, the conflicts table, the resolution choice, and the
    "continue without metadata" box.
  - **"Save or load a mapping profile":** collapsed.
  - **Collapsed:** detected columns, previews, reference and metadata details, the field-by-field
    readiness, and the per-file plan.
- **Step 3 (`templates/_readiness_step.html`).** Status counts, then one row per analysis with a
  badge, the reason, and what it needs.
- **Step 4.** The existing sticky action bar with *Update review* and *Run analysis*. Settings sit
  under "Advanced settings".
- **Results.**
  - A "Download these results" bar above the sections: XLSX, PDF, all CSV/PNG files.
  - Under *Provenance and downloads*, collapsible: the readiness decided before the run, the import
    provenance, and the metadata matching.
- **Unchanged:**
  - the Sprint 2 and Sprint 4 workflows;
  - the legacy Sprint 3 pages;
  - the existing headings tested elsewhere ("Step 2 · Review your dataset", "Analyses available
    with these files").

## 10. Download and export capabilities

| Format | Content | Built by |
|---|---|---|
| CSV | Every result table, including the new canonical, provenance, readiness, and matching tables | `integrated_view.download_files` (unchanged path) |
| PNG | Every figure | unchanged |
| **XLSX** `analysis_results.xlsx` (built on first download) | Summary sheet (findings as on the page, sheet index, limitations), then Readiness, Analyses run, Import provenance, Metadata matching, Samples, Observations, Coverage, Reference check, Paper QC, SNVs per sample, Variant positions, Temporal trends, VAF time series, Uncertainty models, Overdispersion, Mutation eligibility, Annotated mutations, Cross-source, Provenance (when present); tables over 10,000 rows are listed as CSV-only | `genomics/report_exports.build_workbook`, read back from the stored CSVs (openpyxl write-only) |
| **PDF** `analysis_report.pdf` (built on first download) | Disclaimer, dataset, main findings (the dashboard's own sentences), Step 3 readiness, what ran, import choices, limitations, up to 8 main figures (wide figures on landscape pages) | `genomics/report_exports.build_pdf_report` (matplotlib PdfPages; no new dependency) |
| JSON | Mapping profile (no records); `analysis_summary.json` (findings, notes, readiness, figure titles) | `genomics/mapping_profiles`; `report_exports.summary_json` |

- Building on demand keeps the analysis request as fast as in Phase 4 (§13). Values are read back
  from the CSVs with round-trip float precision, so the workbook equals the downloads.
- The report waits if an analysis is running (the same one-at-a-time slot). A failure gives a
  short page saying "the CSV downloads hold the same values"; the analysis is never affected.
- `run_integrated_analysis.py --reports` writes the same two files locally.

## 11. Exact files and modules modified

**New:**
- Genomics modules:
  - `genomics/canonical_model.py`
  - `genomics/table_layouts.py`
  - `genomics/import_options.py`
  - `genomics/metadata_matching.py`
  - `genomics/analysis_readiness.py`
  - `genomics/mapping_profiles.py`
  - `genomics/report_exports.py`
- Templates: `templates/_review_import.html`, `templates/_readiness_step.html`
- Validation: `validation/check_published_datasets.py`
- Tests:
  - `tests/test_phase4a_import.py` (29)
  - `tests/test_phase4a_web.py` (6)
  - `tests/test_phase4a_published.py` (3 optional)
- Docs: `docs/SPRINT5_PHASE4A_PROGRESS.md`, `docs/SPRINT5_PHASE4A_AUDIT_CHECKLIST.md`, this report.

**Modified:**

| File | Change |
|---|---|
| `genomics/input_inspection.py` | Layouts at inspection; `read_raw_table` (header row, file-name sample); `table_for(path, info, role, options)`; `interpretation_review`; `numbering_review`; `interpreted_table`; metadata matching hook; `prepared_pooled_metadata`; coordinates from variant tables; `import_provenance`; `Source_File` per combined record; matched-sample date counts |
| `genomics/column_detection.py` | Optional `Read_Depth` and `Reference_Genome`/`Reference_ID`; new fields (`Ref_Depth`, `Sequencing_Platform`, SNV-role coordinates, positions `Ref`); aliases (§4); `value_warnings(…, frequency_scale)` |
| `genomics/excel_tables.py` | `read_sheet(header_row=…)`, `locate_header` |
| `genomics/file_formats.py` | `read_table(header_row=…)` |
| `genomics/integrated_analysis.py` | `import_details` parameter; canonical tables; readiness and matching tables; `pooled["metadata"]`; `result["input_labels"]` |
| `integrated_view.py` | `download_files` stores `analysis_summary.json`; the dashboard gets report links, findings, readiness, import provenance, matching |
| `app.py` | `build_report_download` (on-demand XLSX/PDF, one at a time, cached); review options (`opt_<i>_*`, conflict resolution, unmatched samples); re-inspection on header-row and file-name changes; mapping-profile load (`use_profile`) and download route `/review/profile`; `import_details` passed to the run |
| `run_integrated_analysis.py` | `--reports` |
| `templates/index.html` | Five step labels, multipart review form, includes, collapsed technical details, results download bar, provenance details, badge styles |
| `genomics/metadata_readiness.py` | Pooled read depth is **Recommended** (it was Required); a stated reference counts as present. This matches the new optional depth: a frequency-only table no longer shows "required information missing" while it can run |
| `README.md` | Formats table rows and a real-data paragraph |
| `tests/test_metadata_readiness.py` | One expectation updated: "Recommended Read depth Detected" (was Required), documented in the test |
| `tests/test_literature_web.py`, `tests/test_upload_enrichment.py` | One expectation each updated for the renumbered step label "4 · Run · Metadata enrichment (Sprint 4, optional)" (was "3 · …"); the Sprint 4 workflow is unchanged |

## 12. Full automated test results

**Command** (project `.venv`; MAFFT on PATH; real files through environment variables):

```bash
FONTENELE_MMC2=<mmc2.xlsx> CRITS_CHRISTOPH_DIR=<folder> SCHUMANN_DIR=<folder> AI_ENRICHMENT_ENABLED=false python -m unittest discover -s tests
```

| Run | Result |
|---|---|
| Baseline at the start of Phase 4A (same command, same real files except the two new folders) | `Ran 636 tests in 448.2 s — OK (skipped=5)` |
| Intermediate run (after Milestones 1–4) | 636 tests, 1 failure: a new "file name as sample ID" hint replaced the existing explanation for a workbook with no frequencies; fixed (§9 of the progress file) |
| Intermediate run (after Milestones 6–7) | 674 tests, 2 failures: two Sprint 4 tests checked the old step label "3 · Metadata enrichment (Sprint 4, optional)"; the label keeps the Sprint 4 wording as step 4 and the two assertions were updated |
| **Final** (final code) | **`Ran 674 tests in 651.6 s — OK (skipped=5)`**: 636 pre-existing + 38 new (29 import, 6 web, 3 published-dataset) |

- **The 5 skips** are the same optional checks as at baseline:
  - 3 SciPy cross-checks (SciPy is not in `.venv`);
  - the NCBI Virus files test;
  - the mpox GenBank test.
- **Real data:** `mmc2.xlsx`, Crits-Christoph, and Schumann ran through `test_phase4a_published.py`
  and `test_real_workbook_exports`.
- **Intentional changes to existing tests** (all documented in the tests themselves):
  - read depth is now *Recommended*, not *Required*;
  - one step label was renumbered.

**What the new tests cover** (goal Milestone 9):
- different column names (CSV, TSV, XLSX), different worksheets, and the header row;
- long and wide layouts;
- separate metadata, and valid and invalid joins (unmatched, near-miss, rows without IDs,
  duplicates, conflicts with all three resolutions);
- ambiguous names ("Site", "sequencing");
- invalid frequencies (above 1 without a choice; above 100 with percent);
- percent against fraction (identical VAF);
- missing depth, missing alternative counts, duplicate rows, missing dates;
- incompatible references (a stated label different from the FASTA, never verified), 0- versus
  1-based positions;
- segmented genomes, sparse SNV tables (no coverage, unknown never 0), full A/C/G/T tables
  (coverage records, Yue–Clayton available);
- verified and unverified replicates, partial compatibility;
- large datasets (a 120 × 1,500 wide matrix under 10 s; a cell-limit refusal);
- downloads (XLSX equal to CSV, a valid PDF), website navigation, mapping-profile round trip,
  conflict choice on the page, session isolation (two uploads keep their own header-row choices);
- the existing Sprint 2 + 3 integration (the full suite, and output stability §8).

## 13. Memory and runtime measurements

Method: the guided web path `/review` → (choices) → `/run`, then a click on each report link.
Windows peak working set, `OPENBLAS_NUM_THREADS=1`, measured with `perf_p4a.py` (session
scratchpad). The Render-like run uses a Windows Job Object limited to 512 MB and 10% of one CPU
core (`render_like.py`, as in Phases 1–4). Phase 4 values are from its report.

| Run | Phase 4 (before) | Phase 4A (final code) |
|---|---|---|
| `mmc2.xlsx`, one run | 295–296 MB, 12.0–12.2 s | **308 MB, 12.2–12.5 s** (canonical tables 0.04 s) |
| `mmc2.xlsx`, four runs in one process | 323 MB | 352 MB (the result cache now also holds the generated workbooks) |
| `mmc2.xlsx`, 512 MB cap, 10% CPU | 288 MB, 124.6 s | **300 MB, 126.4 s**; job peak committed 302 MB |
| On-demand XLSX (mmc2) | — | 4.4–4.6 s, 3.15 MB (42.7 s at 10% CPU) |
| On-demand PDF (mmc2) | — | 0.6 s, 0.75 MB (6.5 s at 10% CPU) |
| Schumann mmc1 S2.1 (38 samples, wide) | could not be imported | 231 MB, 10.9 s; XLSX 1.9 s, PDF 0.7 s (peak 256 MB) |
| Schumann mmc2 Berlin35 (109 × 8,862; 966,000 cells) | could not be imported | 614 MB, 46.4 s locally. **On a 512 MB server it is refused before the run** by `MAX_OBSERVATION_CELLS=700000` (Dockerfile), with an actionable message; nothing is trimmed |
| Schumann mmc3 (4.7 M cells, 14 MB) | — | Refused at inspection (worksheet cell limit, and above the 10 MB upload limit) |
| HTML page (mmc2) | 3.01 MB | 3.02 MB |

**Regression investigated and fixed.** The first Phase 4A version built the XLSX and PDF inside
the analysis request:

- The mmc2 run took 20.7 s and peaked at 365 MB.
- At 10% CPU it took 197.8 s, longer than the 180 s gunicorn timeout in the Dockerfile. On Render
  this would have killed the worker.

The cause was measured with `tracemalloc` and timers:

- XLSX: 6.8 s of openpyxl cell writing, about 1.1 M cells.
- PDF: an 82 MB transient. matplotlib keeps every embedded image until the PDF closes, plus a
  resampled copy of each one.

Fixes:

1. **Build on demand.** Both files are built only when their link is clicked, from the stored
   CSV/PNG downloads plus a small `analysis_summary.json`. This happens one at a time with
   analyses (the same slot), and the result is cached.
2. **Smaller, faster workbook.** Cells are converted in bulk. Detail tables over 10,000 rows (for
   mmc2, the 23,910-row overdispersion sensitivity table) stay CSV-only and are listed on the
   Summary sheet.
3. **Lighter PDF figures.** Figures are embedded as 8-bit images of at most 1,200 px, without
   resampling. The PDF peak fell from 82 to 29 MB.

The remaining single-run difference (+12 MB, about +0.3 s) comes from the canonical tables and
their CSV encoding.

## 14. Scientific limitations

- **The importer reads processed results only.** It cannot correct upstream variant-calling
  choices: thresholds, primer trimming, duplicate removal, reference choice.
- **Base-count derivation.**
  - It turns every non-reference base with at least the user's minimum into an observation. With
    a minimum below the source's own calling rule it can include sequencing errors.
  - The rule should come from the source's methods (Crits-Christoph: 2 reads).
  - Records are labelled calculated, and Sprint 3 depth QC still applies.
- **Wide matrices carry no depth or coverage.** NA conflates "not covered" and "covered but not
  called". Both stay unknown, so absence is never inferred.
- **Stated references are user assertions.** Without a matching FASTA the Ref bases are
  unverified; the review and the outputs say so.
- **The depth definition is recorded, not used.** The Sprint 3 + 5 methods still treat Read_Depth
  as the total usable depth. A source that counts ref + alt only would bias calculated VAF
  slightly upward when other bases exist.
- **Samples are joined only on identical IDs.** Real crosswalks (Crits-Christoph S4 ↔ S1) need
  documented mapping tables from the authors.
- **Readiness predicts; it does not prove.** It is based on mapped columns and fill rates. The
  run still validates every record, and `analysis_status.csv` records what actually ran.
- **QC thresholds.** The Fontenele et al. "paper" QC preset (depth ≥ 20 and so on) is applied to
  other studies by default. It is labelled as that study's preset, and "Exploratory" mode flags
  without excluding.

## 15. Readiness for biological mutation annotation

Phase 4 annotation (GenBank/GFF3, codons, amino acids) now works on any imported layout. The
evidence:
- 5,613/5,613 Fontenele codons;
- 1,830/1,830 Schumann VEP labels, through a wide-matrix import.

What each further extension needs:

| Extension | Data the canonical model already provides | Still required |
|---|---|---|
| Reference-based gene annotation | Reference_ID + 1-based position + Ref/Alt; reference verified by sequence | Annotation of the same accession.version (done in Phase 4) |
| Synonymous / nonsynonymous | As above | Done in Phase 4 for SNVs; indels and MNVs would need the importer to keep them (currently set aside) |
| Amino-acid substitutions | As above; shared-codon flags | Read-level linkage (phased reads) to combine SNVs in one codon |
| Genomic mutation maps | Positions, VAF, depth, sample metadata | Nothing new; figure design only |
| Real-data mutation-group validation | Measured counts (base-count layouts), sample dates | Datasets with known mixtures and replicate libraries (§16) |
| Offline FASTQ preprocessing outputs | iVar / VCF / base-count readers; per-position tables | Running `offline_workflow/` off the web server; per-position A/C/G/T export |
| Other pathogens and segmented genomes | One reference per segment; segments listed and never compared | Validation on a non-SARS-CoV-2 wastewater dataset (influenza, mpox) |

## 16. Requirements for improving experimental MAGOS-inspired clustering

The method stays **disabled by default** and labelled "not validated". Phase 3 failed P1 in 6 of 30
control cells and P2 except without overdispersion. Phase 4A does not change the method. It makes
more inputs usable:

- **Measured counts from base-count layouts.** Crits-Christoph runs now feed the depth-aware
  variant. Frequency-only matrices (Schumann) can feed only ordinary clustering.
- **Readiness says which variant can run.**

Before any promotion, it would need:

1. **Calibrated overdispersion** from real technical replicates: distinct libraries of one
   specimen, with counts. None of the three published datasets documents such replicates in a
   joinable way.
2. **Genome-wide per-position A/C/G/T or depth.** Non-detections become known zeros, which removes
   the "detections only" bias. None of the three studies publishes this table; it needs offline
   reprocessing of SRA reads (for example PRJNA827160 for Schumann).
3. **Ground-truth mixtures with read counts.** For example the UCSD spike-in (Karthikeyan et al.
   2021) reprocessed to counts: Schumann mmc3 has only frequencies. These would test P1/P2 on real
   amplicon noise.
4. **The Phase 3 pre-registered criteria re-run on those data**, with held-out samples, before any
   default change.

The preserved assets are unchanged: the transferable ideas (depth-dependent VAF uncertainty,
extra-binomial variation) and the non-transferred tumor assumptions are as in the Phase 3 report.
Simulations, benchmarks (`docs/sprint5_phase3_results/`), and error-rate disclosures are not
modified.

## 17. Local testing instructions

```bash
.venv\Scripts\activate
set AI_ENRICHMENT_ENABLED=false
python app.py
```

Then open the local address. Use the published files from the retrieval notes (they are not in
Git).

**Schumann S2.1 (wide matrix):**
1. Upload `mmc1.xlsx` with `examples/annotation/NC_045512.2.fasta`, and optionally `.gb`.
2. Choose the worksheet "S2.1 - mutations_found_dataset-" and press *Update review*.
3. Under "How to read each table", tick the layout confirmation and type `NC_045512.2`.
4. Update review → Step 3 → Run.

**Crits-Christoph S4 (base counts):**
1. Upload `mBio.02703-20-st004.xlsx` with the repository's `reference.fasta`, and optionally
   `st001.xlsx`.
2. Set the minimum reads to 2.
3. Set the reference to `hCoV-19_Wuhan_WIV04_2019_EPI_ISL_402124_2019-12-30`.
4. With `st001.xlsx`, the matching report shows 0 of 7 joined. Tick "Continue without metadata"
   (or remove the file).

**Command-line runs:**

```bash
python validation/check_published_datasets.py --output out --fontenele mmc2.xlsx --crits-christoph <dir> --schumann <dir>
```

```bash
python run_integrated_analysis.py --snv mmc2.xlsx --reports --output-folder out_mmc2
```

## 18. Render deployment recommendations (not done; after review only)

1. **Deploy only after your review.** No branch was pushed.
2. **Before deploying,** run the full test suite with the real files, plus the validation script.
3. **On Render,** upload `mmc2.xlsx` twice and Schumann `mmc1.xlsx` (S2.1). Check that Metrics
   shows a peak well below 512 MB and Events shows no "Ran out of memory".
4. **Keep `MAX_OBSERVATION_CELLS=700000`.** Schumann Berlin35 (≈ 966,000 cells) is then refused
   with an explanation. Analyze it with `run_integrated_analysis.py` locally.
5. **The XLSX and PDF are built only when clicked** (42.7 s and 6.5 s at 10% CPU for mmc2, under
   the 180 s worker timeout). Check one download on Render after deploying.

## 19. Remaining limitations and future research requirements

- **The browser preview of the new pages was not visually checked.** The in-app browser refuses
  localhost. The pages were checked through Flask test-client runs (text, headings, all download
  links) and by viewing the PDF pages.
- **Render behaviour is unverified** (no deployment).
- **All three real datasets are SARS-CoV-2 wastewater.** A non-SARS-CoV-2 dataset is still needed
  before claiming general pathogen support (audit checklist §13).
- **Mapping profiles** match FASTA roles by kind and order: a reference needs exactly one record,
  but two individual-sequence FASTA roles are not distinguished. Users see and can correct roles
  on the review page.
- **Crits-Christoph sample crosswalk.** S4 ↔ S1 would need the authors' documented mapping; none
  was invented.
- **No dataset documents verified technical replicates**, so empirical calibration remains
  untested on real data.

## Mandatory completion criteria

| # | Criterion | Status | Evidence |
|---|---|---|---|
| 1 | Existing pipeline inspected, baseline verified | **Met** | Progress §0: 636/5 skips verified (the goal's 585/4 was the Phase 3 baseline); outputs saved for comparison |
| 2 | Reusable canonical sequencing-data model | **Met** | §3; `canonical_samples`, `canonical_variant_observations`, `canonical_coverage`, `import_provenance` |
| 3 | Different column layouts mapped correctly | **Met** | §2, §4, §8; equivalence tests and three real datasets |
| 4 | Ambiguous mappings require confirmation | **Met** | Alias ambiguity (existing) plus wide layout, base-count minimum, percent scale, stated reference, numbering warning (§4) |
| 5 | Multiple-file metadata matching works | **Met** | §5; separate files give identical results to a combined table |
| 6 | Conflicting and missing metadata handled safely | **Met** | §5: explicit resolution or unresolved; unmatched samples blocked until accepted, then unknown |
| 7 | Analysis availability from actual validated inputs | **Met** | §6; four states, no file or paper names used |
| 8 | Original and derived measurements keep provenance | **Met** | Reported versus calculated VAF, `Source_*_Original`, count and frequency provenance, `import_provenance.csv` |
| 9 | Two or more independent published structures end to end, including Fontenele | **Met** | §7: Fontenele, Crits-Christoph (4 structures), Schumann (3 files); 32/32 independent checks |
| 10 | Existing Sprint 3 and Sprint 5 calculations correct | **Met** | §8: 42/48/64 earlier output files identical; all 636 earlier tests pass |
| 11 | Experimental grouping disabled by default | **Met** | Default unchanged and tested; readiness says "Experimental … disabled by default" |
| 12 | Website navigation understandable | **Met (scripted check)** | §9: five labelled steps, only the relevant questions, collapsed details; headless-browser screenshot and test-client flows. A human usability check and a real-browser session remain for the audit (checklist §11) |
| 13 | Results on the dashboard and as downloadable files | **Met** | §10: dashboard, CSV, PNG, XLSX, PDF, JSON |
| 14 | Missing information never becomes fabricated evidence | **Met** | NA → no record; blank depth stays blank; calculated frequencies labelled; no ID crosswalk invented; coverage table never zero-filled |
| 15 | All existing and new relevant tests pass | **Met** | `Ran 674 tests — OK (skipped=5)`, all skips optional |
| 16 | Performance measured, regressions investigated | **Met** | §13: a timeout risk was found and fixed; the final run is within 12 MB and 2 s of Phase 4 |
| 17 | No unrelated Sprint 2 or Sprint 4 functionality broken | **Met** | All Sprint 2 and Sprint 4 tests pass (two label assertions updated); workflows unchanged |
| 18 | No paid API calls, commits, pushes, merges, deployments | **Met** | `AI_ENRICHMENT_ENABLED=false`; `git log` still at `b641533`; nothing pushed or deployed |
| 19 | Implementation report completed | **Met** | This document |

**Not verified (stated, not claimed):**
- hosted behaviour on Render;
- usability with a person who has not seen the application;
- any non-SARS-CoV-2 real dataset.

Scientific validity of the experimental grouping is unchanged: it is **not** validated.

**SPRINT 5 PHASE 4A COMPLETE**
