# Comprehensive Independent Audit: Sprint 3 + Sprint 5 (Phases 1, 2, 3, 4, 4A)

**Audit date:** 2026-11-15 · **Branch:** `sprint4-integration` · **HEAD:** `b641533` with the uncommitted
Sprint 5 work in the working tree.

**This was an audit only.**
- No application code, test, or existing document was changed.
- Nothing was committed, pushed, merged, or deployed.
- No paid API call was made: every run used `AI_ENRICHMENT_ENABLED=false`, and `.env` was never read.
- The original published files were only read. Their md5 checksums are unchanged (§5).

**How evidence is labelled:**
- **Verified** — reproduced in this audit with a command or script.
- **Reported** — stated by an earlier implementation report and not re-checked here.
- **Inferred** — a conclusion from reading code, not from a run.

**Where the audit scripts are.** All audit scripts and outputs are in the session scratchpad
`C:\Users\ALEXAN~1\AppData\Local\Temp\claude\a7271a137b0bdb68\922255f3-0871-46d4-aedd-17d549d92b7a\scratchpad\audit\`.
That folder is written `<audit>` below. It is temporary storage, so copy it somewhere durable if the
evidence must be kept. Each subfolder holds the scripts and outputs for one area:
- `imports/` — import scenarios;
- `exports/` — website outputs;
- `cc_ids/` — Crits-Christoph identifiers;
- `perf/` — memory;
- `math/` — independent scientific checks;
- `pathogen/` — other pathogens;
- `security/` — security and reliability.

**Who did the work.** The scientific, pathogen, and security checks were done by three separate
audit agents. Each was given the code, not the earlier conclusions. Their key claims were spot-checked
by re-running their scripts (§4, §7, §12).

---

## 1. Executive summary

**Verified working:**
- **Full test suite:** `Ran 674 tests in 652.3 s — OK (skipped=5)`. The run used real MAFFT and all
  three published datasets.
- **Published-data validation:** 33 of 33 checks pass for Fontenele 2021, Crits-Christoph 2021 and
  Schumann 2022.
- **Fontenele 2021, independently:**
  - 7,973 records and 52 samples import unchanged; reported frequencies match to within 5×10⁻⁷.
  - All 7,973 reference bases match NC_045512.2.
  - The app's amino-acid predictions agree with the authors' own labels for 5,613 of 5,613
    amino-acid pairs and codons.
- **Statistical methods:** Wilson, beta-binomial, the overdispersion estimator, chi-square, PCoA and
  Yue–Clayton match independent implementations to 10⁻⁹ or better. Yue–Clayton is the exception at
  about 10⁻⁵, because of float32 storage.
- **Annotation:** about 181,000 checks of SNV consequences gave 0 mismatches. They covered NC_045512.2,
  synthetic genomes, real mpox, and real H5N1 influenza.
- **Exports:**
  - The workbook equals the CSV downloads value for value (17 sheets).
  - The PDF text matches the summary JSON.
  - Experimental grouping stays off by default.
- **Missing data:** missing data never becomes 0% anywhere in the scenarios tested.

**Problems found.** No P0 (critical) issue was found. Four P1 issues must be fixed before hosted
validation on a public URL:
1. **A small upload can exhaust server memory during the review step.**
   - A 0.2 MB `.csv.gz` reached a 980 MB peak, and a 4.5 MB wide CSV reached 1,158 MB.
   - File inspection runs outside the one-analysis-at-a-time slot, and text tables have no cell limit.
2. **Spreadsheet formula injection.** A sample ID such as `=HYPERLINK(...)` becomes 14 live formula
   cells in `analysis_results.xlsx`, and stays unescaped in 16 CSV downloads.
3. **"Ready" at review, then failure at run, and the uploads are deleted.**
   - This blocks a real published workflow: Crits-Christoph Table S4 together with the authors' own
     sample metadata file. That file has month-day dates with no year.
   - Other triggers: alternative reads above depth, ALT = N, indel-only or header-only iVar files.
4. **Segments can be silently merged.** In a segmented-genome table, stating one reference name puts
   every segment's records on that one reference.

**Also confirmed (P2):**
- A misleading gene track in the genome-position figure.
- The influenza gene "NA" becomes a blank cell in the XLSX.
- An undeclared genetic code silently defaults to table 1.
- Incomplete validation of mapping profiles.
- Other visitors can evict a user's cached results.

**Not verified:**
- Behaviour on Render itself, which was not tested.
- Usability with real people.
- Any real non-SARS-CoV-2 dataset other than one H5N1 influenza set (§7).

**Conclusion at audit time: READY AFTER SPECIFIC REPAIRS** (§18). The four P1 repairs were then made and
re-verified (§20, 685 tests pass), so the current verdict is **READY FOR HOSTED VALIDATION**.

---

## 2. Architecture and execution flow (verified by reading code and by tracing runs)

```
Browser
 │  POST /review (multipart files)                         app.py:1464 review_upload
 ▼
stage_uploads ─ plain_copy (gzip/zip limit 60 MB) ─ detect_kind ─ inspect_file
 │   app.py:1300          genomics/file_formats.py:93            genomics/input_inspection.py:436
 │   → suggest_roles (input_inspection.py:478) → STAGED[token] (in memory, 30 min, at most 6 sets)
 ▼
review_dataset (input_inspection.py:630)   ← runs on every /review and again on /run
 ├ column_detection.detect_columns (:167)          aliases, ambiguity, Required/Recommended
 ├ table_layouts (wide matrix, A/C/G/T counts)     parse_mutation_header, wide_to_long, derive_from_base_counts
 ├ import_options.interpret_table / frequency_problems / reference_check / numbering_warning
 ├ metadata_matching.match_metadata (:82)          identical Sample_ID only; conflicts need a choice
 ├ metadata_readiness + analysis_readiness (:105)  Available / Partially available / Unavailable / Experimental
 └ mapping_profiles (optional JSON, applied only where it fits)
 ▼
POST /run  app.py:1758 run_reviewed → check_settings → prepare_inputs (input_inspection.py:1350)
 ▼
run_analyses (ANALYSIS_SLOT, app.py:922; one analysis per server) → run_genomic_analyses
 ▼
genomics/integrated_analysis.py:176 run_integrated_analysis
 ├ Sprint 2: align_sequences (MAFFT), compare_genetic_distances, clusters (when FASTA sequences are given)
 ├ Sprint 3: wastewater/processed_snv.py:125 load_processed_snvs (validation; 0% rows set aside)
 │           pooled_observations (Duplicate_Status/Counted_Once :210, observation matrix)
 │           quality_control (paper preset, reference verification), temporal_analysis,
 │           diversity_analysis (Yue–Clayton, PCoA, or exploratory Jaccard PCoA)
 ├ Sprint 5 P1: variant_evidence + depth_aware_analysis (VAF provenance, Wilson intervals)
 ├ Sprint 5 P2: overdispersion + overdispersion_analysis (replicate-calibrated beta-binomial)
 ├ Sprint 5 P3: mutation_profiles, mutation_grouping/mutation_group_analysis (grouping off by default)
 ├ Sprint 5 P4: mutation_annotation.run_mutation_annotation (:469) → genome_annotation, snv_consequences
 └ Sprint 5 P4A: canonical_model.add_canonical_tables (:529), readiness and matching tables
 ▼
integrated_view.download_files (:742) → app.store_result_files (zlib, in memory, 10 min, at most 4 runs)
 ▼
Dashboard (templates/index.html) → GET /results/<id>/<file>
   XLSX/PDF built on demand: report_exports.build_from_stored (:380), app.build_report_download (:666)
```

Uploaded files are deleted after `/run`. This happens even when the run fails (app.py:1801-1804).

---

## 3. Verified functionality

| Area | Result | Evidence |
|---|---|---|
| Full test suite (real MAFFT + three real datasets) | 674 run, 0 failures, 5 skipped | `<scratchpad>/audit_full_tests.txt`: "Ran 674 tests in 652.331s / OK (skipped=5)" |
| Optional mpox annotation test (normally skipped) | Passes with real NCBI NC_063383.1 | `ANNOTATION_MPOX_GENBANK=… python -m unittest tests.test_genome_annotation` → 18 OK |
| Published-dataset validation script | 33 PASS, 0 FAILED | `<audit>/published_run.txt` |
| 30 adversarial import scenarios | 26 behave safely as documented; 4 show "ready, then fail" (§6) | `<audit>/imports/audit_import.py`, `out.txt` |
| Website run of `mmc2.xlsx` + FASTA + GenBank | 200, 56 download links; dashboard, CSV, XLSX and PDF agree | `<audit>/exports/audit_exports.py`, `out.txt` |
| Render-like run of `mmc2` (512 MB, 10% of one core) | 120.6 s, peak 301 MB; XLSX 41.5 s; PDF 5.9 s | `<audit>/render_mmc2.txt` |
| Experimental grouping off by default | Confirmed in code (integrated_analysis.py:77) and in every run's status line | All runs |
| Missing data never becomes 0% | Confirmed: NA cells give no record, unreported variants are "coverage unknown", and blank depth stays blank | §6, §4.5 |

---

## 4. Scientific calculation validation

Every reference implementation was written independently. SciPy is not installed, so these used exact
`Fraction` arithmetic, `lgamma`/`erfc`, closed forms, and hand-written codon tables. Scripts and outputs
are in `<audit>/math/t01…t11*.py` / `*_output.txt`. Each output ends with a SUMMARY line.

| # | Method | Code | Verdict | Largest difference / note |
|---|---|---|---|---|
| 1 | VAF, read depth, base-count derivation | table_layouts.py:135-203, pooled_observations.py:241-264, variant_evidence.py:84-345 | **Correct, with caveats** | 44 checks exact. Caveats: a Ref of `N` is not rejected in base-count derivation (table_layouts.py:154,165). One cell can get two VAF denominators (the table's depth vs A+C+G+T, t01c) |
| 2 | Wilson interval | frequency_uncertainty.py:38-60 | **Correct** | 2×10⁻¹⁶; Newcombe (1998) Table I reproduced |
| 3 | Beta-binomial pmf, variance, intervals | beta_binomial.py:92-231, 290-309 | **Correct** | pmf exact against rational arithmetic; variance = n p(1-p)[1+(n-1)ρ] to 4×10⁻⁹ |
| 4 | Overdispersion estimate and replicate eligibility | overdispersion.py:124-415 | **Correct, with caveat** | Method-of-moments root to 10⁻¹⁶. True ρ of 0.002/0.01/0.05 recovered as 0.00215/0.01004/0.0507. All eligibility rules behave as documented (9 sites → "too few", 10 → estimated). Caveat: under-dispersed replicates give an upper bound of exactly 0, and the "calibrated" interval then equals Wilson (l.268-273) |
| 5 | Temporal trends, non-detection | temporal_analysis.py:43-273, vaf_trends.py:50-156 | **Correct, with caveat** | Exact on a 7-sample series with every status. A missing row is never 0 unless counts show coverage. Inconsistency: 3 points on 2 dates give a slope in `temporal_analysis.weekly_slope` (l.144-149) but NaN in `vaf_trends` (l.146) |
| 6 | Yue–Clayton | diversity_analysis.py:43-100 | **Correct, with caveat** | 9×10⁻⁶ (float32 storage, l.63; 1.5×10⁻⁴ relative over 30,000 near-identical positions). Distances are **summed** over shared positions, so they grow with shared coverage. Whether this matches the 2021 paper's convention was not verified |
| 7 | PCoA | diversity_analysis.py:117-148 | **Correct, with caveat** | Equals PCA (SVD) on Euclidean data to 1.4×10⁻¹⁴. "Explained" is a share of the positive eigenvalues only, which should be stated in the label |
| 8 | Chi-square (central and noncentral, power) | chi_square.py:66-128 | **Correct** | 5×10⁻¹⁴ relative |
| 9 | Reference verification and numbering check | references.py:34-127, import_options.py:121-160, quality_control.py:94-168 | **Correct** | Also verified end to end: an undeclared 0-based table gives "0% vs 100%" warning → "Verified mismatch"; QC keeps 0 of 6 records (`imports/out.txt` scenario 06) |
| 10 | Annotation, codon translation, Syn/Nonsyn | genome_annotation.py:355-477, snv_consequences.py:77-333 | **Correct consequences; one bookkeeping bug** | 0 mismatches across 1,257 synthetic SNVs (minus strand, −1 frameshift, overlapping frames, table 2) and an 89,709-SNV sweep of NC_045512.2 (GenBank and GFF3). **Bug:** a CDS with no `/protein_id` or ID gets an empty CDS_ID. It then loses the same-codon flag (snv_consequences.py:231) and is dropped from the variant summary (l.310). Reproduced: `t10_output.txt` lines 30-31 `[BAD]` |
| 11 | Experimental mutation-frequency grouping | mutation_grouping.py:149-367 | **Software correct on constructed cases; scientifically NOT validated** | Complete linkage, ARI/NMI and the ρ=0 χ² reduction are exact. Planted groups recovered (ARI 1.0). The project's own pre-registered criterion was not met in every condition (mutation_group_analysis.py:11-16) |

**Fontenele real-data biology check (independent of app logic).** The `mmc2.xlsx` columns
`REF_CODON`, `REF_AA`, `ALT_AA` and `Sy/Ns` were compared with the app's NC_045512.2 GenBank
annotation:
- 5,613 of 5,613 amino-acid pairs and reference codons agree.
- Every source "Nsyn" is predicted non-synonymous: 4,368 of 4,368.
- Of the 1,309 source "Syn":
  - 1,243 are synonymous;
  - 64 are untranslated-region positions the authors labelled "Syn";
  - 2 are "Stop retained", that is, synonymous at a stop codon.

Reproduce with `<audit>/exports/audit_exports.py`, then the comparison command in this audit's log.

**Software correctness versus scientific validity:**
- Wilson intervals model read sampling only.
- The beta-binomial adds library-to-library variation only when Library_IDs really mark separately
  prepared libraries.
- Primer bias, mapping error and contamination are not modelled.
- Trends are descriptive least-squares slopes, not tests of change.
- An annotation predicts an amino-acid change, not a phenotype.

**Are the tests independent?**
- **Independent (hard-coded external values):**
  - Wilson against Newcombe;
  - beta-binomial against gamma/beta-mixture forms;
  - Yue–Clayton and slopes derived by hand;
  - SARS-CoV-2 known answers (D614G, P323L, R203K/G204R, C241T);
  - every NCBI protein reproduced.
- **Partly self-referential:** `test_overdispersion.py::test_expectations_hold_exactly` (uses the
  project's own `pmf_all`) and the grouping simulations (use the project's own generative model).
- **Silently skipped without SciPy:** three SciPy cross-checks.
- **Untested:** the empty-CDS_ID path, the slope-rule inconsistency, float32 precision, and
  under-dispersed replicates.

---

## 5. Published dataset comparisons

### Fontenele et al. (2021): verified
- **Records:** 7,973 imported, equal to the workbook rows; 52 samples.
  - An independent openpyxl join on sample|position|alt matches every record.
  - Three repeated records are counted once and flagged, so 7,970 observations are counted.
- **Values:**
  - Reported VAF is kept unchanged.
  - `VAF_Used` differs from `ALT_FREQ` by at most 5.0×10⁻⁷. It comes from reconstructed counts, and
    the gap is the source's 7-digit rounding.
  - Depth equals `TOTAL_DP` exactly.
  - Alternative counts are labelled "Reconstructed" in all 7,973 rows.
- **Reference:** NC_045512.2 "Verified compatible"; 7,973 of 7,973 bases match.
- **Exports:** the XLSX equals the CSVs (17 tables, 0 differences); the PDF findings match
  `analysis_summary.json` (§10).
- **md5:** `c0a75b7fa5bde8db09746613354a9f64`, unchanged.

### Crits-Christoph et al. (2021): verified, with a correction to the earlier report
- **Table S4:**
  - 294 records equal an independent derivation (≥2 reads).
  - Calculated VAF equals the published frequency for 293 of 293 rows.
  - WIV04 (29,840 bases) is kept separate from NC_045512.2.
  - The 0-based GitHub table plus 1 reproduces S4 exactly.
- **Why no IDs match Table S1.**
  - S4 uses per-site/date IDs: `5_19_S`, `5_19_S_2`, `5_28_S`, `6_09_S`, `6_30_S_COL`, `6_30_S_MOS`,
    `MR_7_1_MOS`.
  - S1 and the repository's `data/sample_metadata.tsv` use library IDs: `5_19_F_S3`, `5_19_A_S4`,
    `MR_1_7_1_MOS`, and so on.
  - The inStrain `attributes.tsv` files show the S4 profiles were built from BAMs in a
    `combined_good/bams/<ID>` folder. This suggests the libraries were combined.
  - The article text and the README contain **no** documented crosswalk from combined IDs to library
    IDs. No association was inferred.
- **Correction to the Phase 4A report.** The authors' repository does contain a second metadata file,
  `data/wastewater/sample_metadata.txt` (md5 `340139efce8687709b3017ffba8bfc96`, downloaded to
  `<audit>/cc_ids/`).
  - Its IDs match S4 exactly: 7 of 7, plus `SQ_7_1_MOS`, which has no S4 records.
  - It gives location, and a date as month-day only (`5-19`), with no year.
  - The earlier statement "no sample matched … dates stay unknown" is therefore incomplete.
  - Uploading this file exposes P1-3: the review shows "Collection dates … Invalid: 8 — Warning" with
    Run enabled. The run then fails with "294 metadata collection date(s) are invalid" and the uploads
    are deleted. Reproduce with `<audit>/cc_ids/web_ww_meta.py`.
- **md5 unchanged:** st001 `6041fcdf…`, st004 `6d3fe78d…`, reference `c04b26e5…`.

### Schumann et al. (2022): verified
- **Records:** 5,145 records, one per non-NA single-nucleotide cell, all exact; 67,093 NA cells give
  no record.
- **Samples:** 38 of 38 per-sample totals equal the authors' hidden `total_muts` sheet. Two samples
  have no SNV and are reported.
- **Annotation:** 1,830 of 1,830 VEP amino-acid labels agree.
- **Large table:** UCSD mmc3 is refused with an actionable message.
- **md5 unchanged:** mmc1 `d01b7222…`, mmc2 `a01f839d…`, mmc3 `e25a95ed…`.

**Schumann Berlin35 memory (614 MB reported).** Measured here at **546 MB peak and 25 s** in a direct
run (`<audit>/perf/berlin_plain.py`); the 614 MB figure included the web page.
- The dominant object is `pooled_observation_matrix`: 965,958 rows × 14 columns, **533 MB**.
  - Rows = 109 samples × 8,862 mutations.
  - 941,142 rows (97.4%) are "Coverage unknown" placeholders; only 24,816 hold a reported value.
- **Safe optimization (verified, not applied):** storing its repeated text columns as pandas
  categoricals reduces the table to **50 MB**. Values are identical and the CSV output is
  byte-identical (`csv bytes identical: True`). No calculation changes.
- A second option is to not materialize "coverage unknown" rows when no coverage table exists. It
  needs more care, because downstream code reads the dense matrix.
- **Discarded measurement:** a first profiler that sampled stack frames 20 times a second inflated
  the process to 16 GB and 628 s. That measurement is invalid and is not used.

---

## 6. Metadata mapping and compatibility

All results below come from the app's own `review_dataset` → `prepare_inputs` →
`run_integrated_analysis` path, using `read_file_options` for the form options
(`<audit>/imports/out.txt`).

| Scenario | Outcome | Assessment |
|---|---|---|
| Long table; VCF-like names (`POS`, `REF`, `ALT`, `AF`, `DP`, `CHROM`, `sample`) | Mapped correctly. `AD` (a "ref,alt" list in VCF) is correctly **not** mapped to alt depth | Safe |
| Percent 0–100 without a scale choice; 1.5 declared as a fraction | Refused with a clear request to choose percentages | Safe |
| Percent declared | Exact decimal ÷100; originals kept | Safe |
| Percent values that are all ≤ 1 (e.g. 0.5%), no choice | Cannot be detected: read as fractions, 100× too high. A header like `AF_percent` is not recognized, so the table is refused | **P3 limitation.** Document it, or ask the scale whenever a header mentions "percent" |
| Column `Coverage` holding genome breadth (0.98) | Auto-mapped as Read_Depth with only a note; the run then fails validation | P3. Breadth stored as a percent (e.g. 98) would be accepted silently as depth |
| 0-based positions with a FASTA, undeclared | Warning (0% vs 100%); every record fails the reference check | Safe |
| 0-based positions with no FASTA | No check possible; status "Supplied but unverified" | Expected limitation, disclosed |
| FASTA with the same ID but a different sequence | "Verified mismatch"; 1 of 6 kept | Safe |
| Table reference label ≠ FASTA ID (MN908947.3 vs SYN.1) | Warning: "never assumed to be the same genome"; unverified | Safe |
| Two SNV files on two references | Refused: only one reference FASTA, and repeated sample IDs | Safe |
| Exact or conflicting duplicate observation rows | First row counted; "Repeated record, conflicting measurements" flag and a warning. The QC line still says "7 of 7 records retained" | P3. The conflict is disclosed, but the first value is used without a user choice |
| Frequency inconsistent with counts (0.5 vs 10/100) | Reported value used; `Frequency_Check` = "Differs from counts by +0.4000". No warning at review | P3 |
| Separate metadata: exact match, near-miss IDs (`s1`, `S-3`), duplicate metadata rows, variant table vs metadata conflicts, unmatched samples | Exact joins only; near-misses are hints only; conflicts block until metadata / variant_table / unresolved is chosen; unmatched samples block until accepted, then dates stay unknown | **Safe, as designed** |
| Ambiguous names (`Date` + `sample_date`) | Must be chosen by the user | Safe |
| Ambiguous day/month (`03/04/2026`) | Kept as text, not used, with a warning | Safe |
| Wide matrix: unconfirmed; blank, NA or 0 cells | Must be confirmed. Blank/NA give no record. An explicit 0 is set aside as "not a variant call", with a warning | Safe (conservative) |
| A/C/G/T counts (minimum 2) | Values equal a hand derivation | Correct. P3: the same file's per-position depths are **not** used for coverage-aware non-detection unless it is also uploaded as a positions table |
| Same sample ID in two files | Refused | Safe |
| **Alt depth > depth; ALT = N; indel-only or header-only iVar** | Review says ready → run fails → uploads deleted (`audit/web/ready_then_fail.py`) | **P1-3** |
| **Separate metadata with partial dates** | Review "Invalid: 8 — Warning" → run fails, uploads deleted | **P1-3** |
| **Segment column named `Segment`/`seqid`/`chromosome` + one stated reference** | All segments put on one reference: "27/36 match", 15 conflicting duplicates counted once (`pathogen/test_web.py` W3) | **P1-4** |
| Mapping profile | Applied only to files it fits. But string entries, a bad role, or `header_row: -5` are only validated at the top level. `-5` is accepted silently and a data row is used as the headers | **P2** |

---

## 7. Cross-pathogen generality

Evidence: `<audit>/pathogen/` (scripts, `out/test_*.txt`, md5 lists in `downloads/`).

**What is fixed in the code.**
- Genome length, CDS, strand, phase, `transl_table`, segments and figure axes come from the uploaded
  FASTA/GFF3/GenBank. Nothing in the analysis code assumes 29,903 bases or SARS-CoV-2 genes.
- NC_045512.2 literals exist only for recognizing the Fontenele table format:
  input_inspection.py:112,320-322; processed_snv.py:104-109,161.
- The default QC preset is the 2021 SARS-CoV-2 study's (quality_control.py:29-36,271). It is labelled
  as such, but it is applied to any pathogen.

**Real non-SARS-CoV-2 data (verified): H5N1 influenza from dairy cattle.**
- Source: andersen-lab/avian-influenza, commit `e95ab749…`: 32 iVar tables plus the reference
  (md5 `c22951ae…`, identical to NCBI PP755589.1–PP755596.1).
- **Reference:** the 8-segment FASTA is refused ("exactly one sequence"). With HA as the reference,
  52 of 52 HA reference bases verify.
- **Annotation:** iVar's REGION label (`HA|PP755589.1|…`) does not equal the GFF seqid, so annotation
  is refused, safely. There is no way to map the label.
- **After renaming the seqids** (a derived file):
  - 12 of 12 CDS are supported, including spliced M2/NEP, the PA-X frameshift and PB1-F2.
  - 127 of 127 SNVs agree with iVar's codons.
  - An independent Biopython check gave 0 differences.
- **Real references with synthetic observations:**
  - mpox NC_063383.1: 179 CDS, 96 on the minus strand; 40 of 40 correct.
  - HIV gag-pol, tat and rev, and all 11 RSV-A CDS: supported.
  - Measles V (RNA editing), HIV vpr and nef: refused with a stated reason.
- No documented real processed variant table was found for mpox, measles, RSV or dengue.

**Support matrix**

| Item | Status |
|---|---|
| SARS-CoV-2 | Supported; three published datasets verified |
| Influenza A (segmented), real observations | **Partial**: one segment per run; manual seqid renaming needed for annotation; segment-merge risk (P1-4); the gene "NA" is blanked in XLSX (P2) |
| Mpox, HIV-1, RSV, measles | Annotation verified on real references; **no real observations tested** |
| Multiple chromosomes / plasmids / segments in one run | **Unsupported** (one-sequence reference); coverage and PCoA are unavailable for multi-segment samples |
| Bacteria-like 200 kb genome | Works (0/135 mismatches, 4–21 s). 10 MB upload limit. GTG/TTG start codons are misreported as "Start lost" (snv_consequences.py:80) |
| Minus strand; spliced/joined CDS; ±1 ribosomal slippage | Supported |
| RNA editing, translation exceptions, CDS spanning the origin | Refused with a stated reason |
| Alternative genetic codes | Supported when `transl_table` is declared. **Silently table 1 when missing** (verified: code-4 TGG>TGA reported as "W11\* Stop gained") |
| GenBank (multi-record), GFF3 | Supported; GFF3 `##sequence-region` IDs are not URL-decoded, which gives phantom records |
| GTF | Unsupported (generic refusal message) |
| iVar TSV, Freyja variants TSV, LoFreq/VCF | Supported, with the "ready, then fail" gaps (P1-3); ALT = N refuses the whole file |
| Freyja demix output | Correctly refused |

**Conclusion.** General pathogen compatibility **must not be claimed**. Real observations were
tested for SARS-CoV-2 and one H5N1 dataset only.

---

## 8. Experimental clustering limitations

- **Off by default.** It is off by default (integrated_analysis.py:77) and labelled "Experimental …
  disabled by default" on the review page, the dashboard, the XLSX readiness sheet and the PDF.
  Verified in every audit run.
- **Software correctness** was verified on constructed cases (§4.11).
- **On the project's own synthetic mixture:**
  - Depth-aware mode found 2 of 3 planted groups (purity 1.0, ARI 0.971) and grouped 0 of 10 unrelated
    mutations.
  - Ordinary clustering grouped 4 of 10 unrelated mutations.
  - The simulator is the project's own, so this is not independent evidence.
- **The method is not scientifically validated.**
  - It did not meet its own pre-registered criterion in every condition.
  - It has never been compared with sequencing truth on real wastewater.
  - Groups are not lineages.

---

## 9. Website usability findings

**Evidence.** The real review and results pages were rendered in headless Edge
(`<audit>/exports/mmc2_web/review.png`, `dashboard.png`), and the HTML was read through the test
client.

**What the pages do well:**
- Five labelled steps.
- A sticky "Update review / Run analysis" bar.
- Detected roles shown with evidence.
- Collapsed technical details.
- Required, Recommended and Optional explained in words.
- A download bar at the top of the results.
- The overview and main findings come before the figures.
- Unavailable analyses say why ("Needs individual sequences…").

| # | Finding | Severity |
|---|---|---|
| U1 | A run that fails after "ready" deletes the files and shows a validation message mid-page; the user must re-upload (P1-3) | P1 (part of P1-3) |
| U2 | The review page has two readiness systems (Sprint 3 "Metadata readiness" and Phase 4A "Step 3 · Analyses these data support"), one after the other, partly repeating each other. Step 3 sits inside the Step 2 page, while the step strip shows it as a separate step | P3 |
| U3 | The review page is about 3,000 px tall at desktop width. The Step 3 table sits below several sections, though the sticky bar keeps Run reachable | P3 |
| U4 | The results HTML for `mmc2` is 3.18 MB, and the empty "Cross-source evidence" section takes a full card | P3 |
| U5 | The download bar says "The CSV files keep every row", but `pooled_observation_matrix.csv` (over 8 MB) is in neither the CSV links nor the workbook. The note says "also missing", which hints at this but contradicts the sentence above it | P3 |
| U6 | The review flags dates as "Invalid … Warning" but leaves Run enabled, even though the run will then fail (P1-3) | P1 (part of P1-3) |
| U7 | Raw codec or exception text sometimes reaches users ("utf-8 codec can't decode…", input_inspection.py:466-467). Windows-1252 and UTF-16 CSVs are refused | P3 |

**Human usability checklist (not done; scripted checks are not human validation).** Five to eight
participants who have not seen the app, think-aloud, about 30 minutes each. Record time, errors and
help requests.
1. Without help, upload `mmc2.xlsx` plus the NC_045512.2 FASTA. Can they say what each file was
   recognized as?
2. On the review page, find which analyses will **not** run and explain why, in their own words.
3. Given a wide matrix, find and complete the "confirm layout" and "reference" questions.
4. Given metadata with one conflicting date, resolve it and say which value the run will use.
5. After the run, find (a) the number of samples kept after QC, (b) one mutation's uncertainty
   interval, and (c) the download of all results as one workbook.
6. Explain the meaning of "coverage unknown" versus 0%, and of "Experimental".
7. Rate on a 1–5 scale: knowing which files to upload; clarity of the next action; trust in the
   warnings; amount of scrolling.

Success criteria to agree in advance: for example, ≥80% complete tasks 1–5 unaided, and no participant
reads "coverage unknown" as absence.

---

## 10. Export and report verification

Website run of `mmc2.xlsx` + NC_045512.2 FASTA + GenBank; script `<audit>/exports/audit_exports.py`.

| Check | Result |
|---|---|
| Dashboard ↔ CSV | 7,973 / 52 / 7,970 appear consistently; findings equal `analysis_summary.json` |
| XLSX ↔ CSV | 18 sheets; all 17 tables linked to a CSV compared cell by cell: **0 differences**. Large tables are listed as CSV-only on the Summary sheet |
| PDF ↔ summary | 10 pages. Every main finding, the readiness list, the reference status ("All 7973 record(s) match") and the disclaimer are present in the PDF text |
| Identifiers | Sample IDs and `Variant_Key` (`NC_045512.2:pos:R>A`) preserved; join to the original is complete |
| Reference provenance | Reference ID, length and verification status appear in the dashboard, CSV, XLSX and PDF |
| Original vs derived counts | `Count_Provenance=Reconstructed`, `Frequency_Provenance=Reported by the source` in all 7,973 rows |
| Uncertainty labels | "binomial intervals are the baseline"; overdispersion "scenarios, not estimates" |
| Missing-data explanations | "coverage unknown (never 0%)" in the readiness list, the PDF and the dashboard |
| Experimental warnings | Present on all outputs |

**Problems found:**
- **P2, NA blanked in XLSX.** The workbook is rebuilt from the CSVs with pandas' default missing-value
  parsing (`report_exports.py:146`, `read_csv` without `keep_default_na=False`). Text such as the
  influenza gene **"NA"** therefore becomes a blank cell. So would sample IDs such as "NA", "None" or
  "null". The CSVs are correct. Reproduce with `pathogen/test_na_export.py`.
- **P2, misleading gene track.** The genome-position figure draws each source label from its lowest
  to its highest position (`integrated_analysis.py:153` `annotation_ranges`). Fontenele labels scattered
  intergenic positions "UTR" (25,385–29,557), so one band labelled "UTR" covers ORF3a through N. This
  happens even though a GenBank file was uploaded (`mmc2_web/genome_positions.png`).
- **P2, formula injection** (P1-2 in §14).
- **P3, misleading class label.** `Measurement_Class = "Read counts and depth"` for counts that are
  reconstructed. The readiness text and the PDF say "approximate", but the main finding ("7970 of 7970
  … have a read-sampling uncertainty interval") does not.
- **P3, PCoA percentages.** These are shares of positive eigenvalues; the label does not say so.

---

## 11. Performance and Render readiness

| Item | Finding |
|---|---|
| Dockerfile | python:3.11-slim (floating tag); MAFFT; non-root; gunicorn 1 worker, gthread 4 threads, timeout 180; `MAX_OBSERVATION_CELLS=700000`; single-thread BLAS |
| Limits | 10 MB per request (app.py:79), 20 files, 60 MB decompressed (gzip/zip only), XLSX 4 M cells (not lowered in the Dockerfile) |
| Analysis | One at a time (`ANALYSIS_SLOT`); the slot is released on errors; others get a "busy" message and keep their files (verified, 4 concurrent users) |
| **Review / inspection** | **Outside the slot. No cell limit for text tables. A 0.2 MB gz reached a 980 MB peak; a 4.5 MB wide CSV reached 1,158 MB (P1-1)** |
| Caches | Staged uploads: at most 6 sets / 30 min. Results: at most 4 runs / 10 min, in memory (zlib). Count-based eviction lets other visitors remove a user's results (P2) |
| Temporary files | Upload folder removed after `/run` or on expiry; `TemporaryDirectory` during the run |
| `mmc2`, 512 MB / 10% CPU (measured this audit) | **120.6 s, peak 301 MB** (job committed 303 MB). On-demand XLSX 41.5 s, PDF 5.9 s, both under the 180 s timeout. Run headroom about 33% |
| Schumann S2.1 | Reported 231 MB, 10.9 s (not re-measured) |
| Schumann Berlin35 | 546 MB peak locally (this audit). Refused on a 512 MB server before the run by `MAX_OBSERVATION_CELLS`; a 10× reduction is available (§5) |
| GeoNames | One request per place name, 8 s timeout, no cap, while the slot is held (P2) |
| Render | **Not tested.** The Windows Job Object only imitates Render; Render's shared CPU and memory accounting may differ |

**Needing more than Render Free:**
- matrices above 700,000 sample×mutation cells (e.g. Berlin35);
- XLSX sheets near 4 M cells (527 MB peak);
- large FASTA sets for MAFFT/pairwise comparison, which are O(n²) and bounded only by the 150 s MAFFT
  timeout;
- raw-read processing, which is offline by design.

---

## 12. Security and reliability

Evidence: `<audit>/security/*.py` and `*_output.txt`. All runs were offline: sockets blocked, the
`.env` loader disabled.

**Checked and OK (verified):**
- `.env` is git-ignored and never committed; no key strings are in tracked files.
- Diagnostics never echo keys.
- HTTPS everywhere, with no `verify=False`, and a URL allowlist.
- Paid calls need enablement, a key, confirmation, an access code and a budget.
- Debug is off; the 500 page shows no traceback.
- Uploads: 413 above 10 MB, more than 20 files refused, a gzip bomb refused in 0.07 s.
- `/examples` and `/results` traversal variants all return 404; uploads are stored as `file_<n>`.
- XSS payloads are escaped.
- Tokens come from `secrets` (128-bit).
- The analysis slot is released on errors.
- 20 malformed-file cases returned 200 with a message, and no 500s.
- The container runs as non-root.

**Issues:**
- **P1-1 Memory exhaustion at review** (§11).
- **P1-2 Formula injection.** `report_exports.py:158-163, 213-216` writes strings beginning with `=`
  as formulas (14 formula cells verified in `formula_download.xlsx`). The CSV writers do not escape
  them either (`integrated_view.py:754`, `app.py:396`, `app.py:507`).
- **P2:**
  - The XLSX cell count can be bypassed: rows wider than the header, uncounted scan rows, and a large
    `header_row` from a profile (`excel_tables.py:146-176`).
  - Mapping-profile entries are not validated (`mapping_profiles.py:61-76`, `109-134`).
  - Count-based cache eviction (`app.py:220-231`).
  - Unbounded GeoNames calls (`geocode_locations.py:105-121`).
  - The paid-call gate relies on a loopback address when `AI_ACCESS_CODE` is unset; this is fragile
    behind a same-host proxy (`literature_view.py:99-105`). The budget ledger is in memory only.
- **P3:**
  - raw exception text in messages;
  - process-wide `redirect_stderr` in a threaded server (`app.py:734-745`);
  - staged dicts mutated outside the lock;
  - no CSRF tokens or security headers;
  - download URLs are bearer links not tied to a browser (10 min);
  - `LOCAL_MODEL_URL` allows http;
  - `defusedxml` is not installed;
  - `requests>=2.31` floor (raise to 2.32);
  - no lock file.

---

## 13. Test results

| Run | Result |
|---|---|
| Full suite with real MAFFT, `FONTENELE_MMC2`, `CRITS_CHRISTOPH_DIR`, `SCHUMANN_DIR`, AI off | **`Ran 674 tests in 652.331s — OK (skipped=5)`** |
| The 5 skips | 3 SciPy cross-checks (SciPy not installed), the NCBI Virus real-data test (files not present), and the mpox GenBank test (env var not set) |
| Optional mpox test run separately with real NC_063383.1 | **OK** (18 tests in the module) |
| Published-dataset script | 33 PASS, 0 FAILED |
| Audit-only diagnostics | 30 import scenarios; website export comparison; Render-like run; Berlin35 memory; ~181,000 annotation checks; 11 independent math scripts; 8 security scripts; 16 pathogen scripts |

**Gaps in the tests.**
- The suite does not cover any P1 issue.
- It also does not cover:
  - the empty-CDS_ID path;
  - "NA" text in the XLSX;
  - an undeclared genetic code;
  - partial dates in separate metadata;
  - segment collapse;
  - profile validation.
- Existing tests were not modified.

**Command:**
```bash
FONTENELE_MMC2=<mmc2.xlsx> CRITS_CHRISTOPH_DIR=<folder> SCHUMANN_DIR=<folder> AI_ENRICHMENT_ENABLED=false PATH=/c/Tools/mafft-win:$PATH .venv/Scripts/python.exe -m unittest discover -s tests
```

---

## 14. Prioritized issue register

| ID | Priority | Issue | Evidence |
|---|---|---|---|
| — | **P0** | None found | — |
| P1-1 | **P1** | Review-step parsing has no cell or row limit and runs outside the analysis slot. A 0.2 MB gz upload reached 980 MB, enough to kill a 512 MB server and every user's session | `security/memory_output.txt` |
| P1-2 | **P1** | Formula injection in XLSX (live formulas) and CSV downloads | `security/formula_output.txt` |
| P1-3 | **P1** | Review says ready, but the run fails on checks done only at run time (partial metadata dates, alt > depth, ALT = N, indel- or header-only iVar). The uploads are then deleted. This blocks Crits-Christoph with the authors' own metadata | `cc_ids/web_ww_meta.py`, `web/ready_then_fail.py`, `imports/out.txt` 05/21/27 |
| P1-4 | **P1** | A segment/chromosome column is not recognized; a stated reference merges all segments silently into one | `pathogen/test_web.py` W3 |
| P2-1 | P2 | Genome-position gene track built from min–max source labels ("UTR" band over ORF3a–N) | `exports/mmc2_web/genome_positions.png` |
| P2-2 | P2 | XLSX rebuilt with default NA parsing: "NA" gene or ID text becomes blank | `pathogen/test_na_export.py` |
| P2-3 | P2 | Missing `transl_table` silently uses table 1; alternative start codons are reported as "Start lost" | `pathogen/test_code4_notga.py` |
| P2-4 | P2 | A CDS without an ID loses the same-codon flag and its summary row | `math/t10_output.txt` |
| P2-5 | P2 | Mapping profiles: entries not validated; `header_row: -5` accepted silently | `security/fuzz_followup_output.txt` |
| P2-6 | P2 | XLSX cell-limit bypasses; 4 M-cell default too high for 512 MB | `security/memory_output.txt` |
| P2-7 | P2 | Count-based cache eviction by other visitors | `security/concurrency_output.txt` |
| P2-8 | P2 | GeoNames lookups unbounded while holding the slot | code (Inferred) |
| P2-9 | P2 | One-sequence reference; no mapping from iVar REGION labels to annotation seqids (blocks segmented and multi-chromosome pathogens) | `pathogen/out/test_h5n1_*.txt` |
| P2-10 | P2 | The Phase 4A report's Crits-Christoph metadata statement is incomplete (authors' `data/wastewater/sample_metadata.txt` matches) | §5 |
| P2-11 | P2 | SARS-CoV-2 study QC preset applied by default to any pathogen | quality_control.py:29-36 (Inferred) |
| P3-1 … P3-n | P3 | Percent ≤1 undetectable; "Coverage" breadth auto-mapped as depth; duplicate-conflict first-row use; frequency/count disagreement not raised at review; base-count file not used for coverage; Ref `N` in base counts; two VAF denominators; slope rule inconsistency; float32 Yue–Clayton; PCoA share label; ρ = 0 under-dispersion; Measurement_Class label; "CSV keeps every row" note; duplicate readiness sections; page length/size; raw exception text; stderr redirect; CSRF/headers; bearer download links; `requests` floor; defusedxml; GFF3 URL decoding; empty "not annotated: ." reason; GTF message | §4, §6, §9, §10, §12 |

---

## 15. Exact files and functions affected

| Issue | File:function (line) |
|---|---|
| P1-1 | `genomics/file_formats.py` `read_table` (198-200, `pd.read_csv(dtype=str)`); `app.py` `stage_uploads` (1300-1353), `review_upload` (1464), `run_reviewed` (1771, 1788) |
| P1-2 | `genomics/report_exports.py` `excel_rows` (158-163), `build_workbook` (213-216); `integrated_view.py` `download_files` (754); `app.py` (396, 507) |
| P1-3 | `wastewater/processed_snv.py` `load_processed_snvs` (199, 222-328, 286); `genomics/input_inspection.py` `review_dataset` / `prepare_inputs` (1287); `app.py` `run_reviewed` (1801-1804) |
| P1-4 | `genomics/column_detection.py` reference aliases (53-54); `genomics/import_options.py` `interpret_table` (declared reference, 184-188) |
| P2-1 | `genomics/integrated_analysis.py` `annotation_ranges` (153-160); `genomics/figures.py` `genome_position_figure` (146-215) |
| P2-2 | `genomics/report_exports.py` `StoredTables.get` (146) |
| P2-3 | `genomics/genome_annotation.py` (166, 299); `genomics/snv_consequences.py` (80) |
| P2-4 | `genomics/snv_consequences.py` `linkage_flags` (231), `summarize_variants` (310) |
| P2-5 | `genomics/mapping_profiles.py` `read_profile` (61-76), `apply_profile` (109-134); `app.py` `use_profile` (1420-1448) |
| P2-6 | `genomics/excel_tables.py` (28-31, 146-176); `Dockerfile` (no `MAX_WORKSHEET_CELLS`) |
| P2-7 | `app.py` (220-231, `RESULT_DOWNLOAD_LIMIT`, `STAGED_LIMIT`); `literature_view.py` (31-47) |
| P2-8 | `geocode_locations.py` (105-121) |
| Berlin35 memory | `genomics/pooled_observations.py` (observation matrix construction) |

---

## 16. Recommended repairs

Each repair is small, testable, and leaves calculations unchanged.

1. **P1-1.**
   - Add a cell budget for text tables, using the same idea as `MAX_WORKSHEET_CELLS`; for example
     `nrows`/column-count checks before `read_csv`.
   - Lower `MAX_WORKSHEET_CELLS` for the 512 MB image.
   - Run inspection inside the analysis slot or behind its own `BoundedSemaphore(1)`.
2. **P1-2.**
   - In the workbook, write strings beginning with `= + - @ \t \r` as text (`cell.data_type = "s"`, or
     a `'` prefix).
   - In CSV downloads, prefix them with `'`, and document the prefix in the README.
   - Add a test with a `=HYPERLINK` sample ID.
3. **P1-3.**
   - Move the `load_processed_snvs` validations (alt ≤ depth, single-base ALT, non-empty table, metadata
     date parsing) into `review_dataset`, so Step 2 blocks Run with the file name and row.
   - Treat unparseable or partial metadata dates as "unknown", as variant-table dates already are,
     instead of failing.
   - Keep staged files when a run fails.
   - Add a regression test with month-day dates.
4. **P1-4.** Recognize `segment`, `seqid`, `chromosome`, `contig` and `region` as reference columns.
   When a stated reference would cover more than one distinct value of such a column, refuse it.
5. **P2-2.** Use `read_csv(..., keep_default_na=False, na_values=[""])` in `StoredTables.get`, then
   re-check XLSX equality.
6. **P2-1.**
   - Draw the gene track from the uploaded GenBank/GFF3 CDS when present.
   - Otherwise draw contiguous runs of each source label, not min–max.
   - Label it "from the source table's annotation column".
7. **P2-3.**
   - Warn when `transl_table` is missing ("assumed standard code 1").
   - Report an alternative start codon as Met at codon 1 (e.g. "start retained").
8. **P2-4.** Fall back to `Gene` + coordinates as the CDS key when `protein_id`/ID is missing.
9. **P2-5.** Validate profile entries with the same rules as `read_file_options`: types, role
   allowlist, `header_row` from 1 to 50.
10. **P2-7, P2-8.** Evict cached results by age only, or raise the limit and scope results to a
    browser cookie. Cap GeoNames lookups per run.
11. **Berlin35 memory.**
    - Store the observation matrix's text columns as categoricals: 533 → 50 MB, identical CSV bytes.
    - Then re-measure under 512 MB before raising `MAX_OBSERVATION_CELLS`.
12. **Documentation.** Correct the Crits-Christoph metadata statement in the Phase 4A report and the
    README (P2-10), and state the PCoA share basis and the reconstructed-count basis in the main
    findings.

---

## 17. Remaining scientific validation requirements

1. **Experimental grouping:**
   - validation against sequencing truth, e.g. clinical lineage prevalence or spike-in mixtures, on
     real wastewater;
   - meeting the pre-registered criterion before any default use.
2. **Overdispersion:** a real dataset with documented technical replicates (separately prepared
   libraries). None of the three published datasets has one, so empirical calibration has never run
   on real data.
3. **Yue–Clayton:** confirm the per-position aggregation (sum versus mean, treatment of uncovered
   positions) against the 2021 study's code, and use float64.
4. **Coverage-aware non-detection:** validate on a real dataset with per-position A/C/G/T depth across
   the genome. Crits-Christoph's inStrain covT files are a candidate.
5. **Other pathogens:**
   - real processed observations for at least one DNA virus (e.g. mpox) and one bacterium;
   - multi-segment runs after P1-4 and P2-9 are fixed.
6. **Reconstructed counts:** state that intervals from reconstructed counts (Fontenele) depend on the
   source's rounding precision.

---

## 18. Deployment readiness

| Question | Answer |
|---|---|
| Correct for the supported SARS-CoV-2 processed-table workflows? | **Yes** (verified on three published datasets and through the website) |
| Fits Render Free for the normal datasets? | **Likely** in the local imitation (`mmc2` 301 MB / 120.6 s; reports on demand). **Not tested on Render** |
| Safe to expose publicly? | **Not yet**: P1-1 (memory exhaustion from a tiny upload) and P1-2 (formula injection in downloads) |
| Honest about limits? | Mostly yes (missing ≠ 0%, experimental labels, reference status). Exceptions: P1-3, P1-4, P2-1, P2-2 |

## 19. Recommended next development steps

1. Fix P1-1 to P1-4 and P2-2 as five separate small changes, each with a regression test.
2. Re-run the full suite and `validation/check_published_datasets.py`. Add the authors' Crits-Christoph
   `sample_metadata.txt` case: dates should become unknown, and the run should complete.
3. Deploy to a **non-production** Render service with explicit approval. Repeat the `mmc2` run and the
   on-demand XLSX, watch memory, and try two concurrent users.
4. Run the human usability checklist (§9).
5. Then the P2 items: gene track, genetic-code warning, CDS key, profile validation, cache eviction,
   and the Berlin35 categorical optimization.
6. Validation work from §17, prioritizing a real replicate dataset and a real non-SARS-CoV-2 observation
   table.

---

## 20. Addendum: P1 repairs made after the audit (2026-11-15)

The user approved repairing the four P1 issues. Sections 1–19 describe the code **as audited**; this
section records what changed.
- No scientific calculation, threshold or validation rule was changed.
- Nothing was committed, pushed or deployed.
- Regression tests are in `tests/test_audit_p1_repairs.py` (11 tests).

| Issue | Repair | Re-verified |
|---|---|---|
| P1-1 | `genomics/file_formats.check_table_size`: CSV/TSV refused before parsing above `MAX_WORKSHEET_CELLS` cells or `MAX_TABLE_COLUMNS` (20,000) columns. Worksheet rows counted with their header and scan rows (`excel_tables.py`). Inspection, review and preparation run one at a time (`input_inspection.one_at_a_time`, 60 s wait, then "server busy"). Dockerfile `MAX_WORKSHEET_CELLS=1500000` | Audit probes: 0.2 MB gz 980 → **259 MB**; wide CSV 1,158 → **201 MB**; 4.2 M-cell XLSX refused, peak **252 MB** at the Docker limit. Fontenele, Schumann (mmc1, mmc2) and Crits-Christoph files still inspect normally at 1.5 M |
| P1-2 | `genomics/spreadsheet_safety.py`: text starting with `= @ \t \r`, or `+`/`-` when not a number, gets a leading apostrophe in every CSV and XLSX download (integrated, Sprint 2 relationship files, legacy pages, literature CSVs, CLI). Numbers, `-0.25`, and a lone `+`/`-` strand are unchanged | Security probe: XLSX formula cells 14 → **0**; CSV files with an unescaped payload 16 → **0** |
| P1-3 | `input_inspection.pooled_validation_review` runs the analysis's own `load_processed_snvs` on the prepared files during review. Anything Run would refuse is a Step 2 problem with Run disabled. A run that produces nothing returns to Step 2 and keeps the uploads. The metadata-date message now names an example value | Import scenarios 05/21/27 now stop at review; the other 27 scenarios are unchanged. Crits-Christoph S4 + authors' `sample_metadata.txt`: Run disabled, problem names `'5-19'` and says a missing year is never guessed |
| P1-4 | `segment`, `chromosome`, `chr`, `seqid`, `seq_id`, `seqname` added as reference-column aliases. A stated reference is refused when another column names several sequences | Pathogen W3 table: SEGA.1 verified 12/12; SEGB.1 and SEGC.1 kept apart as unverified (no merge, no conflicting duplicates) |

**Verification after the repairs:**
- **Full suite:** `Ran 685 tests in 598.0 s — OK (skipped=5)`. That is the 674 earlier tests plus 11 new;
  no existing test was modified, and the same five optional skips remain.
- **Published-dataset validation:** 33 PASS, 0 FAILED.
- **Render-like `mmc2`** (512 MB, 10% of one core): 123.5 s, peak 302 MB (before: 120.6 s, 301 MB). On-demand
  XLSX 39.8 s, PDF 5.4 s.

**Still open:**
- All P2 and P3 items. P2-6 is partly addressed: the worksheet count and the Docker limit are fixed,
  but profile `header_row` validation is not.
- Hosted Render testing.
- A human usability test.

---

**READY FOR HOSTED VALIDATION**

**At audit time** (sections 1–19) the verdict was **READY AFTER SPECIFIC REPAIRS**: four P1 problems
blocked a public hosted test. Those specific repairs are now made and re-verified (section 20):
- 685 tests pass;
- 33 of 33 published-data checks pass;
- the audit's own reproductions show bounded review memory, no live formulas, run-time refusals shown
  at review with the uploads kept, and segments kept apart.

The remaining P2 and P3 items are not blockers for a **non-production** hosted validation, but they
should be fixed before wider use. This verdict is about readiness for hosted *testing*. Hosting itself
and real users have not been tested.

The audit-time assessment follows.

The scientific core is sound for the supported SARS-CoV-2 processed-table workflows:
- 674 tests pass;
- 33 of 33 published-data checks pass;
- the statistics match independent implementations;
- the annotation agrees 5,613 of 5,613 with the Fontenele authors' labels and 0 mismatches in about
  181,000 checks;
- the exports agree value for value.

Four specific problems block a public hosted validation:
1. a 0.2 MB upload can push the review step to 980 MB on a 512 MB server;
2. user text becomes live spreadsheet formulas;
3. valid-looking reviews that fail at run time and delete the uploads, which affects a real published
   dataset;
4. silent merging of segments under one stated reference.

Each has a contained repair (§16) that changes no scientific calculation. Render itself and real users
have not been tested, so this audit does not claim hosted or usability validation.
