# Sprint 5 Phase 4 — Progress Record

Reference-based SNV annotation: coding regions, codons, and predicted amino-acid consequences.

- Branch `sprint4-integration`, HEAD `b641533`.
- All Sprint 5 work (Phases 1–4) is local and uncommitted. Nothing is committed, pushed, merged,
  or deployed.
- No paid API calls are made: every run sets `AI_ENRICHMENT_ENABLED=false`.
- Attempt counter: the goal allows 16 attempts. Each milestone records its attempt number.

---

## Milestone 0 — Baseline and architecture (attempt 1) — Done

### 0.1 Reports read

- `docs/SPRINT3_AUDIT_REPORT.md`
- `docs/SPRINT5_PHASE1_IMPLEMENTATION_REPORT.md` and `docs/SPRINT5_PHASE1_PROGRESS.md`
- `docs/SPRINT5_PHASE2_IMPLEMENTATION_REPORT.md` and `docs/SPRINT5_PHASE2_PROGRESS.md`
- `docs/SPRINT5_PHASE3_IMPLEMENTATION_REPORT.md` and `docs/SPRINT5_PHASE3_PROGRESS.md`

The Phase 3 report's §19 ("Readiness for biological mutation annotation") says profiles keep
reference, position, and alleles per `Variant_Key`, so an annotation step can attach consequences
per `Variant_Key` without changing the grouping. Phase 4 follows that plan.

### 0.2 Repository state (verified 2026-11-15)

- `git status`: 55 modified or untracked entries, all Sprint 5 Phase 1–3 work. Nothing is staged
  and nothing was discarded.
- **Full suite**, project `.venv`, MAFFT on PATH (`C:\Tools\mafft-win`), `FONTENELE_MMC2` set to
  the published workbook (md5 `c0a75b7fa5bde8db09746613354a9f64`, unchanged):
  **`Ran 585 tests in 420.6 s — OK (skipped=4)`**. The 4 skips are the 3 optional SciPy
  cross-checks and the optional NCBI files test, the same as the Phase 3 report.
- The reported Phase 3 baseline therefore holds:
  - 585 tests with 4 optional skips;
  - the real `mmc2.xlsx` tests run and pass;
  - experimental grouping is off by default (`grouping_mode = "profiles"` in `DEFAULT_SETTINGS`).
- Memory (297–300 MB peak on `mmc2.xlsx`) is taken from the Phase 3 report and re-measured in
  Milestone 9. Hosted validation has still not been done.

### 0.3 Architecture relevant to annotation (verified by reading the code)

- **Reference FASTA.** `genomics/references.py`:
  - `load_reference_fasta` accepts exactly one A/C/G/T/N sequence and stores its ID,
    accession, version, length, and SHA-256 checksum;
  - `references_compatible` never trusts labels alone;
  - `reference_base(reference, position)` uses 1-based positions (`sequence[position - 1]`).
- **Variant identifiers.** `Variant_Key = "{Reference_ID}:{Position}:{Ref}>{Alt}"` in every
  Sprint 3 and Sprint 5 table, with 1-based reference positions.
- **Pooled pipeline.** `run_integrated_analysis` loads SNVs and runs QC with
  `reference_verification` (the "Reference_Status" column), then the observation matrix, temporal
  trends, and `run_depth_aware_analysis`. That last step builds:
  - `sprint5_vaf_observations.csv`, one row per observation (7,973 rows on `mmc2`);
  - `sprint5_temporal_summary.csv`;
  - the Phase 2 uncertainty tables;
  - the Phase 3 profile and grouping outputs.
- **Existing annotation fields.** Only what an input supplies is used today:
  - `Published_*` columns from `mmc2.xlsx`: `REF_CODON`, `REF_AA`, `ALT_CODON`, `ALT_AA`,
    `Sy/Ns`, `Annotation`;
  - `Source_*` columns from iVar (`GFF_FEATURE`, codons, amino acids);
  - `Annotated_Region` in the evidence table.

  `annotation_ranges` in `integrated_analysis.py` draws bands that span only annotated records;
  its caption says they are not gene boundaries. Nothing in the code translates codons.
- **Uploads.** In the guided flow, `inspect_file` detects the file kind:
  - FASTA, table, Excel, VCF, FASTQ, and alignment are recognized;
  - other kinds are refused as unsupported, so GFF3 and GenBank files are refused today.

  `suggest_roles` and `review_dataset` assign roles, and `prepare_inputs` writes the files.
  `run_genomic_analyses` passes the paths to `run_integrated_analysis`.
- **Web view.** `integrated_view.build_dashboard` builds section 9 (Sprint 5). Phase 2 and Phase 3
  are rendered as their own panels (`uncertainty_panel`, `grouping_panel`) in
  `templates/index.html`.

### 0.4 Plan

1. **New `genomics/genome_annotation.py`.** Parse GenBank (Biopython, already a dependency) or
   GFF3 (a small parser), then:
   - match the annotation to the reference by exact accession.version and identical sequence;
   - verify each CDS (structure checks and translation);
   - mark complex structures unsupported.
2. **New `genomics/snv_consequences.py`.** Annotate each distinct SNV once:
   - reference identity and the reference base;
   - overlapping CDS, codon, strand, and translation;
   - a consequence category;
   - flags for codons shared with other observed SNVs.
3. **New `genomics/mutation_annotation.py`.** Join the annotations to the existing Sprint 5
   evidence and temporal tables by `Variant_Key`, without recomputing VAF. Add tables, one
   figure, and a view summary.
4. **An optional "Genome annotation (GenBank or GFF3)" upload role** in the existing guided flow,
   plus `--annotation` for the command-line runner.
5. **A website panel "Mutation Annotation and Protein Consequences"** in section 9, independent of
   grouping.

---

## Milestone 1 — Verified reference annotation model (attempt 1) — Done

- Official NCBI files for NC_045512.2 are in `examples/annotation/`:
  - `NC_045512.2.gb`: GenBank `gbwithparts`, from NCBI E-utilities `efetch`, 78,517 bytes, md5
    `84cf4fb7d7d061293cafb4e0847106e3`, last updated 18-JUL-2020;
  - `NC_045512.2.gff3`: from the NCBI sviewer GFF3 report, 14,704 bytes, md5
    `2d3c5085fe22f2f9a672ba53964f6ca3`;
  - `NC_045512.2.fasta`: the reference FASTA already used in Phase 1 checks, md5
    `863ee5dba1da0ca3f87783782284d489`.

  All three were downloaded on 2026-11-15 (the FASTA on 2026-11-14). They are public RefSeq
  records, about 120 KB together.
- `genomics/genome_annotation.py` (new). Coordinates are 1-based and inclusive, as in the files.
  Biopython parts are converted with `(start + 1, end)`.
- First load of the real files (verified):
  - **GenBank:**
    - 12 CDS, all supported;
    - for every CDS, translating the reference coding sequence reproduces the record's
      `/translation` exactly. This includes ORF1ab, `join(266..13468,13468..21555)` with
      `/ribosomal_slippage`: 7,096 amino acids;
    - the record's sequence is identical (SHA-256) to the FASTA.
  - **GFF3:** the same 12 CDS with phases 0,0. They pass the internal checks: whole codons, one
    terminal stop, no internal stop, and a start codon.
  - **Mature peptides:** their codon ranges in ORF1ab reproduce NCBI's own protein coordinates.
    For example, nsp2 is 181–818, nsp3 819–2763, and the RNA-dependent RNA polymerase (nsp12)
    4393–5324 across the frameshift.
  - **Fix made:** GFF3 writes the frameshifted nsp12 on two lines. Peptide lines are now grouped
    by ID, as CDS lines are.
- **Other real records** used for verification:
  - **Influenza A/Puerto Rico/8/1934 (H1N1).** RefSeq NC_002016.1–NC_002023.1, 8 segments,
    stored as `examples/annotation/influenza_A_PR8_segments.gb` (48,825 bytes, md5
    `1d2b03c3c0a8b5354663ff9cf4279c39`). All 12 CDS are supported, including spliced M2
    (`join(26..51,740..1007)`), spliced NEP, and the PA-X +1 frameshift
    (`join(25..597,599..784)`).
  - **Mpox.** NC_063383.1 has 179 CDS, 96 of them on the minus strand. All reproduce the
    record's `/translation`. The file (489 KB) is kept outside the repository; an optional test
    reads it from `ANNOTATION_MPOX_GENBANK`.
- Status: **Done** (attempt 1).

---

## Milestone 2 — SNV consequence annotation (attempt 1) — Done

- `genomics/snv_consequences.py` (new).
- **Each distinct SNV is annotated once**:
  - `drop_duplicates("Variant_Key")`, with no per-sample translation;
  - repeated records are tested: 800 records give 2 SNVs and 2 rows.
- **Steps for each SNV:**
  1. Allele check: single A/C/G/T bases that differ.
  2. Reference identity: the exact accession.version, with sequence compatibility.
  3. Reference-base check.
  4. Every overlapping CDS gives one row.
  5. The codon is read on the coding strand; Alt is complemented for minus-strand genes.
  6. Translation uses the CDS's table.
  7. The change is classified with Sequence Ontology categories.
- **Unsupported or ambiguous cases**, each reported with its reason:
  - indels and multi-base changes;
  - ambiguous alleles or reference bases;
  - codons that contain N;
  - unsupported CDS structures;
  - a base read twice in a frameshifted protein (SARS-CoV-2 13468 in pp1ab).
- **Not annotated** (never transferred), with the reason:
  - a different accession or version;
  - a reference mismatch;
  - a position outside the reference;
  - no sequence available.
- **Shared codons.** SNVs in the same codon as another SNV *reported in the same sample* are
  flagged. The combined amino-acid change is never predicted.
- **Mature peptides.** The peptide and its position are given where an in-frame mature peptide
  exists, for example nsp12 P323L for ORF1ab P4715L.

## Milestone 3 — Biological correctness (attempt 1) — Done

- **Synthetic, hand-derived checks.** `tests/annotation_fixtures.py` builds a 2-segment genome
  with these features:
  - a plus-strand CDS;
  - a minus-strand CDS;
  - a +1-frame overlapping CDS;
  - a plus-strand spliced CDS with a codon that spans the intron;
  - a minus-strand spliced CDS;
  - an internal-stop CDS;
  - an RNA-editing CDS;
  - UTR and N bases.

  Every expected codon and amino acid was derived by hand. 19 cases are checked in both GFF3
  and GenBank.
- **Real SARS-CoV-2.** Codons were sliced from the FASTA with plain arithmetic, independently of
  the code:
  - S D614G (GAT→GGT);
  - ORF1ab P4715L, which is nsp12 P323L, past the frameshift;
  - C3037T F924F, which is nsp3 F106F, in both ORF1a and ORF1ab;
  - C241T as 5'UTR noncoding;
  - N R203K, R203R, and G204R, with a shared-codon flag;
  - ORF8 Q27*;
  - N stop lost (*420Q);
  - S start lost;
  - ORF7a/ORF7b overlap: ORF7b start lost while ORF7a changes;
  - base 13468 unsupported in pp1ab and supported in pp1a.
- **Independent published annotation (real).** The `mmc2.xlsx` columns `REF_CODON`, `REF_AA`,
  `ALT_CODON`, `ALT_AA`, `Sy/Ns`, and `Annotation` come from the authors' own pipeline. Result:
  - **5,613 of 5,613 coding SNVs agree exactly**: gene, reference and alternative codon, and
    both amino acids;
  - synonymous and non-synonymous agree 100% (stop retained = Syn; stop gained, stop lost, and
    start lost = Nsyn);
  - the 64 UTR SNVs are labeled "Syn" in the workbook. This application calls them noncoding,
    as the goal requires: noncoding substitutions are not synonymous.

  No expected value was changed to make a test pass; no discrepancy was found.
- **Mature-peptide spans** (ORF1ab, 16 peptides) equal NCBI's own protein coordinates, which
  are encoded in the GFF3 IDs.

## Milestone 4 — Connection to Sprint 5 results (attempt 1) — Done

- `genomics/mutation_annotation.py` (new) joins by `Variant_Key` and reads the existing Phase 1–3
  tables only:
  - `sprint5_vaf_observations.csv`, `sprint5_temporal_summary.csv`, and
    `sprint5_mutation_groups.csv` (experimental only);
  - individual variants.
- **New tables:**
  - `sprint5_mutation_annotations.csv`: one row per SNV per CDS;
  - `sprint5_annotated_mutations.csv`: one row per SNV, with the consequence plus evidence
    (samples, dates, median and maximum VAF, depth range, count basis, Wilson interval at the
    maximum VAF, calibrated beta-binomial count, reference status, QC exclusions, and series with
    2+ dates);
  - `sprint5_annotated_observations.csv`: one row per pooled observation, with the original VAF,
    depth, counts, intervals, and QC copied unchanged;
  - `sprint5_annotation_sources.csv`;
  - `sprint5_annotation_comparison.csv`.
- **A web test checks** that VAF, depth, Wilson bounds, count source, reference status, and QC in
  the annotated observations equal `sprint5_vaf_observations.csv` column for column.

## Milestone 5 — Website section (attempt 1) — Done

- **Section 9 has a new panel, "Mutation Annotation and Protein Consequences".** It shows:
  - the source line;
  - counts: evaluated, annotated, synonymous, missense, other coding, noncoding, unverified or
    unannotatable, overlap disagreements, and shared codons;
  - one figure: positions by category over the CDS map, plus counts per CDS;
  - the main table, 25 rows plus 75 expandable, with plain column headings;
  - collapsible technical details: codons, coding positions, and CDS checks;
  - the limitations and the downloads.
- **Upload.** An optional role, "Genome annotation (GenBank or GFF3)", is added in the existing
  guided upload. The file kind is detected from the content: `LOCUS` or `##gff-version 3`. The
  review step reports whether the annotation matches the variants' reference. The NCBI GenBank
  file can be downloaded as the example.
- **Without an annotation**, the panel explains that coding consequences cannot be determined, and
  every other analysis runs.

## Milestones 6–7 — Safeguards and Phase 3 status (attempt 1) — Done

- The text on the page and in the limitations separates predicted consequences from biological
  effects. No transmissibility, severity, escape, resistance, replication, or outbreak claim is
  made. No AI is involved.
- **Grouping default unchanged:** `DEFAULT_SETTINGS["grouping_mode"] == "profiles"` (tested).
  - Group tables return 404 in default runs with an annotation.
  - Only when the user chooses the experimental view does the annotated table gain
    `Exploratory_Group`, with every value ending "(experimental, not validated)".
  - The limitations say that annotation does not validate groups.
  - No grouping code was changed.

## Milestone 8 — Tests and end-to-end validation (attempt 1) — Done

- **New tests:**
  - `tests/test_genome_annotation.py`: 18 tests, 1 optional (mpox);
  - `tests/test_snv_consequences.py`: 22;
  - `tests/test_sprint5_phase4_web.py`: 10 Flask upload-to-download tests, 1 needing
    `FONTENELE_MMC2`. They cover:
    - count-based pooled data with GenBank;
    - missing annotation;
    - GFF3 with and without a matching FASTA;
    - another pathogen's annotation;
    - a wrong reference base;
    - segmented influenza;
    - combined Sprint 2 + 3 with a test GFF3;
    - experimental grouping with annotation;
    - the unchanged default;
    - real `mmc2.xlsx` with NCBI GenBank (sparse, reconstructed counts).
- **Output regression** against the saved Phase 3 outputs (with MAFFT; no annotation supplied):
  - `--example`: 62 of 64 identical;
  - `mmc2`: 40 of 42 identical, including every PNG.

  The only differences are one new `analysis_status.csv` row ("Skipped: no annotation") and one
  new `provenance.csv` row. With the GenBank file supplied, `mmc2` keeps the same 40 identical
  files and adds 6 Phase 4 files.
- **Refactor for speed.** The per-group pandas operations were replaced by plain dictionary
  passes (see Milestone 9). The output tables are identical before and after; this was checked
  on `mmc2`.

## Milestone 9 — Performance and Render readiness (attempt 1) — Done

Measured with the guided web path `/review` → `/run` on `mmc2.xlsx` (Windows peak working set,
`OPENBLAS_NUM_THREADS=1`):

| Run | No annotation | With NCBI GenBank |
|---|---|---|
| One run | 295–296 MB, 12.0–12.2 s | 308 MB, 12.8–13.0 s (annotation step 0.67–0.88 s) |
| Four runs in one process | 323 MB | 335 MB |
| 512 MB cap, 10% of one CPU core | 288 MB, 124.6 s | 302 MB, 131.9 s |
| HTML page | 3.01 MB | 3.17 MB |

- **First attempt:** the annotation step took 6.8 s and peaked at 325 MB, because of per-group
  pandas operations. They were replaced by dictionary passes; the tables are identical.
- **Downloads added:** 2.2, 2.0, 1.8, and 0.6 MB CSVs, plus an 80 KB PNG. Every file is under
  the 8 MB per-file cache limit.
- No dense genome-wide matrix is created. Each distinct SNV is translated once.

## Milestone 10 — Final report (attempt 1) — Done

- `docs/SPRINT5_PHASE4_IMPLEMENTATION_REPORT.md`, 14 sections plus the completion table.
- **Final suite:** `Ran 636 tests in 445.9 s — OK (skipped=5)`.
- **Attempts used:** 1 of 16 per milestone. One fix during Milestone 1 (grouping GFF3 peptide
  lines) and one performance fix during Milestone 9.
- The in-app browser refused the local preview, so the page was verified through the Flask routes
  and the figure file instead.

**SPRINT 5 PHASE 4 COMPLETE**
