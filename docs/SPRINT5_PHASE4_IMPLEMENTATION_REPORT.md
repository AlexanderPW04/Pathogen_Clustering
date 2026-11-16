# Sprint 5 Phase 4 — Implementation Report

**Reference-based SNV annotation: genes, codons, and predicted amino-acid consequences.**

- Branch `sprint4-integration`, HEAD `b641533`.
- All Sprint 5 work (Phases 1–4) is local and uncommitted. Nothing was committed, pushed, merged,
  or deployed.
- No paid API calls were made: `AI_ENRICHMENT_ENABLED=false` in every run, and no AI is used in
  this phase.

Labels used below:
- **verified**: run or read in this session;
- **real**: the published `mmc2.xlsx` or official NCBI RefSeq records;
- **synthetic**: invented test data with hand-derived expectations.

## Summary

The application can now take an optional genome annotation (GenBank or GFF3) of the same
reference used for variant calling. For each distinct observed SNV it reports:

- the gene or coding region;
- the codon and coding position;
- the reference and alternative amino acids;
- a Sequence Ontology consequence;

together with that SNV's existing sequencing evidence (VAF, depth, count basis, Wilson and
beta-binomial intervals, reference and QC status, and temporal series).

Main results:

- **Independent check on the real workbook (`mmc2.xlsx`).** The prediction agrees exactly with
  the codons and amino acids the authors published for **all 5,613 coding SNVs**. It reports
  their 64 UTR SNVs, which the workbook labels "Syn", as noncoding.
- **Hand-derived checks.** Every hand-derived codon test passes for:
  - plus- and minus-strand genes;
  - spliced and overlapping reading frames;
  - the SARS-CoV-2 −1 frameshift and influenza's +1 frameshift;
  - segmented genomes.
- **Honest failure modes.**
  - Unsupported coding structures are marked unsupported.
  - Incompatible references are never annotated.
  - Each pooled SNV is predicted alone: there is no haplotype claim.
- **Phase 3 grouping is unchanged.** It remains off by default and labelled experimental.
- **Existing outputs are preserved.** Phase 1–3 outputs are identical apart from one status row
  and one provenance row.

## 1. Methods implemented

Three new modules, each small enough to read in one sitting:

| Module | Role |
|---|---|
| `genomics/genome_annotation.py` | Loads GenBank (Biopython, already a dependency) or GFF3 (a small parser with `##FASTA` support) into one model, matches it to the reference, and verifies every CDS |
| `genomics/snv_consequences.py` | Annotates each distinct SNV once against the verified model and classifies the consequence; flags shared codons |
| `genomics/mutation_annotation.py` | Joins the annotation to the existing Sprint 5 tables by `Variant_Key`, writes the Phase 4 tables and one figure, compares with input annotation, and builds the panel summary |

Pipeline position: inside `run_integrated_analysis`, after the cross-source comparison:

1. Phase 1–3 tables already exist.
2. `run_mutation_annotation(...)` reads them.
3. It adds new tables only.

When no annotation file is supplied, it records "not supplied" and returns.

## 2. Reference and annotation requirements

- **Accepted files:**
  - an NCBI-style GenBank flat file (`LOCUS ...`), one record per sequence or segment;
  - GFF3 (`##gff-version 3`), optionally with a `##FASTA` section.

  The kind is detected from the file content, not the name.
- **Gene boundaries come only from the file.** Nothing is inferred from the nucleotide
  sequence.
- **Model fields** (`sprint5_annotation_sources.csv`):
  - `Reference_ID` (accession.version), `Accession`, and `Version`;
  - `Organism` (GenBank `ORGANISM`; for GFF3, an `organism` attribute or the NCBI taxon ID);
  - `Segment` (GenBank `/segment` or GFF3 `segment=`);
  - gene name and NCBI GeneID, locus tag, CDS identifier (`protein_id`), and product;
  - 1-based parts in transcription order, strand, phases (GFF3) or `codon_start` (GenBank), and
    the translation table;
  - structure (single, joined, or programmed ribosomal frameshift), mature peptides (codon
    ranges), status and reason, the verification performed, and the source file and format.
- **Reference compatibility.** Annotations are applied only when **all** of these hold:
  1. The variant's `Reference_ID` equals an annotated accession.version **exactly**. A
     different version, such as `NC_045512.1`, or another pathogen is never annotated.
  2. If a reference FASTA is uploaded with the same ID, the annotation's own sequence must be
     **identical** (SHA-256).
  3. A GFF3 file without a sequence borrows the FASTA only when the ID and the declared length
     match.
  4. The reference base at every SNV position must equal the variant's `Ref`.

  A SARS-CoV-2 annotation therefore cannot be applied to influenza data, or the reverse; this
  is tested both ways.
- **Optional.** The pipeline runs with or without the file. The review step states whether the
  annotation matches the variants' reference.

## 3. Exact consequence-classification rules

Each SNV is evaluated against the reference **alone**, in this order (the first failure stops
it):

| Step | Failing outcome |
|---|---|
| Ref and Alt are single, different A/C/G/T bases | *Unsupported or ambiguous* (indel, multi-base, IUPAC, or Ref = Alt) |
| The annotation describes the variant's accession.version and is compatible | *Not annotated* ("never transferred between different references") |
| A reference sequence is available | *Not annotated* (no sequence) |
| The position is inside the reference | *Not annotated* |
| The reference base is A/C/G/T | *Unsupported* (ambiguous reference base) |
| The reference base equals Ref | *Not annotated* ("Mismatch: reference has X") |

Then every CDS containing the position gives one row:

| Condition | Consequence (Sequence Ontology) |
|---|---|
| The CDS failed verification | *Unsupported*, with the CDS reason |
| The base appears twice in the coding sequence (ribosomal slippage overlap) | *Unsupported* ("read twice") |
| The codon contains an ambiguous reference base | *Unsupported* |
| Codon 1 and the start codon was verified | **Start lost** (SO:0002012: any change to the start codon) |
| Codon 1 and the start codon was not verified | *Unsupported* |
| Reference amino acid `*` and alternative `*` | **Stop retained** (SO:0001567) |
| Reference amino acid `*` and alternative a sense codon | **Stop lost** (SO:0001578) |
| Alternative amino acid `*` | **Stop gained** (SO:0001587) |
| Same amino acid | **Synonymous** (SO:0001819) |
| Different amino acid | **Missense** (SO:0001583) |
| The position is in no CDS | **Noncoding** (with UTR, stem-loop, or gene context); never synonymous |

**Per-SNV category** (`sprint5_annotated_mutations.csv`):
- if every overlapping CDS gives the same consequence, that consequence;
- otherwise "Differs between overlapping coding regions", with every per-CDS consequence kept in
  `Consequences_By_CDS` and in `sprint5_mutation_annotations.csv`. Nothing is chosen
  arbitrarily.

**Annotation status:**
- *Annotated*: every row is supported;
- *Partly annotated*: some overlapping CDS are unsupported;
- *Not annotated*.

**Amino-acid notation:** one-letter, for example `D614G`, `Q27*`, or `*420Q`, plus three-letter
columns. The coding change is `c.<coding position><ref>><alt>` on the coding strand.

**Linkage safeguard:**
- When another SNV at a different position changes the same codon of the same CDS **and both
  were reported in at least one common sample**, both are flagged (`Shared_Codon_With` and
  `Linkage_Note`).
- Each is still predicted alone; the combined codon is never translated.
- Without sample information (individual-only data), every pair in a shared codon is flagged.

## 4. Coding-frame and coordinate handling

- **Files and tables:** 1-based, inclusive (GenBank `266..13468`, GFF3 columns 4–5).
- **Biopython:** `FeatureLocation` is 0-based and end-exclusive; it is converted with
  `(int(part.start) + 1, int(part.end))`. Biopython's `join()` parts are kept in biological
  order. GFF3 parts are sorted by start, ascending for `+` and descending for `-`.
- **Sequence strings:** the base at 1-based position `p` is `sequence[p - 1]`.
- **Coding offsets** (`coding_offsets`): 0-based along the coding strand. For each part
  (start, end) in transcription order:

  `offset = before + (p − start)` for `+`, and `before + (end − p)` for `−`,

  where `before` is the summed length of the earlier parts.
  - The codon is `offset // 3 + 1` and the codon position is `offset % 3 + 1`.
  - Overlapping parts (a −1 frameshift) give two offsets for one base.
- **Coding sequence:** the concatenated parts, each reverse-complemented for `−` genes. The
  alternative base is complemented for `−` genes before it is placed in the codon.
- **Phase checks:**
  - **GFF3:** every part's phase must equal `(3 − before mod 3) mod 3`; otherwise the CDS is
    unsupported.
  - **GenBank:** `codon_start` must be 1. Any other first phase means a partial 5' end, which is
    unsupported.
- **CDS verification** (`verify_cds`). The CDS must have:
  - whole codons;
  - a terminal stop codon;
  - no internal stop codon.

  In addition:
  - **GenBank:** the translation of the reference coding sequence must equal the record's
    `/translation` exactly. An alternative start codon is read as Met.
  - **GFF3:** there is no protein to compare, so only the internal checks are recorded.
  - The start codon is "verified" when the first codon is a start codon of the CDS's
    translation table.

## 5. Files and functions modified

**New:**
- `genomics/genome_annotation.py`:
  - loading: `detect_annotation_format`, `load_genbank`, `load_gff3`;
  - checking: `attach_reference`, `coding_offsets`, `build_coding_sequence`, `translate_codon`,
    `structure_problem`, `verify_cds`;
  - mature peptides: `map_peptides`, `peptide_span`;
  - entry points and output: `finalize_annotation`, `load_annotation`, `source_table`.
- `genomics/snv_consequences.py`: `classify_codon_change`, `cds_row`, `annotate_snv`,
  `linkage_flags`, `combined_consequence`, `annotate_variants`, `summarize_variants`.
- `genomics/mutation_annotation.py`: `distinct_variants`, `joined_counts`, `evidence_summary`,
  `temporal_summary`, `input_annotation_comparison`, `annotated_mutation_table`,
  `run_mutation_annotation`.
- `examples/annotation/`: official NCBI files (real, public RefSeq):
  - `NC_045512.2.gb`;
  - `NC_045512.2.gff3`;
  - `NC_045512.2.fasta`;
  - `influenza_A_PR8_segments.gb`.
- Tests:
  - `tests/annotation_fixtures.py` (synthetic hand-built segments);
  - `tests/test_genome_annotation.py`;
  - `tests/test_snv_consequences.py`;
  - `tests/test_sprint5_phase4_web.py`.
- Docs: `docs/SPRINT5_PHASE4_PROGRESS.md` and this report.

**Modified (additions only):**
- `genomics/file_formats.py`: `ANNOTATION_KIND` and `ANNOTATION_EXTENSIONS`; content detection
  of `LOCUS` / `##gff-version`. Files of other kinds are classified exactly as before.
- `genomics/input_inspection.py`:
  - an `annotation` role, `inspect_annotation`, and `prepare_inputs` copying the file;
  - `annotation_status` and `variant_reference_ids` for the review's compatibility table;
  - a "Mutation annotation and protein consequences" row in the analysis plan.
- `genomics/integrated_analysis.py`: an `annotation_path` argument, the call to
  `run_mutation_annotation`, an `analysis_status` row, and `provenance` rows.
- `app.py`:
  - an optional upload field and named field (`annotation_file`);
  - `annotation_path` passed to the analysis;
  - the NCBI GenBank file offered as an example download.
- `run_integrated_analysis.py`: `--annotation`.
- `genomics/figures.py`: `mutation_annotation_figure` and `cds_lanes`.
- `integrated_view.py`: `annotation_view` and the panel's data.
- `templates/index.html`: the `annotation_panel` macro in section 9.
- `templates/_upload.html`: format and help text.
- `dashboard_guide.py`: the figure guide and a one-sentence interpretation.
- `README.md`: the format table, the requirements table, and real-data instructions.

**Not modified:**
- every Phase 1–3 calculation module: `variant_evidence`, `beta_binomial`, `overdispersion*`,
  `vaf_trends`, `mutation_profiles`, `mutation_grouping`, `mutation_group_analysis`;
- the QC, Sprint 2, and Sprint 4 code;
- `requirements.txt` and the Dockerfile. **No new dependency.**

## 6. Independent biological verification

All results are **verified** (they were run in this session).

1. **NCBI protein sequences (real).**
   - **SARS-CoV-2:** for all 12 CDS of NC_045512.2, the translation of the reference coding
     sequence equals the record's `/translation`. This includes pp1ab across the −1 frameshift
     (7,096 aa).
   - **Influenza A/PR/8/34:** all 12 CDS on 8 segments match. These include spliced M2
     `join(26..51,740..1007)`, spliced NEP, and PA-X with a +1 frameshift
     `join(25..597,599..784)`.
   - **Mpox (NC_063383.1):** all 179 CDS match, **96 of them on the minus strand**. This is an
     optional test; the 489 KB file is not in the repository.
2. **GenBank versus GFF3.** For NC_045512.2 both give identical parts and coding sequences for
   all 12 CDS.
3. **Mature-peptide coordinates.** The 16 ORF1ab peptide codon ranges computed from genomic
   coordinates equal NCBI's protein coordinates encoded in the GFF3 IDs (for example nsp2
   181..818). This includes nsp12, which crosses the frameshift (4393..5324).
4. **Hand-derived codons (synthetic).** 19 cases on hand-built segments, each checked in both
   GFF3 and GenBank. They cover:
   - plus- and minus-strand genes;
   - a minus-strand spliced gene;
   - a plus-strand spliced gene with a codon spanning the intron;
   - a +1-frame overlap;
   - all six coding categories.
5. **Hand arithmetic on the real SARS-CoV-2 FASTA.** Each codon was sliced at
   `start + 3 × (codon − 1)`, independently of the code:

   | SNV | Expected | Result |
   |---|---|---|
   | A23403G | S codon 614 GAT→GGT: D614G | Matches |
   | C14408T | pp1ab codon 4715 (offset 13203 + 940), P4715L = nsp12 P323L | Matches |
   | C3037T | F924F = nsp3 F106F in both ORF1a and ORF1ab | Matches |
   | C241T | 5'UTR noncoding | Matches |
   | G28881A, G28882A, G28883C | N R203K, R203R, G204R; 28881 and 28882 flagged as sharing codon 203 | Matches |
   | C27972T | ORF8 Q27* | Matches |
   | T29531C | N *420Q (stop lost) | Matches |
   | A21563G | S start lost | Matches |
   | A27756G | ORF7b start lost, while ORF7a codon 121 changes | Matches; "differs between overlapping regions" |
   | 13468 | Unsupported in pp1ab (read twice); supported in pp1a | Matches |

   These are widely reported names (for example D614G and nsp12 P323L). The tests check that the
   arithmetic reproduces them, rather than trusting the module.
6. **The authors' own annotation in `mmc2.xlsx` (real, independent pipeline):**
   - all 5,613 coding SNVs agree on gene, reference codon, alternative codon, and both amino
     acids;
   - the synonymous/nonsynonymous split agrees 100%;
   - the 64 UTR SNVs are "Syn" in the workbook and noncoding here. This is a known labeling
     difference, not an error in the prediction.

   No expected value was changed to make a test pass, and no discrepancy needed investigation.

## 7. Real workbook compatibility (`mmc2.xlsx`, md5 `c0a75b7fa5bde8db09746613354a9f64`, unchanged)

**Without an annotation:**
- the panel says coding consequences cannot be determined;
- `analysis_status.csv` has "Skipped" for Phase 4;
- every other output is identical to Phase 3.

**With `NC_045512.2.gb` (verified through the web path `/review` → `/run`):**
- **Reference.** The workbook's reference label (NC_045512.2, assigned from the paper's
  methods) matches the GenBank accession.version. All 5,677 distinct SNV reference bases agree
  with the RefSeq sequence. The pooled-QC reference status stays "Supplied but unverified",
  because no FASTA was uploaded; Phase 4 does not change that column.
- **Counts:**

  | Category | Distinct SNVs |
  |---|---|
  | Evaluated | 5,677 |
  | Annotated | 5,677 |
  | Missense | 4,109 |
  | Synonymous | 1,243 |
  | Stop gained | 243 |
  | Start lost | 12 |
  | Stop lost | 3 |
  | Stop retained | 2 |
  | Differs between overlapping regions (ORF7a/ORF7b) | 1 |
  | Noncoding | 64 |

- **Shared codons.** 348 SNVs share a codon with another SNV reported in the same sample, and
  are flagged.
- **Evidence.** All 7,973 observations are joined with their unchanged VAF and depth. Counts are
  "Reconstructed" because the workbook has only frequencies and depths, as in Phase 1.
- **Grouping.** No groups, no invented zeros, and no ρ, as in Phase 3.

## 8. Website integration

**Upload (step 1).** A new optional field, "Genome annotation (GenBank or GFF3, optional)":
- it is offered for individual, pooled, and combined modes;
- the official NCBI GenBank file is offered as the example;
- the automatic upload box recognizes the file by its content.

**Review (step 2).**
- The file is listed with its references, organism, segments, and CDS count.
- The compatibility table says whether it describes the variants' reference, and warns that
  otherwise nothing will be annotated.
- The analysis plan gains a "Mutation annotation and protein consequences" row.

**Results (section 9, "Mutation Annotation and Protein Consequences").** The panel sits before
the Phase 3 grouping panel and does not depend on it. It contains:
1. A plain-language explanation and a source line: the file, references, and how many coding
   regions passed the checks.
2. The required counts: distinct SNVs evaluated, successfully annotated, synonymous, missense,
   other coding consequences (with breakdown), noncoding, and unverified or unannotatable. Plus
   overlap disagreements and shared codons when present.
3. A boxed reading note: a predicted change is not a biological effect, and pooled reads do not
   show which mutations share a genome.
4. The independent-check sentence, when the input carried its own codons.
5. One figure, "Where the mutations fall and what they are predicted to change":
   - **top:** each SNV at its position, one row per category, over the CDS map, with overlapping
     reading frames on separate lanes and unsupported CDS hatched;
   - **bottom:** counts per CDS by consequence.

   About 80 KB PNG.
6. The main table "Observed mutations: predicted consequence and sequencing evidence", with the
   columns mutation, gene, mature peptide, nucleotide change, predicted amino-acid change,
   consequence, reference verification, pooled samples, median and highest VAF, Wilson interval
   at the highest VAF, depth range, read-count basis, site series with 2+ dates, individual
   carriers, and the exploratory group (only when the user enabled it). It shows 25 rows, 75
   more on request, and the full CSV.
7. Collapsed technical details: per-CDS codons, coding positions, mature peptides, linkage
   notes, and the CDS check table.
8. Collapsed limitations and the 5 downloads.

## 9. Full automated test results

**Command**, from the project folder:

```bash
python -m unittest discover -s tests -v
```

Environment: project `.venv`, MAFFT on PATH (`C:\Tools\mafft-win`), `FONTENELE_MMC2` pointing to
the published workbook, and `ANNOTATION_MPOX_GENBANK` unset (as a reviewer would usually run
it).

| Run | Result |
|---|---|
| Baseline at the start of Phase 4 | `Ran 585 tests in 420.6 s — OK (skipped=4)` |
| **Final** (final code) | **`Ran 636 tests in 445.9 s — OK (skipped=5)`**: 585 pre-existing and 51 new. The 5 skips are the 3 optional SciPy cross-checks (SciPy is not in `.venv`), the optional NCBI files test, and the optional mpox test. All real-`mmc2.xlsx` tests and the real-MAFFT test ran |

New tests:

| File | Tests | Notes |
|---|---|---|
| `test_genome_annotation.py` | 18 | 1 optional (mpox) |
| `test_snv_consequences.py` | 22 | |
| `test_sprint5_phase4_web.py` | 11 | real Flask routes; 1 needs `FONTENELE_MMC2` |

The optional mpox test was also run with the file set: `Ran 18 tests — OK`.

**Output regression** against the saved Phase 3 outputs (verified):

| Input | Identical | Differences |
|---|---|---|
| `--example`, with MAFFT, no annotation | 62 of 64 files, including every PNG | `analysis_status.csv` (+1 "Skipped" row) and `provenance.csv` (+1 row) |
| `mmc2`, no annotation | 40 of 42 files | the same two rows |
| `mmc2` with GenBank | the same 40 files | those two files plus 6 new Phase 4 files |

## 10. Runtime and peak memory

**Measured locally:** Windows peak working set, `.venv`, `OPENBLAS_NUM_THREADS=1`, the guided web
path `/review` → `/run`, and the same method as Phases 2–3. Render was **not** measured.

| Run (`mmc2.xlsx`) | Phase 3 report | Phase 4, no annotation | **Phase 4, with NCBI GenBank** |
|---|---|---|---|
| One run, peak | 297–300 MB | 295–296 MB | **308 MB** (twice) |
| One run, time | 11.7 s | 12.0–12.2 s | 12.8–13.0 s |
| Annotation step (load, annotate 5,677 SNVs, join, compare, figure) | — | 0 s | 0.67–0.88 s |
| Four runs in one process, peak | 318–325 MB | 323 MB | 335 MB |
| 512 MB cap and 10% of one CPU core | 288 MB, 118.8 s | 288 MB, 124.6 s (today's A/B) | **302 MB, 131.9 s** (job committed 304 MB; annotation step 6.6 s at 10% CPU) |
| HTML page | — | 3.01 MB | 3.17 MB |

**Annotation lookup cost:**
- 5,677 distinct SNVs are annotated once each, from 7,973 records.
- Translation is O(SNVs × overlapping CDS) with dictionary passes.
- The traced peak of the step is about +28 MB; it keeps about +13 MB, which are its own tables.
- There is no dense genome-wide matrix: the CDS map is a list of parts, and lookups are
  arithmetic.

**Speed-up made during this milestone.** The first version used per-group pandas operations and
took 6.8 s, with a 325 MB peak. Replacing them with plain dictionary passes cut the step to about
0.8 s and the peak to 308 MB. The output tables are identical before and after, as checked on
`mmc2`.

**Download sizes on `mmc2`:**

| File | Size |
|---|---|
| `sprint5_annotated_observations.csv` | 2.2 MB |
| `sprint5_annotated_mutations.csv` | 2.0 MB |
| `sprint5_mutation_annotations.csv` | 1.8 MB |
| `sprint5_annotation_comparison.csv` | 0.6 MB |
| `sprint5_annotation_sources.csv` | 5 KB |
| `mutation_annotation_overview.png` | 80 KB |

The largest download is still `sprint5_vaf_observations.csv` (6.22 MB), under the 8 MB per-file
cache limit.

## 11. Unsupported organisms or coding structures

- **Organisms.** No organism is hard-coded. Any pathogen works when its GenBank or GFF3
  describes the exact reference used for variant calling. Verified examples are SARS-CoV-2,
  influenza A (segmented), and mpox (both strands).
- **Marked unsupported, with the reason** (tested):
  - partial CDS (`<`/`>` or `partial=true`) or `codon_start` ≠ 1;
  - pseudogenes;
  - `/transl_except`, for example selenocysteine or a stop completed by polyadenylation;
  - any `/exception` other than ribosomal slippage, for example **RNA editing** (paramyxovirus
    P/V) or "unclassified translation discrepancy";
  - parts out of genomic order (origin-spanning or trans-spliced);
  - overlapping parts without documented slippage;
  - inconsistent GFF3 phases;
  - incomplete codons, a missing terminal stop, or internal stops;
  - a GenBank translation that disagrees with the record;
  - an unknown translation table;
  - a base read twice in a frameshifted protein;
  - a codon with an ambiguous reference base;
  - a change to an unverified start codon.
- **Not handled in this phase:** insertions, deletions, multi-nucleotide variants, and IUPAC
  alleles. They are reported as unsupported, never reinterpreted as substitutions.
- **Reference FASTA.** It must still be a single sequence, as in Sprints 2–3. Segmented
  references are therefore annotated through a multi-record GenBank file, or a GFF3 file with
  `##FASTA`. Only one segment can additionally be cross-checked against an uploaded FASTA.

## 12. Biological interpretation limitations

- A **predicted** consequence describes the sequence change only.
  - Synonymous changes can still matter (RNA structure, codon usage, regulatory elements).
  - Missense changes do not by themselves change transmissibility, severity, host range, immune
    or vaccine escape, drug resistance, replication, or outbreak risk.
  - None of these claims is made anywhere in the output.
  - No curated functional descriptions were added, and no AI is used.
- **Pooled data.** Each SNV is predicted alone against the reference. Short reads from a pooled
  sample do not establish which SNVs occur on the same genome. Shared-codon SNVs are flagged,
  and no combined amino-acid change is reported.
- **Start lost** follows the Sequence Ontology definition (any change to the start codon). The
  model does not predict whether translation starts elsewhere.
- **Mature-peptide names and positions** (for example nsp12 P323L) come from the annotation and
  are shown only for peptides that lie in frame within the CDS.
- **The workbook's published annotation** is an independent check, not ground truth for
  function.
- **Annotation does not identify lineages** and does not validate the exploratory
  mutation-frequency groups. Real-world validation of grouping remains a separate research task
  (Phase 3 report §19).

## 13. Render testing instructions (after review; not done here)

1. Deploy the reviewed branch only after approval. `requirements.txt` and the Dockerfile are
   unchanged; `examples/` (with `examples/annotation/`) is already copied into the image.
2. **Pooled with annotation.** Upload `mmc2.xlsx` together with `NC_045512.2.gb` (the upload
   page's "Example file" link). Confirm:
   - the review row "Mutation annotation and protein consequences: Will run";
   - in section 9: 5,677 evaluated, 5,677 annotated, 4,109 missense, 1,243 synonymous;
   - the independent-check sentence: 5,613 agree.
3. **Without an annotation.** Upload `mmc2.xlsx` alone. Confirm that the panel explains the
   missing annotation and that all other sections are unchanged.
4. **Mismatched annotation.** Upload a mismatched annotation, for example an influenza GenBank
   file with SARS-CoV-2 variants. Confirm the review warning and "Unverified or unannotatable"
   equal to all SNVs.
5. **Load.** Watch Render's memory graph and request time. Locally, one run with annotation
   peaks at 308 MB and takes about 1 s longer than without.

## 14. Readiness for the comprehensive Sprint 3 + 5 audit

**Ready:**
- Every Phase 4 output is traceable to its inputs, with:
  - the source file and format;
  - per-CDS checks and reasons;
  - per-SNV reference checks;
  - per-observation evidence copied from Phase 1–2.
- Reproducible commands are in §9 and §13. The real-data numbers in §7 can be regenerated with
  `run_integrated_analysis.py --snv mmc2.xlsx --annotation examples/annotation/NC_045512.2.gb`.

**For the auditor to review:**
1. The shared-codon rule uses co-reporting in the same pooled sample as the trigger. It is a
   conservative flag, not a linkage test.
2. The `mmc2` reference label NC_045512.2 still comes from the paper's methods. Phase 4 checks
   every reference base against RefSeq, but it does not change the Sprint 3 "Supplied but
   unverified" status, which concerns the uploaded FASTA.
3. Hosted (Render) behaviour is unverified.
4. Indel annotation and lineage-level interpretation remain out of scope.

## Mandatory completion criteria

| Criterion | Status | Evidence |
|---|---|---|
| Reference compatibility is checked | **Met** | Exact accession.version, SHA-256 identity with an uploaded FASTA, GFF3 length check, and a per-SNV reference base (§2–3); tests for a version mismatch, another pathogen, a changed sequence, and a wrong Ref |
| Coding annotations come from a supported annotation source | **Met** | GenBank or GFF3 only; no gene inferred from sequence; official NCBI RefSeq files (§2, §5) |
| Synonymous and nonsynonymous consequences correctly distinguished | **Met** | 5,613 of 5,613 real coding SNVs agree with the published codons and amino acids; hand-derived tests for all six coding categories; noncoding never synonymous (§6–7) |
| Strand, phase, and coordinate handling verified | **Met** | Minus-strand and spliced hand cases; 96 mpox minus-strand CDS reproduce NCBI proteins; GFF3 phase checks; frameshift offsets; mature-peptide spans equal NCBI protein coordinates (§4, §6) |
| Unsupported variants marked honestly | **Met** | Indels, ambiguous bases, read-twice frameshift bases, and unsupported structures (partial, editing, internal stop, translation mismatch, order, phase) are reported with reasons (§3, §11) |
| Pooled SNVs not interpreted as complete haplotypes | **Met** | Each SNV predicted alone; shared codons flagged and never combined (tested with N R203K/R203R); page and limitations text (§3, §12) |
| Validated annotations integrated with existing VAF and temporal results | **Met** | Per-SNV evidence and temporal summaries; per-observation join equal to `sprint5_vaf_observations.csv` column by column (web test) (§7–8) |
| Existing Sprint 2–5 outputs remain functional | **Met** | 62 of 64 and 40 of 42 output files identical, including every figure; only one status row and one provenance row added; all 585 pre-existing tests pass (§9) |
| Experimental grouping remains disabled by default | **Met** | `DEFAULT_SETTINGS["grouping_mode"] == "profiles"` tested; group tables return 404 in default runs with annotation; grouping code unchanged (§8) |
| All tests pass without unexplained failures | **Met** | `Ran 636 tests — OK (skipped=5)`, all skips optional and explained (§9) |
| Real workbook analysis completes | **Met** | `mmc2.xlsx` with and without the NCBI GenBank file, through `/review` → `/run`; md5 unchanged (§7) |
| Performance documented | **Met** | §10: +12 MB peak and about +1 s per run; 512 MB / 10% CPU: 302 MB, 131.9 s; download and HTML sizes |
| Implementation report finished | **Met** | This document |
| No commits, pushes, deployments, or paid API calls | **Met** | `git status` shows uncommitted changes only; AI disabled; nothing deployed |

**Not verified:** hosted (Render) behaviour. The browser preview of the page was refused by the in-app browser, so the
page was checked through the real Flask routes (rendered text, counts, and all 47 download links) and by viewing the
figure PNG directly.

**SPRINT 5 PHASE 4 COMPLETE**
