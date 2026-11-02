# Integration Implementation Log (Sprint 2 + Sprint 3)

This log records what was actually implemented to integrate the Sprint 2
individual-sequence pipeline and the Sprint 3 pooled-sequencing work into one
analysis, the research facts it relies on, and the remaining limitations. It
is not Sprint 4 work; Sprint 4 requirements had not been released. Newest
entries are at the bottom of the progress section.

## Research facts used (verified from the sources)

**Fontenele et al. (2021), Water Research 205:117710, PMC8464352** (read from the
Europe PMC full text, not from memory):

- 91 RT-qPCR-positive 24-h composite wastewater samples from 11 states, collected
  7 Apr to 16 Jun 2020, plus one negative control. Arizona location IDs include
  Guadalupe, G2, M1, Rural, and TP01–TP06; Kentucky (Louisville) IDs are S1–S9.
- Illumina HiSeq 2500, 2x100 paired-end reads. Reads were mapped with BWA-MEM
  to MN908947 (RefSeq NC_045512.2). Primers were soft-clipped with `ivar trim`
  (reads < 30 nt removed), with iVar's default sliding-window quality trimming
  (window 4, minimum PHRED 20).
- Coverage QC: positions with at least 10x depth at minimum base quality 20.
  Detailed analysis used samples with **> 90% genome coverage**: 52 of 91.
- Variant calling: iVar, minimum depth **20x**, minimum base quality **20**, **no
  minimum frequency** (population-level data), retained calls with iVar's
  Fisher's exact test **p < 0.05**. Problematic sites suggested by De Maio et al.
  (2020, as of September 2020) were removed; the site list is not in the paper.
- 5,680 unique SNVs (4,372 nonsynonymous, 1,308 synonymous, 246 nonsense, 64
  noncoding). "Novel" meant absent from GISAID clinical genomes on 17 Jun 2020
  (45,836 genomes) and checked again on 20 Nov 2020 (203,741 genomes); SNVs seen
  in only one wastewater sample were removed from the novel set.
- **PCoA (Section 2.6):** each sample is a matrix with one row per reference
  position and four columns (A, C, G, T frequencies). At each position the
  Yue–Clayton dissimilarity is computed between two samples; the per-position
  values are **summed across all genomic positions** into a pairwise distance
  matrix, followed by PCoA (Gower 1966). **If a frequency cannot be calculated at
  a position (zero depth), that position's dissimilarity is set to 0.** The paper
  does not print the Yue–Clayton formula; this project uses the standard
  definition from Yue & Clayton (2005): similarity = sum(a*b) / (sum(a^2) +
  sum(b^2) - sum(a*b)), dissimilarity = 1 - similarity.
- Fig 2A: SNVs per sample colored by state. Fig 2B: regression with 95% CI of
  SNV count against mean depth (model details not stated). Fig 3: novel SNVs by
  genome position, y = number of samples, red = positions with multiple
  variants. Fig 5: PCoA colored by state; Arizona PCoA with a gradient.
- Lineage comparison used GISAID genomes; the paper states that SNV linkage is
  not possible with short reads, so lineage support does not prove presence.
- Raw data: SRA BioProject PRJNA662596 (SRR12618464–SRR12618554, SRR13289969).
- The paper is internally inconsistent in a few places (11 vs 10 states; 10x vs
  ">10 reads"; 52/91 = 57% vs "~55%").

**Supplementary Table 2 (`mmc2.xlsx`)** headers, as documented by the Sprint 3
loader and README (the workbook itself is not in this repository): `SAMPLE`,
`COLLECTION_DATE`, `STATE`, `LOCATION`, `POSITION`, `REF`, `ALT`, `ALT_FREQ`,
`TOTAL_DP`, `REF_CODON`, `REF_AA`, `ALT_CODON`, `ALT_AA`, `BIN`,
`GLOBAL_FREQUENCY`, `USA_FREQUENCY`, `Sy/Ns`, `Annotation`,
`Novel SNVs based on Global frequency`,
`Identified in clinical derived-genomes after June 17th`. It has SNV records
only: no base quality, no Fisher p-values, no per-position coverage, and no
reference-allele counts at unreported positions.

**ZooPhy/radx** (inspected via the GitHub API and raw files): Python wrapper
around `bwa mem`, `samtools`, `ivar trim`, `ivar variants`, LoFreq, and Freyja.
**It has no license file**, so no code was copied. It contains no Yue–Clayton
or PCoA code and no link to the Fontenele paper. It was used only to confirm
the iVar TSV column names (`REGION, POS, REF, ALT, REF_DP, ALT_DP, ALT_FREQ,
TOTAL_DP, PVAL, PASS, REF_CODON, REF_AA, ALT_CODON, ALT_AA, ...`) and the
general mapping/variant-calling stage order.

## Progress

### Phase 0 — Branch integration (complete)

- Created local backup branch `backup/sprint4-before-integration` at `eb72013`.
- Merged `origin/sprint-3` (`1c2d59b`) into `sprint4-integration` with a merge
  commit; excluded the generated `tmp/pdfs/` guide scripts and PNGs and added
  `tmp/` to `.gitignore`.
- Adopted the uncommitted README refresh and changelog entry from the original
  `../sprint_2` working folder; committed `CLAUDE.md` and `PROJECT_HANDOFF.md`.
- Made the suite self-contained: the invented DEMO/NCBI-style example files are
  now in `tests/fixtures/`; the Flask state-only upload test generates its own
  records; the real-NCBI version is kept as an optional, clearly skipped
  integration test (`tests/test_optional_real_ncbi_data.py`).
- Result: 111 tests, all passing, 1 optional test skipped. With the untracked
  real files temporarily present, the optional test also passed.

### Phase 1 — Shared reference-based model (complete)

- `genomics/references.py`: one-record reference FASTA, SHA-256 checksum;
  different labels are accepted as one reference only when both supplied
  sequences are identical. Labels alone (MN908947.3 vs NC_045512.2) are not.
- `genomics/unified_schema.py`: shared sample and variant tables,
  `Variant_Key = Reference_ID:Position:Ref>Alt`, explicit status words
  (variant observed, reference base called, adequate depth not reported, low
  coverage, not covered, coverage unknown, reported but excluded by QC).
- `genomics/pooled_observations.py`: processed SNV rows -> shared table;
  optional per-position depth / A/C/G/T table; unlisted positions are unknown.

### Phase 2 — Reference-coordinate individual profiles (complete)

- `genomics/sequence_variants.py`: MAFFT aligns each sequence to the reference;
  the reference row defines 1-based positions; partial spans, N/IUPAC,
  deletions, insertions handled separately; mapping safeguard (>= 20 callable
  positions, >= 90% identity — prototype rule, not from the paper).
- Backbone changes (documented and tested): `align_sequences.run_mafft()`
  factored out with the identical command; FASTA validation accepts IUPAC
  ambiguity codes, which the p-distance treats exactly like N.

### Phase 3 — Pooled importer and paper QC (complete)

- `wastewater/processed_snv.py` keeps Supplementary Table 2 annotations as
  `Published_*`, reads iVar TSV (indels set aside visibly; `Source_*` columns),
  validates alternate read counts, records its input format.
- `genomics/quality_control.py`: paper preset (>90% coverage at 10x, depth >= 20,
  Fisher p < 0.05, no frequency cut-off, optional masked sites) applied only to
  supplied values; base quality reported as not verifiable; reference-base
  agreement added as a data-integrity check; paper vs exploratory modes; a
  method-reproduction report for every run.

### Phase 4 — Analyses and figures (complete)

- `genomics/temporal_analysis.py`: observation matrix (no unknown -> 0%),
  per-site trends with weekly least-squares slope, SNVs per sample, depth
  regression with seeded bootstrap band and permutation p-value.
- `genomics/diversity_analysis.py`: paper Yue–Clayton PCoA when A/C/G/T counts
  exist (validated on hand-computed examples; negative eigenvalues reported);
  otherwise a separately labeled exploratory Jaccard ordination.
- `genomics/figures.py`: SNVs per sample (Fig 2A analogue), depth vs SNVs
  (Fig 2B), genome positions (Fig 3), trajectories, heatmap, PCoA (Fig 5),
  cross-source categories, cross-source timeline, sequencing quality. The
  categorical palette passed the dataviz colorblind validator (contrast WARN
  for three hues is covered by legends, marker shapes, and table/CSV views).

### Phases 5–6 — Cross-source comparison and synthetic demo (complete)

- `genomics/variant_matching.py`: shared / individual-only / pooled-only /
  unresolved (either source) / incomparable, with dates, date-window pairs,
  carrier locations, precise-coordinate distances, and documented-profile
  overlap (never a lineage call).
- `examples/integrated/`: deterministic, entirely invented paired dataset with
  known answers (`EXPECTED_CATEGORIES`); `tests/test_integrated_example.py`
  checks every category, QC exclusions, trend arithmetic (0.20/week), PCoA.

### Phase 7 — Unified dashboard (complete)

- `templates/index.html` + `integrated_view.py` + `app.py`: one submission for
  individual, reference, and pooled inputs; eight sections; Sprint 2 results
  preserved verbatim; provenance; 10-minute in-memory CSV/PNG downloads;
  legacy Sprint 3 pages still linked. Rendered and inspected in a local
  browser preview with MAFFT mocked.
- `run_integrated_analysis.py`: the same analysis from the command line.

### Phase 8 — Optional offline raw-read workflow (template + converters)

- `offline_workflow/`: BWA-MEM / iVar trim / iVar variants / samtools
  template with the paper's settings (NOT run: tools not installed), and
  tested converters for iVar TSVs and `samtools mpileup` -> Q20 A/C/G/T counts.

### Phase 9 — Verification

- 164 tests: 162 passed, 2 skipped (real MAFFT absent; optional real-NCBI data).
- Mocked: MAFFT (pass-through or `tests/fake_mafft.py`) and GeoNames.
- Not verified here: a real MAFFT run, a live GeoNames call, a Docker build,
  the real `mmc2.xlsx` workbook, and the Render deployment.

### Phase 10 — Professor feedback: usability, formats, metadata, visuals (complete)

- **Checkpoint:** local tag `checkpoint/before-guided-workflow` at `ec1c017`;
  baseline 164 tests (163 passed with real MAFFT on PATH, 1 optional skipped).
- **Guided workflow:** `/review` stages uploads (30 minutes, temporary folder),
  inspects every file by content, and suggests a role only with stated
  evidence (shared IDs, reference names, explicit FASTA annotations, required
  headers). Ambiguous roles and columns must be chosen by the user. `/run`
  writes standardized copies and calls the unchanged loaders and analyses.
- **Formats:** gzip for FASTA/CSV/TSV/VCF (60 MB decompressed limit), TSV/TXT,
  SNV VCF (allele-specific AD, AO/RO, AF, iVar ALT_FREQ, single-sample INFO/AF;
  genotype-only and indel entries set aside and counted), several VCF or iVar
  files combined; FASTQ/SAM/BAM/CRAM refused with a pointer to the offline
  workflow.
- **Metadata:** role-specific header aliases; ambiguity blocks the run; unit
  checks (fraction coverage, percentage frequencies) block instead of
  rescaling; explicit FASTA `[key=value]` annotations fill absent columns with
  a recorded source; individual metadata needs only `Accession` + `Pathogen`
  (a stated pathogen is accepted on the review page); partial dates are not
  given a day.
- **Dashboard:** sections reordered (integrated first), each with purpose,
  statistics, figures, tables, limitations, and downloads; 13 new figures and
  improved genome-position and PCoA figures; per-sample cross-source evidence
  table; an analysis-status table for every run.
- **Performance:** `frequency_trends` and `observation_matrix` vectorized or
  hoisted; the real 7,973-record table runs in about 8 s instead of 106 s.
- **Real-data validation:** the real `mmc2.xlsx` (via the Europe PMC
  supplementary-files service; not in Git) matched the paper's 7,973 records,
  52 samples, 10 states, 24–793 SNVs per sample, D614G in all 52, and 4 novel
  SNVs in more than 8 samples; differences (5,677 vs 5,680 distinct SNVs, 589
  vs 548 multi-sample novel SNVs, 12 records below 20x) are reported in
  `ARTICLE_METHOD_VALIDATION.md`.
- **Real tools:** MAFFT (Windows build in `C:\Tools\mafft-win`) ran the
  synthetic example (same categories as the mocked tests) and the real H5N1
  data through the guided workflow (45 comparisons, 22 candidate pairs among
  9 samples, 0 strict clusters). GeoNames was not called live; Docker was not
  built (daemon not running).
- **Tests:** 204 tests, 202 passed, 2 optional real-data tests skipped by
  default (both passed when their files were supplied).

### Sprint 4 — AI literature and metadata enrichment (November 2, 2026)

Sprint 4 was requested on November 2, 2026 and builds on this integration
without changing its calculations. It adds the optional `literature/`
package and dashboard section 9; see `SPRINT4_AI_NLP_INTEGRATION_REPORT.md`
for the architecture, the extraction contract, cost controls, the synthetic
evaluation, and live-versus-mocked verification. Approved literature
locations enter the Sprint 2 date/distance rescreen only through the existing
`enrich_pairs` / `find_related_pairs` functions and the existing geographic
precision rules (state, county, and administrative areas stay too coarse).

## Remaining limitations

- Exact PCoA needs genome-wide A/C/G/T counts, which `mmc2.xlsx` does not
  contain; the paper's counts would have to be regenerated from SRA reads.
- No GISAID comparison or PANGOLIN lineage support; historical novelty is only
  shown when the source table provides it.
- Base-quality filtering, primer trimming, and the De Maio masked-site list
  cannot be verified from processed records.
- Regression model details are not stated in the paper; the bootstrap band and
  permutation p-value are this project's choices.
- The individual mapping safeguard (20 callable positions, 90% identity) is a
  prototype rule and should be tuned for real genomes.
- Pooled site labels are never geocoded; cross-source geography is contextual.
- Hosted uploads are limited to 10 MB; large tables need the command line.
- Role and column detection relies on headers, shared identifiers, and
  explicit annotations; an unusual table may need its role chosen by hand.
- VCF files without allele-specific depth or frequency values (for example
  genotype-only consensus VCFs) cannot be used as pooled frequencies.
- The exploratory PCoA of the real supplementary table explains little
  variation per axis (4.8% on axis 1) and must not be read as the paper's
  Yue–Clayton result.
- The relationship network uses a simple circular layout and is drawn only
  for 40 or fewer sequences.

## Recommended next steps (candidates for Sprint 4, once its requirements are known)

1. Build and run the Docker image, then repeat the real-MAFFT example checks
   there (done locally on Windows; not yet in Docker).
2. Enter the paper's Table 1 coverage and mean depth as pooled metadata for the
   52 samples to reproduce the coverage filter and the Fig 2B regression (the
   `mmc2.xlsx` import itself is now verified).
3. Assemble SARS-CoV-2 consensus genomes from a comparable place and period
   with NC_045512.2 for a real cross-source comparison.
4. Optionally regenerate A/C/G/T counts from SRA PRJNA662596 offline for the
   exact paper PCoA.
