# Sprint 4 Implementation Log

This log records what was actually implemented for the Sprint 4 integration,
the research facts it relies on, and the remaining limitations. Newest entries
are at the bottom of the progress section.

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
