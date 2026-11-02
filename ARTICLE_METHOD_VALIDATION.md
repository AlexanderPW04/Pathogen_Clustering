# Method-to-implementation validation: Fontenele et al. (2021)

Fontenele RS et al. *High-throughput sequencing of SARS-CoV-2 in wastewater
provides insights into circulating variants.* Water Research 205:117710 (2021).
PMC8464352, doi:10.1016/j.watres.2021.117710.

This table compares what the study did (Methods 2.1–2.6 and Figures 1–5, read
from the full text on November 1, 2026) with what this application implements.
It is part of the Sprint 2 + Sprint 3 integration work, not Sprint 4.

**Status words**

- **Faithful**: the paper's calculation is reproduced exactly when the stated inputs are supplied.
- **Partial**: the same idea is applied, but only to the values the inputs contain, or with a documented
  substitute where the paper does not state the method.
- **Exploratory**: a clearly labeled alternative, not the paper's method.
- **Not reproduced**: outside what processed uploads can support.

None of the application's figures are copies of the paper's figures; they are
drawn only from uploaded data.

## Methods

| Study step (what the paper did) | What this application implements | Status | Inputs required | What has been tested | Remaining limitations |
|---|---|---|---|---|---|
| 2.1–2.3 Sampling, concentration, RT-qPCR, library prep, Illumina HiSeq 2500 2x100 | Nothing (laboratory steps) | Not reproduced | — | — | Laboratory and sequencing steps are outside a data application. |
| 2.4 Read mapping: BWA-MEM to MN908947 / NC_045512.2 | Optional offline command template (`offline_workflow/`); the website imports processed results only | Not reproduced in the website | Raw FASTQ (offline) | Converters for iVar TSV and `samtools mpileup` are unit-tested; the template itself was not run (tools not installed) | Importing processed SNV calls is **not** equivalent to running read mapping and variant calling. |
| 2.4 Primer soft-clipping and quality trimming (iVar trim, window 4, Q20, reads < 30 nt removed) | Recorded in the offline template; marked "Not reproduced here" in every run's method table | Not reproduced | Raw reads and the primer scheme | Method table text checked in tests | Cannot be verified from processed records. |
| 2.4 Genome coverage at ≥10x and Q20; detailed analysis only for samples > 90% covered | `Paper coverage filter`: > 90% at ≥ 10x applied when genome coverage is supplied in pooled metadata or calculable from a complete per-position depth table | Partial (Faithful for the threshold; Q20 not verifiable) | Pooled metadata `Genome_Coverage_Percent`, or a per-position `Depth` table covering every reference position | Synthetic example: WW_S4 (78%) excluded in paper mode, flagged in exploratory mode | The paper's Table 1 coverage values are not loaded automatically; the published SNV table has no coverage column (12 of its records are below 20x, see below). |
| 2.4 iVar variant calling: minimum depth 20x, minimum base quality 20, no frequency cut-off | Depth ≥ 20 applied to supplied `Read_Depth`; no frequency cut-off; base quality reported as "not verifiable" | Partial | `Read_Depth` per record | Unit tests for every QC check; real table: 7,961 of 7,973 records retained | Base quality cannot be checked from processed records. |
| 2.4 Fisher's exact test p < 0.05 (iVar PVAL / PASS) | Applied to supplied p-values or iVar PASS flags; never recalculated | Partial | iVar `PVAL` or `PASS` columns | Unit tests; VCF `FILTER` is recorded but not treated as the Fisher test | The published supplementary table has no p-values, so this step is "not evaluable" for it. |
| 2.4 Removal of De Maio et al. (2020) problematic sites | Applied to an uploaded masked-site list | Partial | A masked-site table | Unit tests | The historical September 2020 list is not bundled or invented. |
| 2.4 SNV annotation (synonymous / non-synonymous / nonsense / non-coding) | Published annotations kept as `Published_*` columns; never recalculated | Partial (display only) | The published table or caller annotations (iVar `GFF_FEATURE`, codons) | Real table: 4,368 Nsyn / 1,309 Syn distinct SNVs read from the published column | No independent annotation is computed. |
| 2.4 Novelty against GISAID (17 Jun and 20 Nov 2020), singletons removed | Published novelty flag shown only when the source table provides it | Partial (display only) | Published table | Real table: 3,701 distinct SNVs flagged "novel"; 589 of them in two or more samples (paper: 548) | GISAID comparison is not possible without GISAID access; an SNV missing from an uploaded comparison set is never labeled "novel". |
| 2.5 PANGOLIN lineage support from GISAID genomes | Optional documented variant profiles compared as profile overlap only | Exploratory | A documented profile table with sources | Synthetic example: overlap counts with "not a lineage assignment" label | Short reads cannot link SNVs, so overlap never confirms a lineage (the paper states the same caveat). |
| 2.6 Per-position A/C/G/T frequency matrices; Yue–Clayton dissimilarity per position, summed over the genome; zero-depth positions contribute 0; PCoA (Gower 1966) | Implemented exactly as described, with the standard Yue & Clayton (2005) formula; negative eigenvalues reported, no correction | Faithful (when counts are supplied) | Per-position A/C/G/T counts for every pooled sample | Hand-computed unit examples; synthetic example PCoA on 7 samples | The paper's own counts would need regenerating from SRA PRJNA662596 (offline workflow). |
| 2.6 when only SNV calls are available | Jaccard distances between reported SNV sets, then PCoA, titled "EXPLORATORY – not the paper's Yue–Clayton PCoA" | Exploratory | Pooled SNV table | Real table: 52 samples ordinated; the label is asserted in tests | Sequencing depth can drive these distances; unreported SNVs are not absences. |

## Figures

| Paper figure (what it shows) | Application figure | Status | Inputs required | Tested | Limitations |
|---|---|---|---|---|---|
| Fig 1A: map of states that supplied samples | `Geographic precision` (precision classes + plain latitude/longitude plot) and `Sampling dates and sites` | Partial (no map tiles) | Locations; coordinates for the plot | Rendering tests; real table (sites shown as named sites without coordinates) | Pooled site labels are never geocoded; broad areas are never drawn as points. |
| Fig 1B: Ct value and coverage uniformity per sample | Not drawn | Not reproduced | Ct values (not in the processed tables) | — | Could be added if Ct values were supplied as metadata. |
| Fig 1C: genome coverage of all samples, > 90% highlighted | `Sequencing quality by source` (coverage bars with the 90% line; excluded samples hatched) | Faithful when coverage is supplied | Genome coverage per pooled sample | Synthetic example | The published SNV table has no coverage, so the real-table run shows "unknown". |
| Fig 2A: SNVs per sample, colored by state | `SNVs per pooled sample` (colored by sampling site) | Faithful | Pooled SNV table | Real table: per-sample counts 24 to 793, as in the paper | Colors follow sites; with more than 8 sites the extra sites share a gray color (legend and CSV give exact values). |
| Fig 2B: regression of SNV count on mean depth with 95% CI | `Sequencing depth versus detected SNVs`: least-squares line, seeded bootstrap 95% band, permutation p-value | Partial (model details not stated in the paper) | Mean depth per pooled sample | Unit tests of the regression; synthetic example | Needs mean depth (paper Table 1), which the published SNV table lacks. |
| Fig 3: novel SNVs by genome position; y = samples; red = positions with several variants; grey = later seen clinically | `Variants along the reference genome`: y = samples, color = alternative base, shape = source, ring = several pooled alternative bases, star = published novelty flag | Partial (all SNVs, not only novel; clinical-follow-up flag not drawn) | Pooled SNV table and/or individual profiles on one reference | Real table: 5,669 retained distinct SNVs drawn; synthetic example | Annotation bands span annotated records, not full gene boundaries. |
| Fig 4: lineage support and days to first clinical observation | `Overlap with documented variant profiles` table and `When variants were seen in each source` (cross-source dates) | Exploratory | Documented profiles; individual sequences for dates | Synthetic example | No GISAID or PANGOLIN; dates are context, not evidence of which source detected a change first. |
| Fig 5A/C: PCoA colored by state; Arizona over time | `Population structure (PCoA)`: color = collection date, marker = site (or state when there are more than 8 sites), eigenvalue panel | Faithful with A/C/G/T counts; Exploratory otherwise | Per-position A/C/G/T counts | Synthetic example (paper method); real table (exploratory) | The real-table ordination is exploratory: its first axis explains 4.8% of the positive eigenvalues. |

## Comparison with the published supplementary data

`validation/check_fontenele_2021.py` reads the real Supplementary Table 2
(`mmc2.xlsx`, obtained through the Europe PMC supplementary-files service; it
is not stored in Git) and compares it with the numbers stated in the paper.
Run on November 1, 2026:

| Quantity | Paper reports | Application calculates | Agreement |
|---|---|---|---|
| SNV records | 7,973 | 7,973 | Exact |
| Samples analyzed (> 90% coverage) | 52 | 52 | Exact |
| States | 10 | 10 | Exact |
| SNVs per sample | 24 to 793 | 24 to 793 | Exact |
| Samples with spike D614G (A23403G) | all 52 | 52 | Exact |
| Low-frequency A23403T (D614V) | sample #147 (described as Tempe) | sample 147, labeled AZ: Guadalupe in the table | Present; site label differs between text and table |
| Novel SNVs in more than 8 samples | 4 | 4 | Exact |
| Distinct SNVs (position, Ref, Alt) | 5,680 | 5,677 | Differs by 3 |
| Non-synonymous / synonymous | 4,372 / 1,308 | 4,368 / 1,309 (published column) | Close, not exact |
| Novel SNVs seen in ≥ 2 samples | 548 | 589 | Differs (the paper's singleton rule may count differently) |
| Records below the paper's own 20x minimum | 0 expected | 12 | Published table contains them; paper-mode QC excludes them |
| Collection dates of the 52 samples | (all 91 samples: 7 Apr to 16 Jun 2020) | 23 Apr to 16 Jun 2020 | Consistent subset |

The differences are reported as found; no value was adjusted to match the
paper. The paper's PCoA cannot be reproduced from this table because it has
no per-position A/C/G/T counts, so the application shows the exploratory
ordination and says so. Reproducing Fig 2B and the coverage filter would need
the paper's Table 1 (coverage and mean depth) entered as pooled metadata.

To repeat the check:

```console
python validation/check_fontenele_2021.py path\to\mmc2.xlsx
set FONTENELE_MMC2=path\to\mmc2.xlsx
python -m unittest tests.test_optional_fontenele_table -v
```
