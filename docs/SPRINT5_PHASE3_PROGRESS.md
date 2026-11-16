# Sprint 5 Phase 3 — Progress Record

Exploratory depth-aware mutation-frequency grouping, inspired by MAGOS. Branch
`sprint4-integration`, HEAD `b641533`. All Sprint 5 work is local and uncommitted. Nothing is
committed, pushed, merged, or deployed. No paid API calls are made: every run sets
`AI_ENRICHMENT_ENABLED=false`.

Attempt counter: the goal allows 18 attempts. Each milestone below records its attempt number.

---

## Milestone 0 — Baseline, research review, plan (attempt 1) — Done

### 0.1 Repository state (verified 2026-11-15)

- `git status`: 41 modified or untracked entries. These are the Sprint 5 Phase 1 and Phase 2
  files, all uncommitted. Nothing is staged and HEAD is still `b641533`.
- **Full suite**, project `.venv`, MAFFT on PATH (`C:\Tools\mafft-win`):
  `Ran 535 tests in 266.4 s — OK (skipped=6)`. The six skips are the optional NCBI files test
  plus five tests that need `FONTENELE_MMC2`, which was not set in this run.
- **The four modules that use `FONTENELE_MMC2`**, rerun with it pointing to the published
  `mmc2.xlsx` (md5 `c0a75b7fa5bde8db09746613354a9f64`, unchanged):
  - `test_optional_fontenele_table`
  - `test_sprint5_depth_aware`
  - `test_sprint5_phase2_web`
  - `test_variant_evidence`

  Result: `Ran 51 tests in 95.4 s — OK`, with no skips.
- **Conclusion.** The reported Phase 2 baseline holds: 535 tests, with only the optional NCBI
  test skipped when the workbook path is set. These are present in the working tree:
  - Wilson intervals;
  - the beta-binomial model;
  - replicate-based ρ;
  - hypothetical sensitivity.

  The constrained-memory numbers come from the Phase 2 report and are re-measured in
  Milestone 8.

### 0.2 Existing data structures reused

| Structure | Built by | What Phase 3 uses |
|---|---|---|
| Evidence table (one row per reported pooled observation) | `genomics/variant_evidence.py` | VAF, counts and their source (Supplied, Base counts, Reconstructed, Frequency only, No depth), Wilson bounds, QC, reference status, `Specimen_ID`, the Phase 2 ρ columns |
| Observation matrix (variant × sample, coverage-aware) | `genomics/temporal_analysis.observation_matrix` | Status per cell: observed, adequate depth but not reported, low coverage, not covered, coverage unknown, or excluded by QC. `Callable_Evidence` marks supported zeros from A/C/G/T counts. It is already dense: 294,788 rows for `mmc2`. |
| Known points | `genomics/vaf_trends.known_points` | Detections plus supported non-detections; unknown cells are never given a value |
| Replicate pairs and ρ | `genomics/overdispersion.py` | `replicate_pairs`, `calibrate_replicates`, ρ per specimen (upper 95% bound) |
| `score_interval`, `frequency_variance` | `genomics/beta_binomial.py` | The variance model p(1−p)/n·[1+(n−1)ρ] |

### 0.3 Source methods (read in full text in this session)

**MAGOS.** Ahmadinejad et al. 2022, *Mol Biol Evol* 39(7):msac136, PMC9260306. The verbatim
Methods were read; the supplement, which holds the J-score formula, was not.

- **VAF model.**
  - A cluster's VAFs follow Beta(α, β), with α+β = ē (mean depth) and α/(α+β) = v̄.
  - This gives var(v) = v̄(1−v̄)/(ē+1).
  - The model is binomial-like with no overdispersion parameter.
  - Depth enters continuously, with no binning.
- **Speed-up.** Variant pairs with a Euclidean VAF distance below 0.01 across samples are
  collapsed first.
- **Agglomerative tree.**
  - Repeatedly merge the closest pair.
  - The distance is d(C1,C2) = w·Σ −log P(v_i; Beta(α,β)), using the pooled Beta fitted to both
    clusters.
  - w combines 1/(m1+m2), 1/var, and 1/range; part of w is reconstructed from a garbled
    extraction.
  - This is a custom likelihood linkage, not Ward.
- **Adaptive partitioning, which sets the number of clusters.**
  - Traverse the tree from root to leaves.
  - For a clade of m variants, draw m values from its fitted Beta 1,000 times. Test whether the
    observed variance is no larger than the simulated variance (one-sided t-test).
  - Reject homogeneity at P < 0.01 and split; otherwise accept the clade.
  - The figure legend says: "A cluster is accepted if the variance is lower than or equal to the
    expected value."
- **Multi-sample mode.**
  - Fit a Beta per sample.
  - The distance is the maximum of the per-sample distances.
  - A clade is accepted only if no sample rejects.
- **Validation.**
  - Simulated tumors at 30×–500×, plus 2–4 samples at 30/50/100/300×.
  - Scored with the J-score, a precision/recall-based measure.
  - The reported J ≥ 0.80 at 30× needs a mean-VAF gap of at least 0.25 for one sample and 0.15
    at 300×.
  - Comparison methods were PyClone and SciClone.
- **Tumor-specific parts.**
  - Diploid SNVs outside CNVs.
  - Cellular prevalence = 2·v̄.
  - Purity is taken from the cluster nearest 0.5.
  - CNV correction (Eq. 6–7).
  - Subclone and evolutionary interpretation.

**Ellmen et al. 2024.** *Sci Rep* 14:20273, PMC11365997, "Reconstructing SARS-CoV-2 lineages
from mixed wastewater sequencing data". The verbatim text was read through the Discussion.

- **Objective.** De novo *lineage definitions* from many samples, which are then classified with
  Pangolin, Scorpio, or Nextclade.
- **Method.**
  1. Build per-sample mutation frequencies.
  2. Drop coverage below 20 reads.
  3. Impute missing cells with k-nearest-neighbour imputation (k = 5).
  4. Run NMF with a user-chosen number of components n, minimising ‖X − WH‖²_F.
  5. Normalise each component by its maximum and keep loadings ≥ 0.25.
- **What it lacks.**
  - No sequencing-noise model.
  - No depth weighting.
  - No time model.
  - No read co-occurrence.
- **Validation.**
  - Benchmark mixtures with known genomes, assessed qualitatively by a tree and by Nextclade
    calls.
  - No precision or recall is reported.
- **Stated limitations.**
  - n is chosen manually.
  - Imputation "likely introduces error".
  - Sub-lineages are "notoriously difficult to distinguish in wastewater".

### 0.4 What transfers and what does not

| MAGOS component | Transfer? | How it is used here |
|---|---|---|
| Hierarchical (agglomerative) grouping | **Yes** | Complete-linkage agglomeration. Complete linkage is valid for any dissimilarity; Ward and centroid linkage need Euclidean distances, which a missing-data distance is not. |
| VAF variability model | **Yes, replaced** | Beta-binomial variance p(1−p)/n·[1+(n−1)ρ] from Phases 1–2. At ρ = 0 this is binomial. MAGOS's ē+1 denominator is a binomial approximation without overdispersion. |
| Depth dependence | **Yes** | Each observation keeps its own depth. Nothing is pooled to a mean depth, because pooled samples have very uneven depth. |
| Adaptive partitioning | **Yes, re-derived** | Top-down acceptance of a clade if a homogeneity test does not reject. A design-effect-corrected Pearson χ² test of equal proportions replaces MAGOS's simulated-variance t-test. A power requirement is added, because "not rejected" is weak evidence at low depth. The threshold is validated here; MAGOS's P < 0.01 is not assumed to carry over. |
| Multi-sample comparison | **Yes, adapted** | Samples are the dimensions of a profile. Each pair is compared only over samples where both values are known, with a minimum overlap. |
| Simulation-based validation | **Yes** | Synthetic pathogen mixtures with known groups, negative controls, and pre-registered criteria. MAGOS's J ≥ 0.80 is not adopted as a criterion. |
| Diploid copy-number expectations, VAF = CP/2 | **No** | Pathogen pools have no diploid expectation. A frequency is the share of reads, not half a cell fraction. |
| Tumor purity | **No** | There is no normal-cell admixture. |
| Cellular prevalence | **No** | The analogue would be a lineage's share of the pool, which is not claimed. |
| Subclone interpretation and tumor evolutionary relationships | **No** | A group is "mutations with similar frequency patterns". It does not claim a genome, lineage, haplotype, or ancestry. |
| CNV correction | **No** | Not applicable. |

**Ellmen et al. compared with this goal.**

| | Ellmen et al. | Phase 3 |
|---|---|---|
| Aim | Reconstructs lineage definitions | Groups co-varying mutations; does not name lineages |
| Missing data | Imputes missing values | Never imputes; unknown stays unknown |
| Sequencing noise | No model | Explicit binomial or beta-binomial model |
| Number of groups | Chosen by the user | Set by a tested criterion |

What carries over from Ellmen et al.: mutations that rise and fall together across samples are
candidates for sharing a source. They also warn that shared mutations and imputation create
errors, and that missing coverage is a central problem.

### 0.5 Feasibility on the published `mmc2.xlsx` (verified with the CLI)

| Quantity | Value |
|---|---|
| Distinct mutations (Reference + position + Ref + Alt) | 5,669 (one reference, NC_045512.2; reference unverified because no FASTA was supplied) |
| Samples, locations, dates | 52 samples, 23 locations, 29 dates |
| Observation-matrix cells | 294,788 |
| Reported observations | 7,958 (4 excluded by QC) |
| **Unknown cells (no coverage information)** | **286,826 (97.3%)**. No per-position depth or A/C/G/T counts, so there are **no supported zeros** |
| Count source | All 7,973 records reconstructed (frequency × depth); none measured |
| Samples per mutation | 1: 4,552; 2: 743; 3: 204; 4: 69; 5: 36; ≥ 6: 65 |
| Mutations in ≥ 3 / ≥ 4 / ≥ 5 / ≥ 10 samples | 374 / 170 / 101 / 24 |
| Pairs among the 374 mutations seen in ≥ 3 samples | 69,751; median shared samples 0. 4,584 pairs share ≥ 3, and 814 share ≥ 5 |
| Technical replicates | None (no `Specimen_ID`) |

**Initial assessment.** In `mmc2` every comparison between two mutations can use only the
samples where **both were reported**. A mutation below the caller's threshold in a sample is
unknown, not zero. These detection-only profiles are conditioned on detection. If a lineage's
frequency falls toward zero, its mutations drop out instead of being observed falling together.
Milestone 5 tests this condition explicitly ("sparse variant-only observations"). Whether `mmc2`
is grouped is then decided by the pre-registered rule P6 below, not by how the results look.

### 0.6 Implementation plan

1. **Milestone 1.**
   - Require an explicit `Library_ID` for technical replicates: distinct documented libraries
     of one `Specimen_ID`.
   - Detect duplicate records of one library.
   - Downgrade undocumented pairs to "unverified" and give them no ρ.
   - Write per-mutation eligibility states.
2. **Milestone 2.** `genomics/mutation_profiles.py`.
   - A long table of known observations, without a new dense matrix over all mutations.
   - Verified replicate libraries are combined into one specimen-level unit.
   - A per-mutation eligibility table.
3. **Milestone 3.** `genomics/mutation_grouping.py`.
   - Complete linkage written in NumPy, because SciPy is not a dependency.
   - The ordinary RMS-difference distance over shared samples.
   - A fixed cut height.
4. **Milestone 4.** In the same module, the depth-aware distance, the χ² homogeneity test,
   the power gate, and top-down acceptance. Hypothetical and calibrated ρ modes.
5. **Milestone 5.** `validation/simulate_mutation_groups.py` plus fast tests.
6. **Milestones 6–7.** Panel, figures, downloads, and settings in section 9; a descriptive
   comparison with Sprint 2 individual sequences.
7. **Milestones 8–9.** Regression, web tests, memory, and the report.

### 0.7 Pre-registered validation criteria (fixed before any simulation is run)

Metrics:

- the adjusted Rand index (ARI);
- normalised mutual information (NMI);
- pairwise precision and recall;
- the false-grouping rate;
- stability;
- the number excluded;
- runtime and memory.

Ungrouped and excluded mutations count as singletons in ARI and NMI. Precision is computed over
pairs placed together; when none are, it is undefined and reported as such. Every value is
reported as the mean over seeds with a 95% interval.

- **P1, false grouping in negative controls.** In a negative control (independent mutations,
  no group structure), the share of eligible mutations placed in any group must be ≤ 0.10 in
  every cell where the method does not abstain, when the method's noise assumption is correct.
- **P2, recovery.** With planted groups, depth ≥ 100×, missingness ≤ 30%, and the correct noise
  model, the depth-aware method must reach mean ARI ≥ 0.80.
- **P3, improvement.** "Improves on ordinary clustering" is claimed for a condition only if the
  paired ARI difference across seeds has a 95% interval above 0. "Worse" is reported when the
  interval is below 0; otherwise "similar".
- **P4, misspecification.** For a wrong ρ (binomial when the truth is overdispersed, and a
  hypothetical ρ too small or too large), report the false grouping, the recall, and the ARI.
  This is informative, not pass/fail. Any condition in which a wrong ρ breaks P1 must be stated
  on the page.
- **P5, abstention.** When no pair has enough comparable, resolvable data, the method must
  return no groups and give a reason. It must never produce a forced partition.
- **P6, decision rule for the website.**
  - The depth-aware method is offered as exploratory only in the input conditions where it
    met P1 and P2.
  - In conditions where it failed P1 (for example, detection-only profiles if they fail), the
    page abstains and explains why.
  - The ordinary baseline is shown only as a labelled comparison, with its own measured false
    grouping.

---

## Milestone 1 — Data and calibration eligibility (attempt 1) — Done

**Audit finding (replicates).** Phase 2 treated any two `Sample_ID`s with the same
`Specimen_ID`, date, location, and reference as technical replicates. Matching metadata does not
show that the two are *separately prepared libraries*. They could be one library listed twice, a
re-run of one library, or duplicate rows.

**Change.** `genomics/overdispersion.replicate_pairs` now also requires each sample to have a
documented `Library_ID`, the two `Library_ID`s to differ, and the records not to be exact
duplicates (identical alleles, counts, and depths). Pairs that fail these checks get the new
status `Not estimated: replicate libraries not documented (unverified)` with one of three
reasons, and receive no ρ.

**Plumbing.**
- `Library_ID` is detected by the aliases in `LIBRARY_ALIASES` (column_detection), in the SNV
  table and in the pooled metadata.
- It is read like `Specimen_ID`, one value per sample (processed_snv).
- It is kept in `PROVENANCE_COLUMNS` and attached to the evidence table.
- User-facing help text names the new requirement.

**Backward compatibility.**
- A Phase 2-style file with `Specimen_ID` but no `Library_ID` still runs. Its pairs are
  reported as unverified and keep the binomial baseline.
- The synthetic replicate demo generator now writes `Library_ID = LIB_<Sample_ID>`. The SNV
  table is identical apart from this one added column (checked with `DataFrame.equals`), and
  the base-count and FASTA files are byte-identical.
- One test helper (`evidence_rows` in `tests/test_overdispersion.py`) now gives each fixture
  sample a `Library_ID`. No assertion was changed.

**Tests.** `tests/test_replicate_verification.py` has 10 tests:
- documented pairs are usable;
- no or blank `Library_ID` → unverified;
- the same library twice → unverified;
- identical records → duplicates;
- conflicting dates still rejected;
- an unverified pair gets no ρ and a reason;
- a duplicated library file is not calibrated;
- loader (table, metadata, conflict);
- a Phase 2-style file is downgraded.

`Ran 29 tests — OK`, together with `test_overdispersion`.

**Eligibility states per mutation** (`genomics/mutation_profiles.mutation_eligibility`). These
are counted per sample from the observation matrix:
- detected;
- supported non-detection;
- counted not called;
- coverage unknown;
- low or no coverage;
- covered but value unknown;
- rejected by QC.

Also kept: the reference status, and the units known by count basis (measured, reconstructed,
or frequency only). A missing variant row is never 0.

**mmc2 feasibility** is in §0.5. Of 5,669 mutations, 170 are known in at least 4 samples, all of
them detection-only with reconstructed counts.

## Milestone 2 — Mutation-profile representation (attempt 1) — Done

`genomics/mutation_profiles.py`:
- A long table with one row per mutation and observation unit with a known value; unknowns are
  not stored.
- The identity is Reference_ID (= segment) + position + Ref + Alt.
- **Units.** A verified replicate pair becomes one unit with summed counts. It keeps
  `Depth_Square_Sum`, so the variance factor c = [n + ρ(Σn_l² − n)]/n² is exact for combined
  libraries.
- Other samples are separate units, with `Same_Site_Date_Units` flagging multiples at one site
  and date.
- **Count basis**: Measured, Reconstructed (approximate), or Frequency only.
- **Intervals**: the Wilson interval (binomial, on the unit's counts) and the calibrated
  beta-binomial interval (only for units of calibrated specimens).
- **Eligibility**:
  - the ordinary method needs at least `minimum_shared` known units;
  - the depth-aware method needs that many units with counts;
  - a reference with mismatches is excluded.

On `mmc2` the module builds 7,958 rows in 0.18 s. Dense arrays are built later, only for the
eligible mutations.

## Milestones 3–4 — Ordinary baseline and depth-aware method (attempt 1) — Implemented

These are implemented in `genomics/mutation_grouping.py` and `genomics/chi_square.py` (no SciPy).
The definitions are in the module docstring. In short:

- **Ordinary.**
  - Complete linkage on the RMS frequency difference over shared units.
  - A pair with fewer than m shared units has distance ∞, never 0.
  - The tree is cut at height 0.10.
- **Depth-aware.**
  - The pairwise design-effect-corrected Pearson statistic X²/df over shared units.
  - A power gate: there must be at least 80% power to detect a 0.10 difference at α = 0.01;
    otherwise the pair is not comparable.
  - Complete linkage.
  - MAGOS-style top-down acceptance of a clade when the group homogeneity test gives P ≥ 0.01.
- **Chi-square helpers.** Checked against SciPy (system Python): chi2 sf, isf, noncentral sf,
  and the power noncentrality agree to about 1e-12.
- **Toy check.** Three planted groups plus five unrelated mutations at 300× gave ARI 1.0 for
  both methods.

Formal tests are added with Milestone 5.

## Milestone 5 — Simulation design (fixed BEFORE running; attempt 1)

Script: `validation/simulate_mutation_groups.py` (its docstring gives the model). Seed 20261115;
30 datasets per cell.

1. **Minimum shared units, chosen first.**
   - Candidates m ∈ {2, 3, 4, 6}, run on structure (planted, negative, look-alike) × depth
     (50, 300) × missingness (0.3, 0.6), with ρ = 0.
   - **Rule:** the smallest m for which the depth-aware method (correct binomial model) meets P1
     in every negative and look-alike cell.
   - The ordinary method uses the same m.
2. **Main grid.**
   - Structure (3) × depth {20, 50, 100, 300, mixed} × true ρ {0, 0.02} × missingness
     {0, 0.3, 0.6}.
   - Methods:
     - ordinary h = 0.10, plus h = 0.05 and 0.15 to test the baseline's own sensitivity;
     - depth-aware with assumed ρ ∈ {0, 0.002, 0.02, 0.05}, so one assumed ρ equals the truth
       and the others are wrong;
     - no grouping.
   - Stability is measured by leave-one-unit-out reruns: planted structure, depth 100 or 300,
     missingness ≤ 0.3, 10 seeds.
3. **Detection-only grid.** The same with missingness 0 and no supported zeros; caller threshold
   ≥ 2 reads and ≥ 2%.
4. **Calibration grid.**
   - Two libraries per unit at per-library depth {100, 300, 1000} × true ρ {0.005, 0.02} ×
     missingness {0, 0.3} × structure (3).
   - Methods: ordinary, binomial, true ρ, and ρ calibrated per unit from the replicates (Phase 2
     estimator, upper bound, 1,000 resamples).
   - Calibrated mode is run only if every unit has an estimate; otherwise it abstains and the
     case is counted.
5. **Units and size study.** Units {3, 6, 12, 24} × scale {1, 3}, at depth 100, ρ = 0, 30% missing.
6. **Null calibration of the homogeneity test.**
   - Two mutations with the same true VAF {0.02, 0.05, 0.2, 0.5} × depth × ρ.
   - Type I error against α = 0.01, from 2,000 replicates.
7. **Runtime and memory.** V, U ∈ {(100, 12), (300, 52), (1000, 52), (1000, 192)}.

P1–P6 (§0.7) apply unchanged. For P3, "improves" is judged on the planted structure. In
negative and look-alike controls ARI degenerates, because the truth is all singletons, so the
false-grouping share is the metric there.

### Attempt 1 result: FAILED P1 (minimum-overlap study, development seeds)

File: `docs/sprint5_phase3_results/attempt1_grouping_minimum_shared.csv`. Means over 30 seeds:

| Control | Depth-aware attempt 1, false-grouping share (m = 2 / 3 / 4 / 6) | Ordinary h = 0.10 |
|---|---|---|
| Negative, 300×, 30% missing | 0.152 / 0.140 / 0.115 / 0.037 | 0.55 / 0.53 / 0.50 / 0.34 |
| Negative, 300×, 60% missing | 0.391 / 0.145 / 0.028 / 0.000 | 0.66 / 0.38 / 0.13 / 0.01 |
| Look-alike, 300×, 30% missing | 0.812 / 0.801 / 0.780 / 0.580 | 0.99 / 0.99 / 0.98 / 0.81 |
| Look-alike, 300×, 60% missing | 0.700 / 0.510 / 0.239 / 0.019 | 0.86 / 0.64 / 0.35 / 0.04 |

**No candidate m met P1.** The main run was stopped before any of its results were viewed, so
the main grid stays untouched confirmation data.

**Diagnosis.** The acceptance step is "the homogeneity test does not reject". The power gate only
guarantees that a 0.10 difference is detectable. Look-alike pairs differ by about 0.04 per unit,
and many negative-control trajectories are nearly flat at similar low levels. Neither can be
distinguished from a true group, because they agree in *level*, not in a *pattern*. Agreement
between two flat profiles carries no information that the mutations vary together.

**Attempt 2 change (method only; simulation conditions unchanged).** A clade must also have an
*informative pattern*: across the units where at least two members are known, its pooled
frequency must span at least 0.20, which is two resolution steps (`PATTERN_RANGE`). This is a
stated heuristic, not derived from a likelihood. Its cost is lower recall for groups with flat
frequencies.

**Ablation.** The same requirement is added to the ordinary method ("Ordinary + pattern
requirement", simulation only) so its effect can be separated from depth-awareness.

### Attempt 2 minimum-overlap study (same development seeds)

File: `docs/sprint5_phase3_results/attempt2_grouping_minimum_shared.csv`.

- The depth-aware method meets P1 in every negative and look-alike cell for every m. The
  largest share is 0.092 (negative, 300×, 60% missing, m = 2).
- **The pre-registered rule selects m = 2**, the smallest candidate, so
  `DEFAULT_MINIMUM_SHARED = 2`.
- The ablation still groups 0.20–0.35 of negative-control mutations, and up to 0.13 of
  look-alikes at 50×. **The pattern requirement alone does not explain the depth-aware
  method's control behaviour.**
- Planted structure, 300×, 30% missing: depth-aware ARI 0.825 versus 0.823 for the ablation and
  0.901 for ordinary (the latter from attempt 1's table; unchanged method). False grouping:
  0.053 depth-aware versus 0.22–0.34 for the two ordinary variants.
- At 50× the depth-aware method mostly abstains: ARI 0.06, abstaining in 43% of datasets at 30%
  missing.

The main grid (held-out seed streams) is now run once with attempt 2.

### Attempt 2: held-out results (main, detection-only, calibration, units, null, runtime)

The run took 282 s; files are `docs/sprint5_phase3_results/attempt2_*.csv`. Full tables are in
report §9–10.

- **P1** (correct noise model):
  - Binomial truth: **fails in 6 of 30 control cells**, at most 0.243.
  - ρ = 0.02 with assumed ρ = 0.02: met, at most 0.059, but mostly because the method abstains.
  - Calibrated mode: fails in 4 of the 8 cells where calibration was possible, at most 0.228.
  - More mutations raise false grouping: at 12 units, 0.19 with 29 mutations and 0.33 with 87.
- **P2:**
  - Met at 100×/0%, 300×/0%, and 300×/30% with binomial truth (ARI 0.91, 0.91, 0.81).
  - Not met at 100×/30% (0.40), nor in any overdispersed cell (≤ 0.26).
- **P3:** ARI is worse than ordinary clustering in 27 of 30 count-based conditions, similar in 2,
  and better in 1 (mixed depth, no missing: +0.18). Precision and false grouping are always
  better, by about 3–10×.
- **P4:**
  - A binomial assumption under ρ = 0.02 gives false grouping up to 0.28 (look-alike) and low
    recall.
  - ρ = 0.05 assumed always abstains.
- **P5:** abstention with a reason was verified (low depth, overdispersion, flat profiles,
  detection-only input).
- **Detection-only input:**
  - ordinary clustering 0.18–1.00;
  - depth-aware 0–0.18 with the correct model, failing P1 at 100× and mixed depth.
- **Null test:** type I error 0–0.0155, against a nominal 0.01.
- **Runtime:** depth-aware 1.3 s at 1,000 × 52 and 22.7 s at 1,000 × 192. `PAIR_BUDGET` was
  lowered to 5×10⁷, so about 1–2 s at most.

**Decision (P6, applied as pre-registered; no further tuning on held-out data).** No input
condition that the website can observe is one where P1 and P2 both hold. The ρ = 0 cases cannot
be verified without calibration, and the calibrated cells are mixed. So:

- the default mode ("profiles") shows profiles, eligibility, and an ungrouped heatmap, with no
  groups;
- groups are available only as an opt-in experimental view labelled "not validated", with the
  measured error rates on the page;
- detection-only input never receives groups.

This follows the goal's instruction to keep an unreliable algorithm experimental and continue
with a defensible baseline. For this purpose the defensible baseline is abstention, not
ordinary clustering, which failed P1 everywhere.

## Milestone 5 — Done (attempt 2; attempt 1 failure documented)

The method is documented as not meeting its criteria and kept experimental.
`tests/test_mutation_grouping_simulation.py` (8 tests) pins the main findings, including the
attempt 1 failure and the overdispersion abstention.

## Milestone 6 — Website integration — Done

- **Panel** "Exploratory Mutation-Frequency Groups" in section 9, after the Phase 2 panel:
  - the required sentence;
  - a mode badge and the reasons;
  - statistics: evaluated, eligible, excluded, groups, method, calibration;
  - the validation summary;
  - figures, a methods table, collapsible details, limitations, and downloads.
- **Settings:** a grouping mode (profiles, experimental, off) and the minimum shared samples
  (2–20). Both are validated and kept on re-runs.
- **Figures:**
  - `mutation_profile_heatmap` and `mutation_group_heatmap` (at most 60 rows, disclosed);
  - `mutation_group_trajectories` (at most 6 groups);
  - `grouping_comparison`.
- **Downloads:** profiles, eligibility, parameters, assignments, summary, trajectories,
  comparison, and the Sprint 2 overlap.
- **Synthetic mixture demo** (`examples/integrated/make_mixture_example.py`; calibrated mode):
  - depth-aware recovers groups A and B exactly with no false members (ARI 0.97, precision 1.0)
    and leaves the flat group C ungrouped;
  - ordinary clustering has ARI 0.91 with 4 non-group mutations in its groups.
- `tests/test_sprint5_phase3_web.py` (11 tests, including `mmc2`) passes.

## Milestone 7 — Integration with existing surveillance — Done

- **Output comparison** with the saved Phase 2 outputs:
  - `--example` (with MAFFT): 58 of 61 files identical;
  - `mmc2`: 35 of 38 files identical;
  - all figures are byte-identical.
- **Differences**, all intended: two status rows, four provenance parameters, and `Library_ID`
  appended as the last evidence column.
- Its first placement shifted the Phase 2 columns; this was caught by the comparison and fixed.
- **Sprint 2 overlap** (`sprint5_group_individual_overlap.csv`) is descriptive only, with a fixed
  disclaimer on every row.
- `mmc2` download sizes: largest 6.2 MB, under the 8 MB limit.

## Milestone 8 — Performance and full regression — Done

**Memory regression found and fixed.** The first measurement of the default mode on `mmc2` peaked
at 343–350 MB, against 300 MB with Phase 3 off. Tracing each step (`tracemalloc` inside the real
run) found two causes:

- the `imshow` heatmap: +64 MB of full-canvas resampling buffers. It now uses `pcolormesh`, which
  draws the same image.
- string copies in `mutation_eligibility`: +13 MB. It now counts with integer codes and
  `bincount`, and the output is identical on `mmc2`.

After the fixes:

| Measurement | Result |
|---|---|
| Single run, default mode | 297–300 MB (off: 299–300 MB) |
| Single run, experimental mode | 300 MB |
| Four runs | 318–325 MB on, 322 MB off; the growth is the download cache, not Phase 3 |
| 512 MB cap, 10% CPU | 288 MB, 118.8 s |
| Extra time | about +1.3 s per `mmc2` run |

**Final output comparison** (final code):

- `--example`: 58 of 61 files identical.
- `mmc2`: 35 of 38 files identical.
- The differences are only the intended ones.
- `mmc2` md5 is unchanged.

**Final full suite:** `Ran 585 tests in 371.8 s — OK (skipped=4)`. The skips are 3 optional SciPy
checks, which pass under system Python, and the optional NCBI test.

## Milestone 9 — Report — Done

`docs/SPRINT5_PHASE3_IMPLEMENTATION_REPORT.md` (20 sections plus the completion table) ends with
**SPRINT 5 PHASE 3 COMPLETE**. The uncertainty-aware method is reported as not meeting its
pre-registered criteria and kept experimental.

Attempts used: 2 of 18 for the method (attempt 1 failed P1; attempt 2 partly failed). Nothing was
committed, pushed, or deployed, and no paid calls were made.
