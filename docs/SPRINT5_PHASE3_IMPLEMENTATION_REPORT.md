# Sprint 5 Phase 3 — Implementation Report

**Exploratory depth-aware mutation-frequency grouping, inspired by MAGOS.**

- Branch `sprint4-integration`, HEAD `b641533`.
- All Sprint 5 work (Phases 1–3) is local and uncommitted.
- Nothing was committed, pushed, merged, or deployed.
- No paid API calls were made: `AI_ENRICHMENT_ENABLED=false` in every run.

Labels used below:
- **verified**: run or read in this session;
- **synthetic**: simulated data with known truth;
- **real**: the published `mmc2.xlsx`.

## 1. Summary of the implementation

The research question was: *does accounting for sequencing depth and allele-frequency uncertainty
improve the grouping of mutations with similar frequency patterns, compared with ordinary
clustering?* The short, verified answer, on synthetic data:

- **It greatly reduces false grouping.** With the correct noise model, 0–24% of unrelated
  mutations were placed in a group, against 10–100% for ordinary clustering. The groups it does
  form are 69–100% pure.
- **It does not improve overall recovery.** Its adjusted Rand index (ARI) was lower than ordinary
  clustering in 27 of 30 count-based conditions, similar in 2, and higher in 1. It trades recall
  for precision: at low depth, with missing coverage, or with library-to-library overdispersion
  it mostly declines to group.
- **It did not meet the pre-registered criteria in every condition.**
  - P1, false grouping ≤ 10% in negative controls: failed in 6 of 30 binomial-truth control
    cells (5 negative, 1 look-alike).
  - P2, ARI ≥ 0.80 at depth ≥ 100× with ≤ 30% missing: met only without overdispersion and with
    little missing coverage.
  - The first version (attempt 1) failed P1 badly. It was diagnosed and revised once (§7);
    attempt 2 was then confirmed on held-out simulation seeds.

By the pre-registered rule P6, the method is therefore **not offered by default**. The website
section **"Exploratory Mutation-Frequency Groups"**:

- by default shows coverage-preserving mutation profiles, eligibility and exclusion counts, and an
  ungrouped frequency heatmap, with the reason no groups are shown;
- can show the depth-aware groups as an **opt-in experimental view (not validated)**, with
  ordinary clustering as a labelled comparison;
- never shows groups for detection-only input, such as `mmc2.xlsx`.

Other changes:

1. **Technical-replicate audit (Milestone 1).**
   - Phase 2 inferred replicates from a shared `Specimen_ID`, date, and location alone. Pairs now
     also need distinct, documented `Library_ID`s, and exact duplicate records are rejected.
   - Undocumented pairs are reported as unverified and receive no ρ. There are 10 regression
     tests for false-positive replicates.
2. **Mutation-profile representation (Milestone 2).** `genomics/mutation_profiles.py`:
   - long-form, known values only;
   - verified replicate libraries are combined into one specimen unit;
   - count provenance kept (measured, reconstructed, frequency-only);
   - per-mutation state counts;
   - unknown is never 0.
3. **Two grouping methods (Milestones 3–4).** `genomics/mutation_grouping.py` and
   `genomics/chi_square.py`, using NumPy only, with complete linkage.
4. **A pre-registered simulation study (Milestone 5).**
   `validation/simulate_mutation_groups.py`. Summary tables are in
   `docs/sprint5_phase3_results/`.
5. **Website, downloads, settings, and a synthetic mixture demo (Milestones 6–7).** Phase 1–2
   outputs are identical apart from documented additions (§15).

## 2. Review of MAGOS and the pathogen-specific method

Both papers were read in full text in this session (PMC9260306, PMC11365997).

**MAGOS** (Ahmadinejad et al. 2022, *Mol Biol Evol* 39:msac136):
- **VAF model.**
  - A cluster of m variants follows Beta(α, β), with α+β = ē (mean depth) and α/(α+β) = v̄.
  - This gives var(v) = v̄(1−v̄)/(ē+1). It is binomial-like, with no overdispersion parameter.
- **Speed-up.** Variant pairs closer than 0.01 in Euclidean VAF distance are collapsed first.
- **Agglomerative tree.** A custom Beta-likelihood distance, with weights 1/(m1+m2), 1/var, and
  1/range. It is not Ward linkage.
- **Adaptive partitioning, root to leaf.**
  - For a clade, simulate m Beta draws 1,000 times.
  - Run a one-sided t-test of whether the observed variance is ≤ the expected variance.
  - Split if P < 0.01; otherwise accept the clade.
- **Multi-sample mode.** The distance is the maximum over samples, and a clade is accepted only if
  no sample rejects.
- **Validation.** J-score at 30–500×. A J ≥ 0.80 at 30× needs a mean-VAF gap of 0.25.
- **Tumor-specific parts.** Diploid SNVs, cellular prevalence = 2·v̄, purity, and a CNV
  correction.

**Ellmen et al. 2024** (*Sci Rep* 14:20273), "Reconstructing SARS-CoV-2 lineages from mixed
wastewater sequencing data":
- **Aim and inputs.** De novo *lineage definitions*, built from per-sample mutation frequencies
  (coverage ≥ 20 reads).
- **Method.**
  - Missing values are imputed with k-nearest neighbours (k = 5).
  - Non-negative matrix factorization with a user-chosen number of components.
  - Loadings ≥ 0.25 after max-normalization define a lineage.
  - The lineages are then classified with Pangolin, Scorpio, or Nextclade.
- **What it lacks.** No sequencing-noise model, no depth weighting, no time model, and no read
  linkage.
- **Validation.** Qualitative only: a tree and Nextclade calls. No precision or recall is
  reported.
- **Stated limitations.** Manual choice of the number of components; imputation "likely
  introduces error"; sub-lineages are "notoriously difficult to distinguish in wastewater".

**What transfers, and what does not** (detail in `docs/SPRINT5_PHASE3_PROGRESS.md` §0.4):

| Transferred (adapted) | Not transferred |
|---|---|
| Hierarchical agglomeration: complete linkage, valid for non-Euclidean, overlap-based distances | Diploid copy-number expectations and VAF = CP/2 |
| A depth-dependent VAF variance, replaced by the Phase 2 beta-binomial p(1−p)/n·[1+(n−1)ρ] | Tumor purity |
| Adaptive top-down acceptance against expected measurement variation, re-derived as a χ² test | Cellular prevalence |
| Multi-sample profiles, compared only over shared known samples | Subclone and evolutionary interpretation; CNV correction |
| Simulation-based validation, with our own criteria and not MAGOS's J ≥ 0.80 | Lineage assignment and NMF imputation from Ellmen et al. |

**How Ellmen et al. differs from this project:**
- They name lineages; we only group co-varying mutations.
- They impute missing values; we never impute.
- They have no noise model; we use an explicit binomial or beta-binomial model.
- Their number of groups is user-set; ours comes from a stated criterion.

## 3. Exact mathematical definitions

Notation:
- f_is is the observed VAF of mutation i in observation unit s; n_is is its depth.
- A unit's variance factor is c_is = [n + ρ_s(Σ_l n_l² − n)]/n², so Var(f) = p(1−p)·c. For one
  library this is the Phase 2 beta-binomial; with ρ = 0 it is binomial, c = 1/n.
- S_ij is the set of units where both mutations are known, and k_ij = |S_ij|.
- m is the minimum number of shared units (default 2; §8).

**Ordinary clustering (baseline; frequencies only).**

- d_ij = sqrt( (1/k_ij) Σ_{s∈S_ij} (f_is − f_js)² ) if k_ij ≥ m and both mutations are on the same
  reference. Otherwise d_ij = ∞, never 0.
- Complete linkage: D(A,B) = max_{i∈A, j∈B} d_ij. A cluster can never contain a non-comparable
  pair.
- The groups are the clusters after every merge at height ≤ h = 0.10. Every pair in a group
  differs by at most 0.10 RMS.

**Depth-aware grouping (MAGOS-inspired; read counts required).**

- **Pair statistic.** In each shared unit, with the inverse-variance pooled frequency
  p_s = (f_i/c_i + f_j/c_j)/(1/c_i + 1/c_j):

  z_s² = (f_is − f_js)² / [p_s(1−p_s)(c_is + c_js)]

  Units with p_s ∈ {0, 1} are dropped; there both values are equal, and dropping them can only
  make the test stricter. X²_ij = Σ_s z_s², with df_ij = the number of units kept. At ρ = 0 this
  equals the Pearson 2×2 χ² summed over units (tested exactly).
- **Distance.** d_ij = X²_ij / max(df_ij, 1), which is about 1 under "same frequency".
- **Power gate.** λ_ij = Σ_s δ² / [(c_is + c_js) p*_s(1−p*_s)], where δ = 0.10 and p*_s is the
  value within δ/2 of p_s with the largest variance. The pair is comparable only if
  λ_ij ≥ λ₈₀(k_ij), the noncentrality giving 80% power for a level-0.01 χ² test on k_ij degrees
  of freedom. Otherwise d_ij = ∞.
- **Group homogeneity test** (the design-effect-corrected, Rao–Scott-type Pearson test of equal
  proportions):

  X²_G = Σ_s Σ_{i∈G, known in s} (f_is − p_Gs)² / [p_Gs(1−p_Gs) c_is]

  with df = Σ_s (m_Gs − 1) over units where at least 2 members are known and p_Gs ∈ (0, 1).
  P = P(χ²_df > X²_G).
- **Informative-pattern requirement** (attempt 2). Range_G = max_s p_Gs − min_s p_Gs over units
  with at least 2 known members, and it must be ≥ 0.20.
- **Adaptive partitioning.** Walk each complete-linkage tree from its root. Accept a node if
  P ≥ 0.01 and Range_G ≥ 0.20; otherwise split it into its children. Leaves are ungrouped.

**Agreement metrics.**
- ARI is Hubert–Arabie, and NMI uses the arithmetic-mean normalization. Both match scikit-learn
  on test cases.
- Pairwise precision and recall are computed over mutation pairs.
- The false-grouping share is the share of truly ungroupable mutations placed in any group.
- Stability is the share of leave-one-unit-out reruns that keep a group's members together.

**Chi-square functions without SciPy.** These are the regularized incomplete gamma by series and
continued fraction, a Poisson-mixture noncentral χ², and bisection inverses. They match SciPy to
about 1e-11 (tested when SciPy is available).

## 4. Input eligibility and missing-data rules

- **Mutation identity.**
  - Reference_ID (which is one segment) + position + Ref + Alt.
  - Mutations on different references or segments are never compared (distance ∞). This is
    tested at unit level and through the website.
- **Sample-level states** come from the Sprint 3 observation matrix and are never inferred:
  - detected;
  - supported non-detection (A/C/G/T counts show 0 reads at adequate depth);
  - counted but not called;
  - coverage unknown;
  - low or no coverage;
  - covered but value unknown (no base counts);
  - rejected by QC.

  Only the first three give values. A missing variant row is never 0.
- **Count basis per unit.**
  - *Measured*: supplied or A/C/G/T counts.
  - *Reconstructed (approximate)*: frequency × depth, with the Phase 1 rule.
  - *Frequency only*: a reported VAF without usable counts. It is used by ordinary clustering and
    never by the depth-aware method.
  - Reconstructed counts are labelled approximate everywhere.
- **Eligibility.**
  - Ordinary clustering needs at least m known units.
  - The depth-aware method needs at least m units with counts.
  - A reference with Ref-base mismatches is excluded.
  - Unverified references are allowed but flagged, as in Phase 1.
  - Every exclusion has a written reason in `sprint5_mutation_eligibility.csv`.
- **Detection-only input** (no A/C/G/T counts) never receives groups (§9).
- **Size limit.** The pairwise step is O(V²·U). If V²·U > 5×10⁷ (`PAIR_BUDGET`), or V > 1,000, only the mutations
  with the most known values are grouped. This is disclosed on the page, and the rest are
  listed as "not analyzed", not as ungrouped.

## 5. Technical-replicate validation

**Finding.** A shared `Specimen_ID`, date, and location do not show that two records are
separately prepared libraries. They could be one library listed twice, a re-run, or duplicate
rows.

**Rule now used** (`genomics/overdispersion.replicate_pairs`). Exactly two samples share a
`Specimen_ID`, with the same date, location, and reference, and in addition:

1. each has a documented `Library_ID`;
2. the two `Library_ID`s differ;
3. the records are not exact duplicates (identical alleles, alternate counts, and depths).

Failing 1, 2, or 3 gives the new status **"Not estimated: replicate libraries not documented
(unverified)"**, with the reason, and no ρ. Contradicting metadata (dates, locations, references)
still gives "not a usable replicate pair".

`Library_ID` is read from the SNV table or the pooled metadata, one value per sample, and is
detected by aliases such as `library_id` and `replicate_id`.

**Backward compatibility.**
- Phase 2-style files without `Library_ID` still run. Their pairs are now *unverified* and keep
  the binomial baseline; this is the intended downgrade.
- The synthetic replicate demo now includes `Library_ID`. Apart from that one column the files
  are identical, so its Phase 2 calibration results are unchanged: ρ̂ = 0.0105 for SPEC_A and
  0.0021 for SPEC_B.

**Regression tests** (`tests/test_replicate_verification.py`, 10 tests):
- documented pairs are usable;
- no or blank `Library_ID` → unverified;
- the same `Library_ID` → unverified;
- duplicated records → unverified;
- dates still checked first;
- an unverified pair gets no ρ and an explanation;
- a duplicated library file is not calibrated;
- loader behaviour (table, metadata, conflict);
- a Phase 2-style file is downgraded.

**In the profiles.** Only verified pairs are combined into one unit. Their counts are summed, and
the variance uses Σn_l², which is exact for two independent libraries with the same ρ. Unverified
samples at the same site and date stay separate units and are flagged `Same_Site_Date_Units`.

## 6. Ordinary clustering baseline

Defined in §3. It was validated:

- numerically: the distance over shared units uses no zero filling and gives ∞ for insufficient
  overlap; complete-linkage heights match SciPy exactly;
- on deterministic synthetic profiles: three planted patterns at 500× give ARI 1.0.

Known weaknesses, measured in §9–10:
- it groups unrelated mutations with similar levels (look-alike controls: up to 100%);
- it joins independent trajectories (negative controls: 10–67%);
- detection-only input makes this worse.

Sensitivity of the baseline to its own threshold:
- h = 0.05 halves the false grouping but loses most recall (ARI 0.40 at 100×).
- h = 0.15 raises false grouping to 0.46–0.80.

No choice of h gives ordinary clustering both low false grouping and good recall.

## 7. Depth-aware method and its development

Defined in §3. The development history is recorded in full in the progress file.

**Attempt 1** had no informative-pattern requirement. It **failed P1** in the pre-registered
selection of the minimum overlap, on the development seeds:
- look-alike controls at 300×: 58–81% of mutations grouped;
- negative controls at 300×: 4–39% grouped, depending on m and missingness.

Diagnosis: "the homogeneity test does not reject" is absence of evidence. The power gate only
guarantees that a 0.10 difference is detectable. Flat profiles at similar levels agree in
*level*, not in *pattern*.

**Attempt 2** adds the informative-pattern requirement: the pooled profile must span ≥ 0.20, two
resolution steps. This is a stated heuristic, not a likelihood result, and it costs recall for
groups with flat frequencies.

**Ablation.** The same requirement was added to the ordinary method. It still grouped 20–35% of
negative-control mutations, so the depth-aware method's control behaviour is not explained by
the pattern requirement alone.

The main grid was run **once**, on seed streams not used during development. Further tuning on
it was deliberately avoided, to prevent overfitting to the simulation.

**Safe behaviour.**
- The method returns no groups, and says why, when:
  - no pair is comparable (low depth or overlap);
  - every candidate fails the homogeneity or pattern test;
  - the input is detection-only.
- A larger assumed ρ makes pairs less precise. It can join a pair that differs by less than the
  added noise, and when the noise prevents resolving 0.10 it makes pairs non-comparable
  (abstain). Both behaviours are unit-tested.

## 8. Simulation design and predefined criteria

The design was written to the progress file before running. Seed 20261115; 30 datasets per cell.

- **Model.**
  - 2 sites × 6 dates, giving 12 units.
  - 3 latent sources with logistic shares.
  - Planted groups of 10, 6, and 3 mutations (unequal sizes).
  - 2 "shared" mutations carrying the sum of two sources (partly overlapping patterns; truth:
    ungroupable).
  - 8 background mutations with their own trajectories.
- **Controls.**
  - Negative: 29 independent trajectories.
  - Look-alike: 14 pairs of unrelated mutations at the same level, L + N(0, 0.03).
- **Depth.** 20, 50, 100, 300, or mixed, each with unit and cell variation.
- **Reads.** Beta-binomial with true ρ 0 or 0.02.
- **Missing coverage.** 0, 30%, or 60%.
- **Detection-only variant.** A caller threshold of 2 reads and 2%.
- **Technical replicates.** Two libraries per unit, with ρ estimated by the Phase 2 estimator.
- **Units and size.** 3–24 units, and 29 or 87 mutations.
- **Null calibration.** Two mutations with the same frequency.
- **Runtime and memory.**

**Methods compared:**
- ordinary (h = 0.05, 0.10, 0.15);
- the ordinary + pattern ablation;
- depth-aware with assumed ρ ∈ {0, 0.002, 0.02, 0.05}, so one is correct and the others are wrong;
- depth-aware with ρ calibrated from replicates;
- no grouping.

**Metrics.** ARI, NMI, pairwise precision and recall, false-grouping share, pair false-positive
rate, stability, excluded count, abstention, and runtime. Unassigned and excluded mutations count
as singletons. Each metric is reported as the mean over seeds with a 95% normal interval (in the
CSV files).

**Criteria:**
- **P1.** In negative and look-alike controls, at most 10% of eligible mutations are grouped, with
  the correct noise model.
- **P2.** ARI ≥ 0.80 for planted structure at ≥ 100× with ≤ 30% missing, with the correct model.
- **P3.** "Improves" only if the paired ARI difference has a 95% interval above 0.
- **P4.** Misspecification is reported.
- **P5.** The method abstains with a reason when data are insufficient.
- **P6.** The method is offered only in the conditions where P1 and P2 hold.

**Minimum shared units.** Chosen by the pre-registered rule: the smallest of {2, 3, 4, 6} meeting
P1 in the development control cells. For attempt 2 every m met it, so **m = 2**. It can be
raised on the website (2–20) but not lowered.

## 9. Group recovery and false-positive results

All results in this section are synthetic. The files are in `docs/sprint5_phase3_results/`
(`attempt2_*`); means are over 30 seeds.

**False grouping in the controls, with the correct noise model** (P1, limit 0.10):

| True ρ | Control | Depth-aware, range over depth × missingness | Cells failing P1 | Ordinary h = 0.10 |
|---|---|---|---|---|
| 0 | Negative | 0.000–0.201 | 5 of 15: 50×/0%, 100×/30%, 300×/60%, mixed/0%, mixed/30% | 0.13–0.67 |
| 0 | Look-alike | 0.000–0.243 | 1 of 15 (mixed/0%) | 0.35–1.00 |
| 0.02 | Negative | 0.000–0.059 | 0 | 0.10–0.66 |
| 0.02 | Look-alike | 0.000 | 0 | 0.28–0.87 |

**P1 result.** It holds whenever overdispersion is present and modelled, but this is mostly
because the method then abstains (see P2). It fails in 6 of 30 binomial-truth control cells.

**Planted structure, ARI** (P2, target ≥ 0.80 at ≥ 100× with ≤ 30% missing):

| Depth / missing | Ordinary | Depth-aware, ρ = 0 truth | Depth-aware, ρ = 0.02 truth (assumed 0.02) |
|---|---|---|---|
| 100× / 0% | 0.94 | **0.91** (precision 0.97, recall 0.90) | 0.03 (abstains) |
| 100× / 30% | 0.79 | 0.40 | 0.00 |
| 300× / 0% | 0.94 | **0.91** (precision 0.98) | 0.26 |
| 300× / 30% | 0.92 | **0.81** (precision 0.98) | 0.01 |
| 50× / 0% | 0.75 | 0.49 | 0.00 |
| 20× (any missingness) | 0.16–0.34 | 0.00 (abstains) | 0.00 |
| mixed / 0% | 0.64 | **0.82** | 0.04 |

**P2 result.** It is met in 3 of 4 binomial-truth cells (100×/0%, 300×/0%, 300×/30%) and in none
with overdispersion. False grouping of shared and background mutations in planted data was
0–12% for the depth-aware method, against 4–56% for ordinary clustering.

**Calibrated mode** (two libraries per unit, ρ estimated from them):
- Calibration is impossible at 100× per library, because sites must be above 100× in both
  replicates.
- At 300× and 1,000×, ρ was estimable for all units in 67–100% of datasets.
- Planted ARI with calibrated ρ: 0.80/0.89 at 300×/1,000× when ρ = 0.005 and there is no
  missingness; 0.10–0.60 when ρ = 0.02.
- Negative-control false grouping with calibrated ρ: 0.04–0.23. **It fails P1 in 4 of the 8 cells
  where calibration was possible.**
  Combining two libraries gives enough power to compare pairs, while ρ widens what counts as
  similar.

**More mutations means more false grouping.** At 12 units, 100×, and 30% missing,
negative-control false grouping rose from 0.19 with 29 mutations to 0.33 with 87. At 24 units it
rose from 0.04 to 0.10. Coincidental similarity grows with the number of pairs compared. This
matters for real datasets with hundreds of mutations.

**Detection-only input** (no supported zeros):
- ordinary clustering grouped 18–100% of unrelated mutations;
- the depth-aware method grouped 0–18% with the correct noise model, failing P1 at 100× and mixed
  depth.

Hence detection-only input never receives groups.

**Null calibration of the homogeneity test** (2,000 replicates per cell, nominal α = 0.01):
- type I error 0.000–0.0155;
- conservative at low frequency and low depth;
- 0.0155 at 50×, p = 0.5, ρ = 0.02, about 2 Monte Carlo standard errors above nominal.

**Stability** (leave one unit out, planted, 100–300×). Depth-aware groups: 0.57–1.00.
Ordinary: 0.72–0.97.

## 10. Comparison between methods

| Question | Answer (synthetic) |
|---|---|
| When does depth-aware grouping improve accuracy? | By ARI, only at mixed depth with no missing coverage (+0.18, 95% CI 0.07 to 0.30). It always improves **precision and false grouping**, often 3–10-fold. |
| When is it similar? | At 100× and 300× with no missing coverage and ρ = 0 (ARI difference −0.03; CI includes 0). |
| When is it worse? | By ARI, in 27 of 30 count-based conditions: low depth, missing coverage, and any overdispersion, where it abstains. |
| When does missing coverage prevent reliable grouping? | At 60% missingness every method's ARI is ≤ 0.29. With 30% missingness the depth-aware method keeps its precision but loses recall. Detection-only profiles mislead both methods. |
| How do wrong ρ assumptions behave? (P4) | **Binomial assumed when the truth is ρ = 0.02:** false grouping 0–14% (negative) and 0–28% (look-alike); recall collapses at depth. **Assumed ρ too small (0.002 vs 0.02):** similar, up to 23% false grouping in look-alike controls. **Assumed ρ too large (0.05):** never groups, because the power gate makes every pair non-comparable. A hypothetical ρ is therefore a sensitivity scenario, never a calibration. |
| Ablation | The pattern requirement alone (ordinary + pattern) still grouped 8–39% of negative-control mutations. Depth-awareness is responsible for most of the reduction. |
| Robust enough for exploratory use? | **Not by default** (P6). It is offered only as an opt-in experimental view with the measured error rates shown. |

## 11. Real workbook eligibility and findings (`mmc2.xlsx`, real)

| Quantity | Value |
|---|---|
| Distinct mutations | 5,669 (NC_045512.2; reference unverified because no FASTA was supplied) |
| Samples, locations, dates | 52, 23, 29 |
| Observation-matrix cells | 294,788, of which 286,826 (97.3%) have **unknown coverage** |
| Supported non-detections | **0** (no per-position depth or A/C/G/T counts) |
| Count basis | All 7,958 profile values reconstructed (approximate) |
| Mutations known in ≥ 2 / ≥ 3 / ≥ 4 / ≥ 5 / ≥ 10 samples | 1,117 / 374 / 170 / 101 / 24 |
| Pairs among the 374 mutations in ≥ 3 samples | 69,751; median shared samples 0; 4,584 share ≥ 3 |
| Technical replicates | None |

**Result.** `mmc2` is detection-only. In that condition both methods failed P1 in simulation, so
the website **shows no groups**, even when the experimental view is requested. It gives the
detection-only explanation instead.

It still provides:
- the 7,958 profile rows (all "Detected", no invented zeros, no ρ);
- the eligibility table for 5,669 mutations, with 1,117 eligible at m = 2;
- an ungrouped heatmap of the 60 mutations with the most known values.

Every Sprint 3 and Phase 1–2 output is unchanged (§15).

## 12. Website integration

Section 9 now has, after the Phase 2 panel, a panel titled **"Exploratory Mutation-Frequency
Groups"**. It contains:

- the required sentence: "Mutations with similar frequency patterns across samples may be grouped
  together for comparison. These groups do not necessarily represent the same viral genome or
  lineage.";
- a badge: *Profiles only: no groups shown*, *Experimental (not validated)*, or *Off*, followed
  by the reason;
- key statistics:
  - mutations evaluated;
  - mutations eligible for grouping (and with read counts);
  - mutations excluded because of missing information;
  - calibration status;
  - in the experimental view, the number of groups and the method used;
- the one-sentence validation summary;
- the first figure;
- a methods table;
- collapsible group details, comparisons, the remaining figures, and the limitations;
- downloads.

**Figures.** Each is size-limited, and every caption states what was selected:
- `mutation_profile_heatmap` (default; up to 60 mutations, unknown cells hatched and never 0);
- `mutation_group_heatmap`;
- `mutation_group_trajectories` (up to 6 groups);
- `grouping_comparison` (a contingency of ordinary versus depth-aware groups).

**Settings.** Advanced settings → "Exploratory mutation-frequency groups":
- *Profiles only* (default), *Also show experimental groups (not validated)*, or *Off*;
- "Minimum shared samples", from 2 to 20.

Re-runs keep both settings, and invalid values give a message, not an error.

**Downloads:**
- `sprint5_mutation_profiles.csv` (long; coverage state and provenance);
- `sprint5_mutation_eligibility.csv` (state counts and exclusion reasons);
- `sprint5_grouping_parameters.csv` (method, model, thresholds, calibration, mode, validation
  status);
- in the experimental view: `sprint5_mutation_groups.csv`, `sprint5_mutation_group_summary.csv`,
  `sprint5_mutation_group_trajectories.csv`, and `sprint5_grouping_comparison.csv`, which includes
  hypothetical-ρ scenarios labelled "Hypothetical; not calibrated";
- `sprint5_group_individual_overlap.csv`, when individual sequences are also uploaded (§12a).

All Phase 1–2 downloads are kept. Nothing is added to the PCoA section, and PCoA coordinates are
never relabelled as clusters.

**Status and provenance.** A status row "Exploratory mutation-frequency groups (Sprint 5
Phase 3)" reads *Limited* (profiles only, or experimental) or *Skipped*. There are 4 new
provenance parameters.

**12a. Comparison with Sprint 2 individual sequences (Milestone 7).** When both sources are
present, `sprint5_group_individual_overlap.csv` counts, for each group and each individual
sequence, how many of the group's mutations it carries, calls as reference, or cannot call. The
counts come from the existing cross-source evidence, matched by reference-aware `Variant_Key`.
Every row says: "Descriptive only. Carrying some or all of a group's mutations does not place
this sequence in the group, and a single shared mutation is not evidence of a relationship." No
membership, linkage, lineage, or transmission is inferred.

**Synthetic mixture demo** (`examples/integrated/make_mixture_example.py`, seed 20261115):
- 12 specimens × 2 documented libraries, about 300× each, true ρ = 0.005, and full A/C/G/T
  counts;
- `synthetic_mixture_truth.csv` holds the planted groups; it is not a website input.

Verified results:
- Phase 2 calibrated all 12 specimens (ρ̂ 0.0009–0.0097).
- The experimental view runs in **empirically calibrated** mode.
- It recovered planted groups A (10 mutations) and B (6) exactly, with no false members
  (ARI 0.97, precision 1.0). It left the flat group C (3) and all 10 non-group mutations
  ungrouped.
- Ordinary clustering also found A, B, and C, but formed groups containing 4 non-group
  mutations (ARI 0.91, precision 0.86).
- On the replicate demo, which has no group structure, the depth-aware method formed 0 groups
  while ordinary clustering grouped 32 of 47 mutations.

## 13. Files changed

**New files:**

| File | Contents |
|---|---|
| `genomics/mutation_profiles.py` | Units, known observations, combining verified replicates, variance factor, eligibility |
| `genomics/mutation_grouping.py` | Distances, power gate, complete linkage, homogeneity and pattern tests, top-down acceptance, ARI, NMI, precision and recall, stability, size limit, assignment table |
| `genomics/chi_square.py` | χ² and noncentral χ² tail functions, inverses, power noncentrality |
| `genomics/mutation_group_analysis.py` | Website analysis: modes, tables, figures, Sprint 2 overlap, page summary |
| `validation/simulate_mutation_groups.py` | Pre-registered simulation study |
| `examples/integrated/make_mixture_example.py` and `synthetic_mixture_{snvs,base_counts,truth}.csv`, `synthetic_mixture_reference.fasta` | Synthetic demo |
| `tests/test_replicate_verification.py` (10), `tests/test_mutation_grouping.py` (21), `tests/test_mutation_grouping_simulation.py` (8), `tests/test_sprint5_phase3_web.py` (11) | 50 new tests |
| `docs/SPRINT5_PHASE3_PROGRESS.md`, this report, `docs/sprint5_phase3_results/*.csv` | Records and result tables |

**Modified files:**

| File | Change |
|---|---|
| `genomics/overdispersion.py` | `replicate_pairs` needs documented, distinct `Library_ID`s and rejects duplicates; new `UNVERIFIED_PAIR` status |
| `wastewater/processed_snv.py`, `genomics/column_detection.py`, `genomics/pooled_observations.py` | Optional `Library_ID` (aliases, loader, provenance column) |
| `genomics/depth_aware_analysis.py` | Attaches `Library_ID` as the last evidence column; calls the Phase 3 step |
| `genomics/integrated_analysis.py` | Default settings, status row, provenance parameters, Sprint 2 overlap call |
| `genomics/figures.py` | Three new figure functions |
| `genomics/overdispersion_analysis.py`, `integrated_view.py`, `templates/index.html` | Replicate wording names `Library_ID`; Phase 3 panel, view, and settings |
| `templates/_toolbar.html` | Re-runs keep the new settings |
| `app.py` | New settings validated and passed |
| `dashboard_guide.py` | Guides for the four new figures |
| `examples/integrated/make_replicate_example.py` and `synthetic_replicate_snvs.csv` | `Library_ID` column added (draws unchanged) |
| `tests/test_overdispersion.py` | The fixture helper gives each sample a `Library_ID`; no assertion changed |

## 14. Numerical verification

- **χ² functions.** Closed forms (k = 1, 2), a critical value of 3.8415, and requested power hit
  to 1e-6. They match SciPy for the central and noncentral tails to about 1e-11.
- **Pair statistic and group test.** Equal to hand-computed Pearson 2×2 and k×2 χ² at ρ = 0, to
  1e-9.
- **Variance factor.** Equals the Phase 2 design-effect variance for one library and the exact
  two-library combination.
- **Ordinary distance.** Hand computation over shared units; ∞ (not 0) for insufficient overlap
  and for different references.
- **Complete linkage.** Hand example (forest, heights) and identical merge heights to SciPy on
  random data.
- **ARI, NMI, and precision/recall.** Equal the scikit-learn values for the test partitions
  (4/7 and 0.8).
- **Profiles.**
  - Replicate counts summed exactly, with Σn² kept.
  - Single-library Wilson bounds identical to Phase 1.
  - Unknown cells have no row.
  - Frequency-only values are excluded from the depth-aware method.
  - State counts sum to the number of samples per mutation (web test).

## 15. Full regression results

**Command**, from the project folder. MAFFT is on PATH (`C:\Tools\mafft-win`) and
`FONTENELE_MMC2` points to the published workbook (md5 `c0a75b7fa5bde8db09746613354a9f64`,
unchanged):

```bash
python -m unittest discover -s tests -v
```

| Run | Result |
|---|---|
| Baseline at the start of Phase 3 | `Ran 535 tests — OK (skipped=6)`: five mmc2 tests skipped because `FONTENELE_MMC2` was unset. The four mmc2-dependent modules were then rerun with it set: `Ran 51 tests — OK` |
| Intermediate full run (before the memory fixes in §16) | `Ran 585 tests in 378.7 s — OK (skipped=4)` |
| **Final full run** (final code; MAFFT on PATH; `FONTENELE_MMC2` set) | **`Ran 585 tests in 371.8 s — OK (skipped=4)`**: 535 pre-existing and 50 new. The 4 skips are the 3 optional SciPy cross-checks (SciPy is not in `.venv`) and the optional NCBI files test. All real-`mmc2.xlsx` tests and the real-MAFFT test ran |
| SciPy cross-checks (system Python, final code) | `test_mutation_grouping.py`: `Ran 21 tests — OK`; `test_beta_binomial.py`: OK. No skips, so the three `.venv` skips pass where SciPy exists |

**Output regression** (verified, final code). `run_integrated_analysis.py` was rerun and
compared file by file with the saved Phase 2 outputs:
- **Synthetic integrated example** (`--example`, with MAFFT): 58 of 61 files are identical,
  including every figure (byte-identical PNGs).
- **`mmc2.xlsx`**: 35 of 38 files are identical, including every figure.

The only differences in both are intended:
- `analysis_status.csv`: one Phase 2 detail now names `Library_ID`, and one Phase 3 row was added;
- `provenance.csv`: 4 Phase 3 parameters were added;
- `sprint5_vaf_observations.csv`: all original columns are identical, with `Library_ID` appended
  last;
- the new Phase 3 files.

This covers MAFFT alignment, Sprint 2 distances and screening, Sprint 3 QC, Yue–Clayton, PCoA,
cross-source matching, Phase 1 VAF and temporal trends, and Phase 2 calibration and sensitivity.
Sprint 4 literature code was not modified.

## 16. Performance and memory measurements

**Measured locally.** Windows peak working set, project `.venv`, `OPENBLAS_NUM_THREADS=1`, the
guided web path `/review` → `/run` on `mmc2.xlsx`, with the same profiling scripts as Phase 2.
Render was **not** measured.

| Run (`mmc2.xlsx`) | Phase 2 (reported) | Phase 3 off (same code, A/B) | **Phase 3 default (profiles)** | Phase 3 experimental |
|---|---|---|---|---|
| One run, peak | 296 MB | 300 / 299 MB | **300 / 297 MB** | 300 MB |
| One run, time | 10.4 s | 10.3 / 10.4 s | 11.7 / 11.7 s | 11.5 s |
| Four runs in one process, peak | 306 MB | 322 MB | 325 / 318 MB | — |
| 512 MB cap and 10% of one CPU core | 285 MB, 110.4 s | — | **288 MB, 118.8 s** (job committed 290 MB) | — |
| Largest download | 6.27 MB | — | 6.22 MB (`sprint5_vaf_observations.csv`) | — |

**A memory regression was found and fixed during this milestone.**

- **First measurement:** the default mode peaked at 343–350 MB, about 50 MB above "off". Tracing
  each step (`tracemalloc` inside the real run) located the cause:
  - the profile heatmap drawn with `imshow` added +64 MB, because `imshow` resamples into
    full-canvas floating-point buffers;
  - `mutation_eligibility` made string copies of the 295,000-row matrix and built a wide
    integer frame, adding +13 MB.
- **Fixes:**
  - the heatmap is drawn with `pcolormesh` (0.6 MB in isolation; the same image);
  - the eligibility counts use integer codes with `np.bincount`. The output table is identical
    to the earlier version (checked on `mmc2`).
- **After the fixes**, Phase 3's traced peak is +16 MB, and the process peak equals "off" within
  run-to-run variation (297–300 MB, against 299–300 MB).

**Interpretation.**

- **Time.** About +1.3 s on `mmc2`, mostly drawing one extra figure and building the profile
  tables. The constrained run was 8.4 s slower than the Phase 2 report (118.8 s against
  110.4 s), still well inside gunicorn's 180 s timeout.
- **Memory over four runs.** It rises with the bounded download cache
  (`RESULT_DOWNLOAD_LIMIT = 4`), in both modes: 322 MB with Phase 3 off and 318–325 MB on. Phase 3
  adds no detectable growth. Today's machine gives higher four-run peaks than the Phase 2 report
  did, for the identical "off" configuration.
- **Experimental mode on `mmc2`** does no pairwise work, because the input is detection-only.
- **Pairwise cost, from the simulation runtime study** (`attempt2_grouping_runtime.csv`):

  | Mutations × units | Ordinary | Depth-aware | Traced memory (depth-aware) |
  |---|---|---|---|
  | 100 × 12 | 0.001 s | 0.02 s | 0.8 MB |
  | 300 × 52 | 0.003 s | 0.9 s | 7 MB |
  | 1,000 × 52 | 0.03 s | 1.3 s | 76 MB |
  | 1,000 × 192 | 0.04 s | 22.7 s | 76 MB |

  The last row is why `PAIR_BUDGET` (V²·U ≤ 5×10⁷) limits large inputs to a disclosed subset.
  1,000 × 192 would be cut to about 510 mutations.
- **Not measured:** Linux and Render memory, and Render's proxy timeout.

## 17. Known limitations

1. **The method did not meet its pre-registered criteria in every condition (§9).** It is
   therefore not offered by default. The experimental view is for exploration only, and its
   groups can include unrelated mutations, especially:
   - with many mutations;
   - with mixed depths;
   - with missing coverage;
   - with calibrated ρ when there is a lot of data.
2. **Synthetic validation only.** No real pooled dataset with known group structure was
   available. The simulation is one plausible model (logistic source shares, independent cells,
   beta-binomial reads); real data have amplicon dropout, primer bias, and correlated sites.
3. **The informative-pattern threshold (0.20) is a heuristic** introduced after attempt 1 failed.
   It is not derived from a likelihood.
4. **Recall is low** whenever the measurement noise is large compared with the 0.10 resolution
   (low depth, overdispersion). The method abstains rather than guesses.
5. **Detection-only data cannot be grouped**, and this includes the published `mmc2.xlsx`.
6. **False grouping grows with the number of mutations compared.** No multiple-comparison
   correction across the whole tree is applied.
7. **No read-level linkage, no lineage definitions, and no phylogeny.** A group is never a
   genome, haplotype, lineage, outbreak, or transmission cluster.
8. **Stability uses leave-one-unit-out reruns**, not a parametric bootstrap. It is skipped above
   a size budget.
9. **The Phase 1–2 limitations still apply.** Reconstructed counts are approximate, and ρ is per
   library pair.

## 18. How to test locally

From the project folder:

1. Regenerate the synthetic demos (deterministic; the files are already present):

   ```bash
   python examples/integrated/make_mixture_example.py
   ```

   ```bash
   python examples/integrated/make_replicate_example.py
   ```

2. Run the simulation study (about 5 minutes):

   ```bash
   python validation/simulate_mutation_groups.py validation_output
   ```

3. Run the Phase 3 tests:

   ```bash
   python -m unittest tests.test_replicate_verification tests.test_mutation_grouping tests.test_mutation_grouping_simulation tests.test_sprint5_phase3_web -v
   ```

4. In the website:
   - Choose **Pooled sequencing**.
   - Upload `synthetic_mixture_snvs.csv`, `synthetic_mixture_base_counts.csv`, and
     `synthetic_mixture_reference.fasta`, and run.
   - In section 9 the panel says *Profiles only: no groups shown*.
   - Then open **Advanced settings → Exploratory mutation-frequency groups → Also show
     experimental groups** and re-run. The panel shows calibrated experimental groups G1 and G2;
     compare them with `synthetic_mixture_truth.csv`.
   - Upload `mmc2.xlsx` with the experimental view on to see the detection-only explanation.
5. **To document technical replicates:** give each library its own `Sample_ID` and `Library_ID`,
   and give both libraries of one specimen the same `Specimen_ID`, date, and location.

## 19. Readiness for biological mutation annotation

Not ready for biological interpretation of groups. What exists and what is needed:

- **Ready.**
  - Profiles keep reference, segment, position, and alleles, plus `Annotated_Region` when the
    input supplies one.
  - A later gene/protein annotation step could attach consequences per `Variant_Key` without
    changing the grouping.
- **Needed before annotating groups.**
  1. Real pooled data with known composition: spike-in mixtures or replicate libraries from the
     same protocol, to test false grouping outside simulation.
  2. Calibrated ρ for the specimens being grouped.
  3. A correction for the number of mutations compared.
  4. Ideally, read-level co-occurrence for nearby mutations, the only direct evidence that they
     sit on one molecule.
- **Out of scope.** Annotating a group with lineage-defining mutations must not be presented as
  lineage identification.

## 20. Render deployment considerations

- **No new dependency.** `requirements.txt` and the Dockerfile are unchanged. SciPy appears only
  in optional cross-check tests, which are skipped without it.
- **Default cost.** Only profiles, eligibility, and one heatmap are computed: +1.3 s and no
  detectable memory on `mmc2` (§16).
- **Experimental cost.**
  - The pairwise step is O(V²·U) and is capped at V²·U ≤ 5×10⁷ (with disclosure).
  - At V = 1,000 and U = 52 it took 1.3 s and 76 MB traced in the simulation runtime study.
  - Stability reruns are skipped above V²·U² > 2×10⁸.
  - Detection-only input such as `mmc2` does no pairwise work.
- **Gunicorn timeout.** Unchanged at 180 s.
- **To verify on Render after review** (not done here):
  1. Upload `mmc2.xlsx` and confirm section 9 shows both panels and the downloads.
  2. Upload the mixture demo with the experimental view and confirm the calibrated groups.
  3. Watch the memory graph and the request time.
- **Nothing was pushed or deployed**, so Render remains **unverified**.

## Required completion conditions

| Condition | Status | Evidence |
|---|---|---|
| Existing Phase 1–2 functionality preserved | **Met** | Output comparison (§15); all pre-existing tests pass |
| Input and calibration eligibility rules validated | **Met** | §4, §5; unit and web tests |
| Technical replicates not inferred from insufficient metadata | **Met** | `Library_ID` rule; 10 regression tests (§5) |
| Mutation profiles preserve coverage and provenance | **Met** | §4; profile and web tests (unknown never 0; count basis kept) |
| Conventional baseline implemented and tested | **Met** | §6, §14 |
| Uncertainty-aware method implemented and validated, or documented as unsupported | **Met (documented as not meeting its criteria; experimental only)** | §7–§10; attempt 1's failure and attempt 2's partial failure are reported |
| Synthetic comparisons use known ground truth | **Met** | §8 |
| False-grouping and missing-data behaviour measured | **Met** | §9 |
| Real data get no fabricated zeros or calibrated parameters | **Met** | §11; web test on `mmc2` |
| The website separates exploratory groups from lineages | **Met** | §12; required sentence, labels, and badges |
| Unsupported datasets continue through the original analyses | **Met** | Web tests: detection-only, frequency-only, single date, two references, off mode, invalid setting, `mmc2` |
| Regression tests pass | **Met** | Ran 585 tests — OK (skipped=4, all optional); outputs identical apart from intended additions (§15) |
| Performance evaluated | **Met** | §16: one-run peak equal to Phase 3 off (297–300 MB); a +50 MB regression was found and fixed; 512 MB / 10% CPU: 288 MB, 118.8 s |
| The report accurately records what was and was not achieved | **Met** | This document, including negative results |
| No commits, pushes, deployments, or paid API calls | **Met** | `git status` shows uncommitted changes only; AI disabled in every run |

Every mandatory condition is met. The uncertainty-aware method is implemented, numerically verified, and validated on synthetic data with known truth. It **did not meet** its pre-registered false-grouping and recovery criteria in every condition (attempt 1 failed outright; attempt 2 failed in part). It is therefore documented as not validated for general use, kept as an opt-in experimental view, and never shown by default, as the goal requires for an algorithm that fails validation. Real data receive no groups, invented zeros, or calibrated parameters. Render deployment remains unverified because deploying was out of scope.

**SPRINT 5 PHASE 3 COMPLETE**
