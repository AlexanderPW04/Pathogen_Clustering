# Independent audit checklist: complete application after Sprint 5 Phase 4A

This checklist is for a later, **independent** audit of the whole application (Sprints 2–5). It
lists what to check and how. It does not record that any check passed.

Passing unit tests shows only that the code runs as its tests expect. It does **not** show that a
method is scientifically valid. Mark each item with the kind of evidence actually obtained:

| Code | Kind of evidence | What it shows |
|---|---|---|
| **C** | Code execution correctness | The code runs and does what its tests expect |
| **M** | Mathematical correctness | Formulas reproduce hand or independent calculations |
| **S** | Scientific method validity | The method is appropriate and calibrated for pathogen pooled sequencing (simulation with known truth, literature, expert review) |
| **R** | Real-data compatibility | Published or real files import and analyze correctly, checked against independent values |
| **H** | Hosted website reliability | Works on the deployed Render instance under its memory and CPU limits |

---

## 1. Input-format compatibility (C, R)

- [ ] **Each documented format imports end to end.** Fontenele Table 2, iVar TSV, VCF (AD/DP and
  AF), standard long CSV/TSV/XLSX, wide sample-by-mutation matrices, base-count (A/C/G/T + Ref)
  tables, per-position tables, GenBank/GFF3.
- [ ] **Raw formats are refused with a clear message:** FASTQ, BAM, CRAM, unknown files.
- [ ] **Header-row detection and override.**
  - Tables with title rows (Fontenele, Schumann).
  - A wrong detected header, overridden by the user.
- [ ] **Worksheet choice.** Workbooks with several data sheets ask; nothing is guessed.
- [ ] **Aliases.** For every field in `genomics/column_detection.py`:
  - each alias is a documented name for that field;
  - "Site", "sequencing", inStrain `var_freq` / `var_base` are **not** mapped.
- [ ] **Wide-matrix parsing.**
  - Every accepted header pattern parses to the right Ref/Position/Alt.
  - Insertions, deletions and multi-base headers are listed, never analyzed.
  - NA cells never become 0.
- [ ] **Base-count derivation.**
  - The minimum read count comes from the source's methods.
  - Derived records are compared with the study's own calls (Crits-Christoph: 293/293 rows).
  - The frequency is labelled "calculated".
- [ ] **Percent-to-fraction conversion is exact.** Decimal arithmetic; significant digits kept for
  the Phase 1 count-reconstruction rule.
- [ ] **0-based to 1-based conversion** only after the user's choice; the warning is checked
  against a reference FASTA.
- [ ] **Mapping profiles.**
  - They contain no records.
  - They apply only to matching schemas.
  - Incompatible files are reported.
  - Note the known limitation: two FASTA roles are matched by order.

## 2. Metadata matching (C, R)

- [ ] **Joins use identical Sample_ID text only.** Near matches are listed and never joined.
- [ ] **Unmatched samples** block the run until the user accepts them. Their values stay unknown
  (Crits-Christoph S4 against S1: 0 of 7).
- [ ] **Conflict resolutions** behave exactly as labelled: metadata / variant table / unresolved.
  Duplicate-ID conflicts accept "unresolved" only.
- [ ] **Rows without IDs** (summary rows) are reported, not joined.
- [ ] **Specimens and libraries.**
  - Specimens are never merged with libraries.
  - Replicate candidates need distinct documented Library_IDs (Phase 3 rules).
- [ ] **Dates.** Ambiguous day/month dates are never guessed. Date-time cells keep their original
  text.

## 3. Scientific calculations (M, S)

- [ ] **Reported versus calculated VAF.** Reported VAF is never overwritten. Calculated VAF is
  alternative reads / total depth.
- [ ] **Wilson intervals:** recompute by hand for a sample of records.
- [ ] **Beta-binomial.**
  - Moment estimator, ρ bounds, and calibration only from verified replicates (Phase 2 report).
  - Re-run `validation/simulate_overdispersion.py`.
- [ ] **Depth QC thresholds** (Fontenele presets) are labelled as that study's settings.
- [ ] **Yue–Clayton / PCoA.**
  - Full-count mode only with per-position A/C/G/T.
  - Otherwise the plot is labelled an exploratory Jaccard ordination.

## 4. Reference verification (C, M, R)

- [ ] **Labels are never assumed equal.** A stated reference is recorded as user-stated. Every Ref
  base is checked when a FASTA is uploaded.
- [ ] **WIV04 versus NC_045512.2:** treated as different sequences (29,840 against 29,903 bases).
- [ ] **Segmented genomes.** Each segment is its own reference; segments are never compared with
  each other.

## 5. VAF and sequencing depth (C, M, R)

- [ ] **Missing depth stays missing.** Frequency-only tables are labelled "Frequency only (no
  depth)" (Schumann: 5,145 of 5,145).
- [ ] **Count provenance:**
  - measured (supplied or base counts);
  - reconstructed (approximate; Fontenele: 7,973);
  - frequency only.
- [ ] **Depth definitions.** Documented formats carry their definitions. User-stated definitions
  are provenance only.

## 6. Beta-binomial calibration (S)

- [ ] **Calibration needs documented technical replicates.** It is never inferred from shared
  dates or sites.
- [ ] **Hypothetical ρ scenarios** are labelled "scenarios, not estimates".

## 7. Missing coverage (C, S)

- [ ] **An unreported variant is "coverage unknown"** without a per-position table, never 0%. Check
  `pooled_observation_matrix.csv`.
- [ ] **`canonical_coverage.csv`** is empty without a per-position table. It is never zero-filled.

## 8. Experimental clustering (S)

- [ ] **Off by default.** `DEFAULT_SETTINGS["grouping_mode"] == "profiles"`, and readiness shows
  "Experimental".
- [ ] **When enabled**, groups are labelled "not validated". The Phase 3 error rates are shown.
- [ ] **Groups are never called lineages or haplotypes.**

## 9. Temporal trends (C, S)

- [ ] **Only dated samples are used.** Detections-only time series are labelled as such.
- [ ] **Single-date and undated data are reported as unavailable**, not drawn.

## 10. Download consistency (C)

- [ ] **The XLSX sheets equal the CSV downloads** value for value (tests do this for Observations,
  Samples and Readiness; extend it to every sheet).
- [ ] **The PDF** findings, readiness and limitations equal the dashboard text.
- [ ] **Size limits.** Files over 8 MB are listed as skipped. The XLSX row limit is reported, not
  silently cut.

## 11. Browser usability (H)

- [ ] **A new user** can go upload → review → Step 3 → run → results without help.
  Use Schumann mmc1 S2.1 and Fontenele mmc2.
- [ ] **Primary buttons stay visible** (the sticky action bar). Technical details sit in collapsed
  sections.
- [ ] **Accessibility:** labels on every input; the Step 3 statuses are readable without color.
- [ ] **Check in a real browser:** the in-app browser preview was refused during Phase 4A.

## 12. Performance on Render (H)

- [ ] **Deploy only after review.** Upload `mmc2.xlsx` twice and Schumann mmc1 S2.1 once; record
  peak memory in Render Metrics.
- [ ] **`MAX_OBSERVATION_CELLS=700000`** (Dockerfile) refuses Schumann Berlin35 (8,862 × 109 ≈
  966,000 cells) with an explanation. Confirm that the message is clear.
- [ ] **Uploads over the 10 MB request limit** (Schumann mmc3, 14 MB) get the 413 page.

## 13. Multiple real datasets (R)

- [ ] **Re-run `validation/check_published_datasets.py`** with all three datasets and confirm 0
  FAIL.
- [ ] **Add at least one non-SARS-CoV-2 dataset**, for example influenza or mpox wastewater, before
  claiming the import is general for pathogens.

## 14. Readiness for future biological annotation (S)

- [ ] **Phase 4 annotation** against independent labels: Fontenele codons (5,613/5,613) and
  Schumann VEP labels (1,830/1,830).
- [ ] **Overlapping CDS, frameshifts, and segmented genomes** (Phase 4 report §4–6).
- [ ] **Shared-codon SNVs** are never combined, because pooled reads do not show linkage.
