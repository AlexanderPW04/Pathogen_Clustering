# Final Polish: Implementation Report

**Date:** 2026-11-16
**Working branch:** `fix/final-polish-20261115`. Local commits only: nothing pushed, merged, or
deployed.
**Resume record:** [FINAL_POLISH_PROGRESS.md](FINAL_POLISH_PROGRESS.md)

**What was not done:**
- No paid API calls.
- `.env` was not read or changed, and no secrets are in any change.
- Experimental mutation grouping stays disabled by default.
- Live Claude NLP is **not** validated.
- Pathogen support is claimed only for the datasets actually tested.

---

## 1. Backup branch and verified recovery commit

| Item | Value |
|---|---|
| Validated application commit | `5a20fb073aa828a5dfd411309a08f83b2e039baf` (staging, `render-staging`) |
| **Backup branch (GitHub)** | **`backup/sprint5-staging-validated-20261115`** |
| **Backup commit** | **`24f81fe3d3a5d1b87e5bbe08930315eae4f7f693`**: `5a20fb0` plus the hosted-validation report, its progress record, and `validation/hosted_smoke_test.py` |
| Verification | `git ls-remote`, `git fetch`, and the GitHub API branch endpoint all give `24f81fe`; 231 files, identical to local |
| Local tag | `checkpoint/sprint5-staging-validated-20261115` (annotated; not pushed) |
| Untouched | `render-staging` = `5a20fb0` (staging), `sprint4-integration` = `6e122ff` (production), `main` = `eb72013`. No Render service was triggered |

Restore steps are in the progress record, Milestone 0. They are manual; nothing restores automatically.

## 2. Working branch and current commit

`fix/final-polish-20261115`, created from `24f81fe`. Commits:
1. `5d04408`: checkpoint and baseline.
2. `457c3c6`: Milestone 1.
3. `ef67eb7`: Milestone 2.
4. `b0286d4`: Milestone 3.
5. `2f6f787`: Milestone 4.
6. The final commit: the cache regression fix and these reports (§5, §14).

## 3. Issues repaired

| Audit item | Repair |
|---|---|
| **P2-1** "UTR" band across ORF3a–N | The gene track comes from the uploaded GenBank/GFF3 CDS coordinates (`annotation_file_bands`; overlapping CDS of one gene merged). Without a file, a source label is drawn only if its records form one uninterrupted run, so a recurring label such as Fontenele's intergenic "UTR" stays unannotated. The caption names the source |
| **P2-2** influenza gene "NA" blank in XLSX | `StoredTables.get` reads with `keep_default_na=False, na_values=[""]`, so only empty cells are missing. Sprint 2 readers (`read_pair_table`, cluster and region tables) also keep "NA" segment and gene labels |
| **P2-3** silent standard genetic code | Each CDS records `Translation_Table_Source`: stated, verified (table 1 reproduces the GenBank `/translation`), or confirmed by the user (new Advanced setting / `--genetic-code`). Otherwise the CDS is unsupported with the reason, and no amino-acid consequence is predicted. A stated table is never overridden |
| No feedback during long runs | Banner with the step's message, a spinner, and the real elapsed time (no invented %); repeated submissions ignored; reset on Back |
| No feedback for on-demand XLSX/PDF | "Preparing Excel workbook / PDF report…", then "ready" or the server's error. On the server, one build per file: concurrent clicks wait and reuse it |
| Upload controls 1,780 px down; file help late | "Which files do I need?" leads Step 1; file fields before the optional extras; "How this works" collapsed; wider upload card |
| Review problems at the bottom | A status box at the top ("Ready to run" with a link to Run, or the required fixes); numbered parts |
| Results about 21 screens, weak navigation | Jump bar with Downloads and the Sprint 5 parts (Uncertainty, Mutation annotations, Experimental grouping); empty sections marked "(n/a)"; sticky one-line bar on phones |
| **P2-5** mapping-profile values (`header_row: -5`) | Strict validation: keys, roles, field names, text, structure, and option values under the review form's own rules. Refused as a whole with a specific message |
| **P2-7** cache eviction between users | Memory budget (48 MB compressed by default, `RESULT_CACHE_MB`), at most 12 sets (was 4), least-recently-used eviction where creating and downloading both count as use; 10-minute lifetime and random IDs unchanged. Expired links explain themselves |

## 4. Source files changed

Compared with `24f81fe`: 21 files, about 1,207 lines added and 120 removed.

**Application code:**
- `app.py`: settings, report builds, result cache, expired-link page, `TEXT_LABEL_COLUMNS` readers;
- `compare_genetic_distances.py`: "NA" labels;
- `genomics/genome_annotation.py`: genetic-code provenance;
- `genomics/input_inspection.py`: review warning for an unstated code;
- `genomics/integrated_analysis.py`: gene track;
- `genomics/mapping_profiles.py`: validation;
- `genomics/mutation_annotation.py`: reuses the loaded annotation;
- `genomics/report_exports.py`: "NA";
- `run_integrated_analysis.py`: `--genetic-code`.

**Templates:**
- `templates/index.html`: feedback script, CSS, review status, navigation, genetic-code setting;
- `_upload.html`: guide and order;
- `_how_it_works.html`: collapsed;
- `_review_import.html`: `data-no-busy`;
- `_toolbar.html`: keeps the setting.

**Tests:**
- new: `test_final_polish_annotation.py` (12), `test_final_polish_feedback.py` (4),
  `test_final_polish_navigation.py` (3), `test_final_polish_robustness.py` (9);
- changed on purpose (documented in each test): `test_genome_annotation.py` (9 calls),
  `test_snv_consequences.py` (2), and `test_sprint5_phase4_web.py` (4 uploads). These now confirm the
  standard genetic code for GFF3 files, which state none; before, it was assumed silently.

**Docs:** this report and the progress record. `INTEGRATION_IMPLEMENTATION_LOG.md` was not changed.

## 5. Scientific and numerical regression results

**Test suites:**
- **Baseline at `24f81fe`:** `Ran 685 tests in 616.7 s — OK (skipped=4)`, with real MAFFT,
  `FONTENELE_MMC2`, `CRITS_CHRISTOPH_DIR`, `SCHUMANN_DIR` and `ANNOTATION_MPOX_GENBANK`.
- **First full run after Milestones 1–4:** 712 tests with **1 error**,
  `test_upload_enrichment.test_genetic_outputs_are_unchanged_by_enrichment`.
  - The first cache design ranked a just-created, not-yet-downloaded result below recently downloaded
    ones, so a full cache evicted the newest result.
  - Fixed by plain least-recently-used eviction, with a new regression test
    (`test_a_just_finished_analysis_is_never_the_one_evicted`).
- **Final full run:** **`Ran 713 tests in 634.5 s — OK (skipped=4)`**. That is the 685 baseline tests
  plus 28 new, with the same 4 optional skips: 3 SciPy cross-checks and the NCBI Virus files test.

**Published datasets** (`validation/check_published_datasets.py`): **33 PASS, 0 FAILED**.
- Fontenele: 7,973 records and 52 samples; NC_045512.2 verified.
- Crits-Christoph: 294 records; S4 equals the independent derivation; WIV04 kept separate.
- Schumann: 5,145 records; 38 of 38 per-sample totals; **1,830 of 1,830 VEP amino-acid labels**.
- Original files unchanged (md5).

**Biological checks:**
- **Genetic code**, synthetic CDS ATG TGG AAA TAA (hand-derived):
  - TGG>TGA is **W2\*** under code 1 and **W2W** under code 4;
  - with no stated code, no protein, and no confirmation, no amino-acid change is given (the old
    false "W2\* Stop gained" is gone);
  - a stated `transl_table=4` is used even if the user chooses 1;
  - a GenBank `/translation` verifies table 1, and a disagreeing protein makes the CDS unsupported.
- **NC_045512.2 GenBank:** all 12 CDS supported and `TABLE_VERIFIED`; D614G unchanged.
- **Gene track:** bands equal the NCBI CDS (S 21563–25384, ORF3a 25393–26220, N 28274–29533,
  ORF1ab 266–21555); 25385–25392 unannotated; no "UTR" band, with or without the file.
- **Real H5N1** (andersen-lab iVar tables, confirmed code 1): 127 annotated mutations, 10 on gene NA
  (for example NA:P3A). **"NA" appears in 10 of 10** XLSX cells (it was blank before); no blank
  genes; the PDF builds.

Verified methods are unchanged: MAFFT and distances, temporal/geographic screening, SNV observations,
VAF and depth, Wilson and beta-binomial, replicate eligibility, Yue–Clayton and PCoA, consequences
(where the code was verifiable), exports, missing-data handling, and grouping off by default. This
is shown by the suite and the output comparison in §9.

## 6. UI changes

See §3. No page was redesigned and no figure was duplicated; every result stays on the page and in
the downloads.

## 7. Browser test findings

**Scripted only.** These checks used the built-in browser against the local server at 1024×768,
1280×800, and 375×812. They are not a human usability test.

**Feedback banner and downloads:**
- The banner appears at once with the message and a running clock. The 2nd and 3rd clicks were
  blocked, and Back reset the page.
- The XLSX and PDF show "Preparing…", ignore repeated clicks, then show "ready". An expired link shows
  a clear failure message.

**Homepage** (1024×768):

| Element | Before | After |
|---|---|---|
| File guide | 2,010 px (after the button) | **374 px** |
| File picker | 1,778 px | **956 px** |
| Review button | 1,939 px | **1,332 px** |
| Page height | 2,737 px | 2,130 px |

**Phone:** no sideways page scrolling; the results bar is a single 40 px line that stays at the top.

**Review page:** the status box is on the first screen (504 px), and the sticky Run button is visible.

**Results page:** Downloads at 696 px and Major findings at 1,165 px. All jump links resolve.

**Needs a human test:**
- At 1024×768 the file picker is still one scroll down, because the mode chooser comes first.
- The review page (about 6.7k px) and results page (about 15.8k px) are long, but navigable.
- Whether first-time users understand the genetic-code setting and the "(n/a)" sections.

## 8. Memory and performance results

Local imitation of Render Free (512 MB, 10% of one core), `mmc2.xlsx` through the website:

| Measure | Before (`5a20fb0`) | After |
|---|---|---|
| Run time | 123.5 s | 127.1 s (measured while the full suite ran in parallel) |
| Peak memory | 302 MB | **300 MB** |
| Committed memory | 304 MB | 302 MB |
| XLSX | 39.8 s | 41.0 s |
| PDF | 5.4 s | 5.9 s |
| HTML | 3.02 MB | 3.03 MB |

**Local full speed:**
- `mmc2` + FASTA + GenBank 13.1 s (12.5 s), integrated example 14.3 s (14.6 s), Sprint 2 sets unchanged.
- The cache can now hold up to 48 MB compressed (before: 4 sets, about 3–5 MB each).
- No memory or timeout regression was found.
- **Hosted performance has not been re-measured** for this branch.

## 9. Output-equivalence comparisons

The browser-like smoke test (`validation/hosted_smoke_test.py`) was run against the new code and
compared with the saved baseline of the validated version (same inputs and settings, GeoNames unset):

| Scenario | CSVs identical in content |
|---|---|
| Sprint 2 synthetic | 11 / 11 |
| Sprint 2 NCBI H5N1 ten-record set | 11 / 11 |
| Integrated example | 44 / 44 |
| Fontenele `mmc2.xlsx` + FASTA + GenBank | 33 / 34 |
| Formula-safety table | 30 / 30 |
| **Total** | **129 / 130** |

- **The only difference** is the intended new column `Translation_Table_Source` in
  `sprint5_annotation_sources.csv`.
- **Also identical:** main findings, readiness badges, the Sprint 4 H5N1 result (the 2022 article is
  rejected, nothing is changed), the formula check (0 formula cells), unsupported-input handling, and
  the busy-server message.
- **Figures:** the genome-position figure changed on purpose (correct gene track, new caption).

## 10. Remaining P2/P3 issues

- **Schumann Berlin35 memory: not implemented.**
  - Categorical storage would save 533 → 50 MB with identical CSV bytes.
  - But at least 9 consumers group or join on those columns, and pandas 2.x categorical groupby
    (`observed=False`) creates every category combination.
  - It needs a consumer-by-consumer audit. `MAX_OBSERVATION_CELLS` keeps refusing it clearly on the
    free plan.
- **Start codon:** a change to another start codon (ATG>CTG, or GTG>ATG under table 11) is still
  reported as "Start lost".
- **Remaining audit P3s:**
  - percent values all ≤ 1 cannot be detected;
  - "Coverage" breadth can be auto-mapped as depth;
  - a frequency that disagrees with counts is not raised at review;
  - base-count files are not used as coverage unless also uploaded as positions;
  - two VAF denominators;
  - the slope-rule inconsistency between two trend modules;
  - float32 Yue–Clayton;
  - the PCoA explained-share label;
  - raw exception text in a few messages;
  - CSRF tokens and security headers;
  - `requests` minimum version;
  - `defusedxml` not installed;
  - GFF3 URL decoding;
  - GTF message.
- **Pathogen support:** segmented and multi-chromosome genomes still need one reference sequence per
  run.

## 11. Known scientific limitations

- Wilson intervals model read sampling only. Beta-binomial calibration needs real technical
  replicates, and none of the three published datasets has them.
- Trends are descriptive.
- An annotation predicts an amino-acid change, not a phenotype. Consequences are given only where the
  genetic code is stated, verified, or confirmed by the user.
- Experimental mutation-frequency grouping is **not validated** and stays off by default.
- Real observations were tested for SARS-CoV-2 (three studies) and one H5N1 influenza dataset only.
  No broader pathogen compatibility is claimed.
- Live Claude NLP has never been run or validated.

## 12. Render staging test instructions (when approved)

1. Push only this branch, without force, after review:
   `git push origin fix/final-polish-20261115`.
2. In the **staging** service only, change the branch from `render-staging` to
   `fix/final-polish-20261115`, or push the branch's commit to a new staging-only branch, and deploy.
   Production and `render-staging` stay as they are.
3. Read the build log for the commit SHA, `mafft`, and gunicorn `gthread`.
4. Run the browser-like test against staging:
   `python validation/hosted_smoke_test.py --base https://pathogen-clustering-2.onrender.com --output <dir> --data <dir with mmc2.xlsx, sequences.fasta, sequences.csv>`.
   Compare the CSVs with the hosted baseline. Only `Translation_Table_Source` should differ, plus the
   GeoNames-dependent Sprint 2 geography already seen on staging.
5. In the browser:
   - check the progress banner during the about 3-minute Fontenele run;
   - check that "Preparing Excel workbook" appears (about 70 s on staging);
   - check the homepage guide and the review status box;
   - check the jump bar on a phone.
6. Read Render Metrics → Memory during the Fontenele run, the one number the earlier hosted
   validation could not see.

## 13. Steps required before production promotion

1. Done: **`Ran 713 tests in 634.5 s — OK (skipped=4)`**, the 685 baseline tests plus 28 new, in the
   same environment as the baseline.
2. Staging revalidation of this branch (§12), with explicit approval.
3. A short human usability test (hosted report §12 checklist), focusing on the new upload guide,
   review status, progress banner, and genetic-code setting.
4. Explicit approval to promote. Then make a normal (non-force) push of the validated commit to
   `sprint4-integration`, and repeat the smoke test against production.
5. Keep `backup/sprint5-staging-validated-20261115` unchanged as the recovery point.

## 14. Completion assessment

| Required criterion | Status |
|---|---|
| Remote backup created and verified | **Met** (`24f81fe`, GitHub API) |
| Staging and production unchanged | **Met** (`render-staging` `5a20fb0`, `sprint4-integration` `6e122ff`; no deploy) |
| Separate working branch | **Met** |
| Three scientific-correctness issues | **Met** (UTR band, "NA", genetic code), with independent tests |
| Long-run and download feedback | **Met** (browser-verified) |
| Upload and result navigation | **Met** (browser-measured; human test pending) |
| Mapping-profile validation | **Met** |
| Safer download cache | **Met** (bounded by memory and count, LRU, isolated, tested) |
| All relevant regression tests pass | **Met** (713 OK, 4 optional skips; 33/33 published checks) |
| Existing scientific outputs correct | **Met** (129/130 CSVs identical; 33/33 published checks) |
| Performance regressions evaluated | **Met** (local imitation; hosted not re-measured) |
| Report documents changes and limitations | **Met** |

**FINAL POLISH COMPLETE — READY FOR STAGING REVALIDATION**

Every required criterion is met with local evidence. This is not a hosted result: this branch has not
been deployed. Staging revalidation (§12) and a human usability test are the next steps before any
production promotion.
