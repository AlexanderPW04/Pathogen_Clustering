# Final Polish: Staging Deployment and Hosted Validation Report

**Date:** 2026-11-16
**Staging URL:** https://pathogen-clustering-2.onrender.com

**How each result is labelled:**
- **Hosted:** observed on the Render staging service.
- **Local:** observed on this computer.
- **Requested:** only the Render dashboard holds the evidence; it is listed in §15.

Expected results from earlier reports are quoted only as comparison values.

---

## 1. Recovery branch and backup commit (verified before any change)

- **Recovery branch:** `backup/sprint5-staging-validated-20261115` on GitHub at
  **`24f81fe3d3a5d1b87e5bbe08930315eae4f7f693`** (`git ls-remote`, after `git fetch`). It is unchanged at
  the end of this task.
- **History:** `5a20fb0` (the previous staging version) and `24f81fe` are ancestors of everything
  deployed below, so the validated version stays recoverable.
- **Secrets:** no tracked `.env`, `.pem`, `.key` or credential files, and no key patterns in the
  branch. `.env` is git-ignored; the local `.claude/launch.json` is excluded through
  `.git/info/exclude`.

## 2. Previous staging commit

**`5a20fb073aa828a5dfd411309a08f83b2e039baf`** on `render-staging`, validated on 2026-11-15
(`docs/RENDER_HOSTED_VALIDATION_REPORT.md`).

## 3. New staging commits

All three updates were normal fast-forward pushes of `render-staging` only. Git confirmed each
previous tip was an ancestor; there was no force push and no rewrite.

| Push (UTC) | `render-staging` | Content | Live (new code served) |
|---|---|---|---|
| 08:46:50 | `5a20fb0 → 8287576` | The final-polish branch as tested locally | 08:47:43 (homepage file guide present) |
| 09:13:53 | `8287576 → 84d4f64` | Staging repair 1 (§17), plus local test tools (`validation/`, not in the image) | 09:14:43 |
| 09:19:39 | `84d4f64 → 0119813` | Staging repair 2 (§17) | 09:20:14 |
| 10:03:58 | `0119813 → 8690561` | Staging repair 3 (§17), caption text only | confirmed by a hosted run (caption names the uploaded file) |

**Now live:** **`86905613f49be24243da318be72d0062e1ddba83`**, branch `render-staging`.

**Unchanged on GitHub (verified at the end):**
- `sprint4-integration` = `6e122ff` (production);
- `main` = `eb72013`;
- the backup branch = `24f81fe`;
- `fix/final-polish-20261115` remote = `8287576`. Its local copy also has the three later commits.

Render deployment identifiers: **Requested** (§15).

## 4. Render service configuration

- **Staging builds from `render-staging`:** shown by the user-supplied build log of 2026-11-16 06:03
  ("Checking out commit … in branch render-staging"). Auto-deploy is on: each push went live within
  30–55 s, without any manual deploy.
- **Production builds from `sprint4-integration`:** as stated by the user. **Not independently
  verified:** the Render dashboard stayed on its SSO login page in the built-in browser.
- **Paid Claude is off:** the hosted status page reports "Claude API configuration: Disabled … No paid
  requests are possible".
- **External services:** GeoNames is configured and works; PubTator works.
- **Limits** (from the Dockerfile, observed in behaviour): a 1.6 M-cell CSV is refused with "the limit
  is 1,500,000", and an upload over 10 MB gets a 413.

## 5. Docker build and startup findings

- **Hosted:** each push was served within about a minute. `/health` returned 200 in 0.1 s after each
  deploy, and the homepage served the expected new markup:
  - the "Which files do I need?" guide;
  - the `genetic_code` setting;
  - the progress banner;
  - the scroll-offset script;
  - the old "up to two minutes" text gone.
- **Not inspected by Claude:** the build logs, startup lines (MAFFT install, gunicorn `gthread`) and
  runtime logs of the three new deploys, because Render could not be signed in. **Requested** (§15).
- **Indirect evidence for MAFFT:** the hosted status page reports MAFFT v7.505, and Sprint 2
  alignments ran (§7).
- **Indirect evidence against restarts:** no request failed, and no error page appeared across about
  45 minutes of hosted tests.

## 6. Full local test results

| Run | Commit | Result |
|---|---|---|
| Full suite before deployment (this task) | `8287576` | **`Ran 713 tests in 652.4 s — OK (skipped=4)`** |
| UI modules after each staging repair | `84d4f64`, `0119813` | 28 tests OK (1 optional skip) each |
| Full suite before the last staging repair | `8690561` | **`Ran 714 tests in 1307.4 s — OK (skipped=4)`** (one new regression test) |

The 4 skips are the 3 SciPy cross-checks and the NCBI Virus files test. The run used real MAFFT plus
`FONTENELE_MMC2`, `CRITS_CHRISTOPH_DIR`, `SCHUMANN_DIR` and `ANNOTATION_MPOX_GENBANK`.

- Repairs 1–2 change only `templates/index.html` and one test file.
- Repair 3 changes caption text in `genomics/integrated_analysis.py`. That is why the full suite ran
  again; the slower 1,307 s reflects other load on the computer.

## 7. Sprint 2 hosted results

Summarised in §14 (Sprint 2 paragraph). Same inputs and settings as the 2026-11-15 hosted baseline.

## 8. Sprint 3 + 5 hosted results

All results come from the hosted run on `0119813` (09:22–09:36 UTC), with the same inputs and settings as
the 2026-11-15 baseline. Repair 3 (`8690561`) changes only a caption (§17).

**A. Synthetic integrated example** (bundled files; uploaded through the page, and also run from a real
browser click three times):
- Readiness: 15 available, 1 partially available, 2 unavailable, 1 experimental.
- 5 individual records, 4 placed on the reference; 9 pooled samples, 8 kept after QC.
- Cross-source: 2 variants in both sources, 1 individual only, 1 pooled only, 1 unresolved (pooled),
  1 unresolved (individual), 2 incomparable.
- Population structure: the paper-method Yue–Clayton PCoA of A/C/G/T frequencies.
- Sprint 5: 22 of 22 pooled observations have intervals.
- **All 44 CSVs are identical to the previous staging version.**
- The grouping panel says groups are not shown by default.

**B. Fontenele `mmc2.xlsx` + NC_045512.2 FASTA + NCBI GenBank:**

| Check | Hosted result |
|---|---|
| Samples / observations | **52 / 7,973**; QC paper preset keeps 7,961 |
| Reference | NC_045512.2 **Verified compatible** |
| Provenance | Frequencies "Reported by the source" (all rows); alternative counts **"Reconstructed"** (all rows) |
| Uncertainty | 7,970 of 7,970 observations have a 95% read-sampling interval (from reconstructed counts) |
| Overdispersion | **Calibration skipped: "No documented technical replicates … cannot be estimated"**. Only rho = 0.002/0.02/0.05 "scenarios, not estimates". No fabricated estimate |
| Missing coverage | "No per-position depth: unreported variants stay 'coverage unknown' (never 0%)" |
| Temporal | 52 of 52 samples dated; 474 mutation-site series |
| PCoA | Exploratory Jaccard PCoA (no A/C/G/T counts), labelled as not the paper's method |
| Annotation | 5,677 of 5,677 distinct SNVs annotated; 12 of 12 CDS supported; all 12 **"Not stated; the standard code (table 1) reproduces the record's own protein sequence"** |
| Grouping | Mode `profiles` (default): "groups not shown by default because the method did not meet its validation criteria" |
| Genome figure | Gene track from the GenBank CDS (ORF1ab, S, N labelled; the smaller genes narrow), **no "UTR" band**. The caption names the source file (after repair 3) |
| Without optional files (`mmc2.xlsx` alone) | Runs (164.3 s). Readiness: 4 available, 5 partially available, 5 unavailable, 1 experimental. The caption explains that source labels are drawn only as uninterrupted runs and that "UTR" is left unannotated |

**Unavailable analyses** are explained on the page and in `analysis_status.csv`: "Needs individual
sequences…", "Genome coverage known for 0 of 52…", "Needs mean depth for at least three…".

## 9. Biological annotation verification (hosted)

Script: `validation/hosted_polish_checks.py`. It ran against `8287576`; the later commits change no
analysis code.

**Genetic code.** Synthetic CDS ATG TGG AAA TAA, SNV TGG>TGA, GFF3 with no `transl_table`:

| Setting | Result on staging |
|---|---|
| Not confirmed (default) | CDS **Unsupported**, source "Not stated in the annotation"; consequence "Unsupported or ambiguous", **no amino-acid change**; the review page warns |
| Code 4 confirmed | **codeA:W2W, Synonymous**; source "Confirmed by the user … (table 4)" |
| Code 1 confirmed | **codeA:W2\*, Stop gained**; source "Confirmed … (table 1)" |

**Gene "NA"** (real andersen-lab H5N1 iVar tables for the NA segment, NCBI annotation, code 1
confirmed):
- 10 mutations on gene NA (for example NA:P3A, NA:Q5Q, NA:G11G).
- **"NA" survives in 10 of 10 CSV rows, 10 of 10 XLSX cells, and the website**, with no blank genes.
- 7 of the 8 NA tables were used. One table holds only insertions/deletions and is refused at review;
  see §16.

**Genome-position figure (Fontenele):** see §14.

**Existing predictions kept:**
- Fontenele annotation CSVs are identical to the previous staging run.
- NC_045512.2 GenBank coding regions are verified through `/translation`, so no valid prediction was
  removed.

**Known limitation:** a change from one start codon to another (for example ATG>CTG, or GTG>ATG under
table 11) is still reported as "Start lost". Alternative initiation codons are not claimed to be
fully supported.

## 10. Long-running operation feedback (real browser on staging)

**Run banner** (built-in browser, "Run Integrated Example", the real hosted run):
- It appeared immediately with "Running the example analysis. This can take a minute or more." and
  "(elapsed 0:00)".
- No percentage was shown; the help text says "No progress percentage is shown because the server
  does not report one".
- Of three rapid clicks (Run, Run, another example), the first was **sent** and the next two were
  **blocked**.
- Completion went straight to the results page.

**Report downloads:**
- **XLSX:** "Preparing Excel workbook… the download starts by itself", the link marked busy, and three
  clicks produced **one** request. About 4 s later it read "Excel workbook ready: the download has
  started."
- **PDF, failure then retry:** an expired link gave "The PDF report could not be prepared: These
  results are no longer on the server … kept in memory for 10 minutes … You can try again". The retry
  on the real link gave "PDF report ready" (about 9–15 s).
- **Duplicate builds on the server:** three concurrent requests for the same workbook all returned 200
  with identical bytes (scripted). They were served from the cached copy, so concurrent *building* was
  covered by the local test (`tests/test_final_polish_feedback.py`, one build for three requests), not
  on staging.
- **Server health:** `/health` was 200 after each operation and during a run.

**Limitations of this check:**
- The browser pane was hidden, which slows page timers. The clock was read at 0:00 and not watched
  over a full run on staging; it advanced correctly (0:01, then 0:02) in the earlier local browser
  check.
- No screenshots were taken on staging; layout was measured from element positions.

## 11. Browser navigation findings (staging)

**Homepage:**

| Viewport | File guide | Mode choice | File picker | Review button | "How this works" |
|---|---|---|---|---|---|
| 1024×768 | 374 px | 595 px | 956 px | 1,332 px | collapsed |
| 375×812 | 595 px | — | 1,490 px | 1,954 px | collapsed; no sideways scrolling |

**Review page** (1024×768, bundled example files uploaded through the page):
- The status box is at 504 px: "Ready to run. 5 files recognized. 1 note about optional information
  is listed below…".
- The sticky Run button is visible at 707 px.
- With all roles set to "Do not use", the box reads "Required information to fix before running — No
  file has been assigned a role yet". The heading appears once and Run is disabled.

**Results** (integrated example):
- 14 jump links, including Downloads, Uncertainty, Mutation annotations and Experimental grouping;
  all targets exist.
- Downloads is at 920 px and Major findings at 1,414 px (1280×800). No figure is wider than the phone
  screen, and the page never scrolls sideways.
- **Two defects found here and repaired** (§17): headings landed under the sticky bar after a jump.
  After the repairs, **0 of 14 jumps land under the bar** at 1024×768 (bar 129 px) and at 375×812
  (bar 40 px). The 1280×800 case was verified locally (bar 101 px).
- The grouping panel states that groups are "not shown by default because the grouping method did not
  meet its pre-set validation criteria".

**Remaining for a human test:**
- At 1024×768 the file picker is one scroll down.
- Results pages are long (about 16k px) but navigable.

## 12. Sprint 4 findings (hosted)

Summarised in §14 (Sprint 4 paragraph).

## 13. Downloads and cache tests (hosted)

| Check | Result |
|---|---|
| Invalid profile, `header_row: -5` | "The mapping profile's entry 1 option 'header_row' is -5; it must be a whole number from 1 to 50." The review stays open |
| Invalid JSON profile | "The mapping profile is not valid JSON." |
| Profile role `../../app` | "unknown role '../../app'" |
| Old profile (saved by `5a20fb0` on staging) on `mmc2.xlsx` | "Profile entry 1 (snv) applied to mmc2." |
| Formula injection (`=HYPERLINK` sample ID) | 0 XLSX formula cells; 0 CSVs with an unescaped formula |
| Three concurrent identical XLSX requests | All 200, identical size, valid workbooks |
| PDF twice | 200, 200, identical bytes, valid `%PDF-` |
| Another token, unknown file in valid results, `..%2f` traversal (×3) | All 404 |
| Binary garbage `.csv` | Accepted as a table with 0 records; Run disabled; no error |
| Broken `.xlsx` | "Could not read this file: File is not a zip file" |
| FASTQ / 1.6 M-cell CSV / 11 MB upload | Refused with clear messages (200 / 200 in 0.3 s / 413) |
| Two runs at once | The second gets "Another analysis is running … Your files are still here"; the first completes |
| `/health` after all failures | 200 |

- **Cache bounding** (memory budget, least-recently-used eviction, isolation, 30 concurrent stores):
  verified by the local tests `tests/test_final_polish_robustness.py`. It was not stress-tested on
  Render Free, as instructed.
- **Cold start:** a cold start was observed once: the first request after idling took **23.3 s**
  (08:46 UTC, before the push).
- **Not tested:** a deliberate restart. Results are in memory by design, so they are lost on restart,
  and the 404 page says so.

## 14. Local-versus-hosted scientific output comparisons

Comparisons of the CSV downloads by content (numeric tolerance 1e-9 relative):

| Comparison | Identical CSVs | Differences |
|---|---|---|
| Previous staging `5a20fb0` against final staging `0119813` | **129 / 130** | Only `sprint5_annotation_sources.csv`: new column `Translation_Table_Source` (intended) |
| First new deploy `8287576` against final `0119813` | **130 / 130** | none |
| New code locally (GeoNames unset) against new staging | 119 / 130 | Python version in `provenance.csv` (3.11.9 locally, 3.11.18 on staging); GeoNames configured on staging gives one place approximate coordinates in the Sprint 2 synthetic set (distances filled in). Strict pairs and clusters unchanged |

**Main findings and readiness badges are identical** to the previous staging version in all six
scenarios.

**Sprint 2:**
- Synthetic: **4 strict related pairs, 2 clusters** (0 pending geography, because GeoNames is
  configured).
- H5N1 ten-record set: **0 strict pairs; 2 genetic-and-date candidates blocked by state-level
  geography**. Both compare the same segment (Segment 1 with Segment 1, Segment 2 with Segment 2), with
  0 differences, 0 days apart and no kilometre distance.
- MAFFT alignment, distances, dates, geography precision and downloads were all exercised.

**Sprint 4** (hosted, free search only):
- "10 accession records searched. One publication found. No verified collection locations could be
  assigned … Ten accessions represent two distinct isolates. One associated publication was retrieved.
  The publication is from 2022; the samples were collected in September 2026 … cannot be applied to
  your samples. Claude was not used…".
- The 2022 article (PMID 35302933) is **not** used for the 2026 samples, and no metadata is changed.
- Continuing with the original metadata works (verified on 2026-11-15; the Sprint 4 code is unchanged
  since).
- No paid call was made. Live NLP is not validated.

## 15. Actual memory and runtime observations

Hosted run times (scripted client; "previous" = the 2026-11-15 run on `5a20fb0`):

| Operation | Previous staging | First new deploy (`8287576`) | Final staging (`0119813`) |
|---|---|---|---|
| Homepage | 0.18 s | 0.09 s | 0.09 s |
| Cold start (first request after idle) | — | 23.3 s (before the push, old code) | — |
| Review (Fontenele upload) | 13.1 s | 14.2 s | 13.2 s |
| Integrated example run | 71.2 s | 74.9 s | **70.8 s** |
| Fontenele without optional files | 163.9 s | 168.3 s | **164.3 s** |
| Fontenele + FASTA + GenBank | 178.8 s | 183.2 s | **177.6 s** |
| Fontenele XLSX (on demand) | 68.6 s | 69.1 s | **68.2 s** |
| Fontenele PDF | 9.4 s | 9.8 s | 8.4 s |
| Busy-server reply to a second run | 11.1 s | 11.5 s | 19.0 s |
| Sprint 4 search | 7.0 s | 9.7 s | 8.9 s |

- **No significant regression:** final-staging times are within about 1% of the previous version.
- **Errors and workers:** no request failed, and every page returned 200 (or the intended 404/413).
  The worker never stopped answering: `/health` was 200 during and after the long runs.

**Memory:** the Render Metrics memory graph has **not** been seen; Claude has no dashboard access, so
no hosted memory value is claimed. Locally, under the 512 MB / 10% CPU imitation, the same commit's
Fontenele run peaked at 300 MB.

**Requested from the Render dashboard** (staging service):
1. **Events:** the four deploys (commits `8287576`, `84d4f64`, `0119813`, `8690561`) with their
   identifiers.
2. **Logs:** the build lines for MAFFT and pip, `Listening at`, and `Booting worker` (one per deploy),
   plus any traceback or `WORKER TIMEOUT` between 08:47 and 09:45 UTC.
3. **Metrics → Memory:** the peak during the Fontenele runs.

## 16. Remaining limitations

- The hosted memory, build and runtime logs for the new deploys are not seen (§15). Production's
  branch link is user-stated, not independently verified.
- **Unnamed refusal:** when one of several uploaded iVar tables contains only insertions/deletions, the
  review says "The processed SNV table has no records" without naming the file. It is safe, because
  nothing runs, but it is unhelpful.
- **Start codon:** see §9.
- **Schumann Berlin35:** still refused on the free plan by `MAX_OBSERVATION_CELLS`. The memory
  optimization is deferred. It was not run on staging, as instructed.
- **Live Claude NLP** is not validated. Pathogen support is shown only for SARS-CoV-2 (three published
  studies) and one H5N1 dataset.
- **Not observed on staging:** a full-length progress-clock observation (hidden browser pane), a
  deliberate restart, and a human usability test.
- **Remaining audit P3 items:** listed in `docs/FINAL_POLISH_IMPLEMENTATION_REPORT.md` §10.

## 17. Staging-only repairs

Both repairs were found on staging, fixed with the smallest change, tested locally, deployed by
fast-forward, and re-verified on staging. Neither touches analysis code.

1. **`84d4f64`: Sprint 5 jump targets.** "Uncertainty", "Mutation annotations" and "Experimental
   grouping" pointed at panels without the scroll offset that sections have. Their headings landed
   under the sticky bar on phones (top 0 px against a 40 px bar). One CSS selector fixed it, and a
   regression assertion was added.
2. **`0119813`: offset from the bar's real height.** At 1024 px the jump bar wraps to three lines
   (129 px), taller than the fixed 6.5rem (104 px) offset, so **every** heading landed about 25 px
   under it. This was partly caused by the links added in Milestone 3. A small script now stores the
   bar's height in a CSS variable. Verified at 1024, 1280 (local) and 375 px: 0 of 14 jumps hidden.
3. **`8690561`: caption provenance.**
   - The new gene-track caption named the server's copy ("annotation.txt") instead of the uploaded
     file, because the track reads the annotation before the upload name is attached.
   - Two lines now set the uploaded name. A regression test was added, and the full suite (714) passes.
   - Confirmed on staging: "…from the uploaded annotation (NC_045512.2.gb)".
   - Caption text only; no table or calculation changes.

## 18. Production-readiness assessment

**Ready for production promotion after the §19 steps, not before.**

**Met:**
- The scientific outputs are unchanged on staging (129/130, with the one intended column).
- All three annotation corrections work on staging.
- The safeguards hold: no fabricated overdispersion, grouping off, coverage unknown never 0%,
  segments kept apart, the 2022 article rejected, formula-safe downloads, and strict profiles.
- The feedback and navigation improvements work in a real browser, after the three staging repairs.
- No paid calls were made, and production and `main` are untouched.

**Missing before production:**
1. The Render evidence that only the dashboard holds (§15): deploy identifiers, build and runtime
   logs, and the memory peak.
2. Independent confirmation that production builds from `sprint4-integration`.
3. A human usability test.

## 19. Steps for an eventual production promotion

1. Collect the §15 evidence: deploy events, build and runtime logs, and the memory peak.
2. Run a short human usability test (hosted report §12 checklist).
3. With explicit approval only, push the validated commit `8690561` to `sprint4-integration` with a
   **normal fast-forward push**, if `6e122ff` is an ancestor; otherwise stop and review. Do not merge
   into `main` unless separately decided.
4. Watch the production deploy, then run
   `python validation/hosted_smoke_test.py --base <production URL> …` and
   `python validation/hosted_polish_checks.py --base <production URL> …`. Compare with this report.
5. Keep `backup/sprint5-staging-validated-20261115` (`24f81fe`) unchanged.
6. The deployed commit `8690561` exists only on `render-staging` and on the local
   `fix/final-polish-20261115`; the remote final-polish branch is still at `8287576`. Before promotion,
   push the final-polish branch normally (fast-forward) so the reviewed branch matches what was
   validated.
   - **To roll back staging:** fast-forwarding is not possible backwards. Point the staging service's
     branch at the backup branch in the dashboard; never force-push.
   - **To roll back production:** production still runs `6e122ff` until it is promoted.

---

**FINAL POLISH HOSTED VALIDATION PASSED WITH LIMITATIONS**

**Passed (hosted, observed):**
- The verified recovery checkpoint is preserved (`24f81fe`).
- Only `render-staging` was updated (normal fast-forwards), and the new code is confirmed live.
- Health and homepage work.
- Sprint 2 and Sprint 3 + 5 match the baseline, apart from the intended column.
- The three annotation corrections work.
- The progress banner, repeated-click blocking, and the XLSX/PDF preparation, failure and retry all
  work in a real browser on staging.
- Profile validation, isolation and formula safety hold.
- Grouping stays off by default, and Sprint 4 free discovery works.
- No paid Claude requests were made, and production and `main` are unchanged.

**Limitations:**
- The Render build and runtime logs and the memory graph for the new deploys were not seen (§15).
- Production's branch link is user-stated.
- The progress clock was not watched over a full hosted run, because the browser pane was hidden.
- No deliberate restart and no human usability test were done.
- Three small staging-only UI/provenance repairs were needed (§17).
