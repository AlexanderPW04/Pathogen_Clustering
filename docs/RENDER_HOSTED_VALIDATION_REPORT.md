# Render Staging Hosted Validation Report

**Date:** 2026-11-15
**Staging URL:** https://pathogen-clustering-2.onrender.com
**Git:** commit `5a20fb0`, branch `render-staging`

**Scope.** Every hosted test below ran against the staging service on Render. No paid Claude request
was made; AI enrichment was off throughout. The production service and its branch were not changed.

**How each result is labelled:**
- **Hosted (verified):** observed on the staging service.
- **Local:** observed on this computer.
- **Requested:** evidence that only the Render dashboard holds; it is listed in §14.

The resume record is [RENDER_HOSTED_VALIDATION_PROGRESS.md](RENDER_HOSTED_VALIDATION_PROGRESS.md).

---

## 1. Git branch and deployed commit

**Production deploys from `sprint4-integration`** (stated by the user). The goal's step "push
`sprint4-integration`" would therefore have redeployed production.

**The user chose a new branch instead:**
- Commit `5a20fb0` was pushed to a new GitHub branch, `render-staging`. The user ran the push.
- The staging service builds from `render-staging`.

**GitHub after the push** (`git ls-remote origin`, verified):

| Branch | Commit | Changed? |
|---|---|---|
| `render-staging` | `5a20fb073aa828a5dfd411309a08f83b2e039baf` | new |
| `sprint4-integration` | `6e122ff…` (production source) | unchanged |
| `main` | `eb72013…` | unchanged |

Nothing was merged or force-pushed.

**The deployed commit is confirmed by the Render build log** that the user supplied:

> `==> Checking out commit 5a20fb073aa828a5dfd411309a08f83b2e039baf in branch render-staging`

The site's behaviour agrees (hosted, verified):
- The new upload cell limit refuses a 1.6 M-cell CSV with "the limit is 1,500,000". Both the check
  and the `MAX_WORKSHEET_CELLS=1500000` Dockerfile setting were added in `5a20fb0`.
- Downloads hold no live spreadsheet formulas.
- The homepage is byte-identical to the local `5a20fb0` homepage (54,382 bytes).

## 2. Render staging service configuration (no secrets)

| Setting | Value |
|---|---|
| Type / runtime | Web Service, Docker, repository `Dockerfile` |
| Branch | `render-staging` |
| Instance | Free (as instructed; plan not seen in the dashboard by Claude) |
| Health check | `/health` (`app.py:256`, returns `OK`) |
| `AI_ENRICHMENT_ENABLED` | `false`. The status page reads "Disabled … No paid requests are possible" (hosted, verified) |
| `ANTHROPIC_API_KEY`, `AI_ACCESS_CODE` | not set |
| `GEONAMES_USERNAME` | set by the user. The status page reads "GeoNames: Operational" (hosted) |
| `NCBI_EMAIL` / `NCBI_API_KEY` | optional; the free NCBI/PubTator lookups worked (hosted) |
| From the Dockerfile | `MAX_OBSERVATION_CELLS=700000`, `MAX_WORKSHEET_CELLS=1500000`, 1 gunicorn worker × 4 threads, `--timeout 180`, non-root user |
| Flask session secret | **Not needed**: the app sets no cookies and has no `SECRET_KEY` setting |

## 3. Docker build outcome

**Hosted. The build succeeded.** Render build log, 2026-11-16 06:03 UTC, supplied by the user:
- **Base image:** `python:3.11-slim-bookworm@sha256:a83ba18a…`.
- **MAFFT:** `apt-get install mafft` installed `mafft 7.505-1`, with Perl as a dependency.
- **Python packages:** `pip install -r requirements.txt` succeeded: pandas 2.3.3, numpy 2.4.6,
  openpyxl 3.1.5, Flask 3.1.3, biopython 1.88, requests 2.34.2, gunicorn 23.0.0, matplotlib 3.11.2,
  anthropic 2.8.0.
- **Copy steps:** all 16 build steps finished, including every application `COPY` and the
  `appuser` creation.
- **Image:** exported and pushed.
- **Warnings:** pip's notices about running as root and a newer pip release only.

The status page reports **MAFFT v7.505** aligning a test pair successfully, so the Linux MAFFT
install works. Python is 3.11.18 (from `provenance.csv`).

**Local Docker.** **Not run.** Docker Desktop launched, but its Linux engine never started on this
computer (the `docker-desktop` WSL distribution stayed stopped). There is no local image test.

**Preflight (local, verified):**
- The Dockerfile copies every module `app.py` imports.
- There is no static folder.
- No Windows paths are used at runtime.
- `.env` is git-ignored and excluded from the image.
- No credentials are tracked in Git.

## 4. Runtime and health-check results (hosted)

| Request | Result |
|---|---|
| `GET /health` | 200 `OK`, 0.3–0.5 s |
| `GET /` | 200, 54,382 bytes, 0.2 s |
| `POST /system-status` | MAFFT Operational (v7.505); GeoNames Operational; PubTator3 Operational; Claude "Disabled … No paid requests are possible" |
| Startup exceptions, missing modules or templates | None: all routes exercised returned 200 (or the intended 404/413), and the runtime log has no traceback |

**Runtime log** (supplied by the user, 06:03–06:32 UTC; it spans every hosted test):
- **Startup:** `Setting WEB_CONCURRENCY=1`; `Starting gunicorn 23.0.0`;
  `Listening at: http://0.0.0.0:10000`; `Using worker: gthread`; `Booting worker with pid: 25`;
  "Your service is live" 12 s after boot.
- **Restarts:** there is **only one "Booting worker" line**, so the worker was not restarted during
  the 178.8 s Fontenele runs (around 06:13–06:27).
- **No failures:** no `WORKER TIMEOUT`, `Out of memory`, traceback or `ERROR` line.
- **Sprint 2 progress lines** confirm segment handling:
  - "Pathogen/gene groups (8): Influenza A virus / Segment 1 (2 samples); … Segment 8 (1 sample)";
  - "Pairwise comparisons: 2";
  - the synthetic set reports 4 related pairs, 2 clusters, and "GeoNames coordinates added: 1".
- **No secrets** appear in the log. Only temporary `/tmp/tmp…` paths are printed.

**Not observed:** the Render Metrics memory graph. It is **still requested**; under the local 512 MB
imitation, the same run peaked at 302 MB.

## 5. Sprint 2 hosted results

Script: `validation/hosted_smoke_test.py`. It uploads files, reads the review form exactly as a
browser shows it, presses Run, and downloads every file.

**Synthetic example** (`examples/synthetic_sequences.fasta` + `synthetic_metadata.csv`):
- Run 12.6 s (2.8 s locally).
- 4 strict related pairs, 2 potential clusters, as locally.
- 16 downloads.
- One difference, caused by GeoNames being configured on staging:
  - One named place received approximate GeoNames coordinates.
  - So "pairs with unresolved geography" is 0 on staging and 1 locally, and `Km_Apart` is filled for
    7 comparisons.
  - Strict pairs and clusters are unchanged.

**Ten-record NCBI Virus H5N1 example** (A/Turkey/SD/26G12570-002, segments, state-only location "USA: SD"):
- Run 13.8 s.
- All 10 records read; 0 strict pairs, 0 clusters.
- 2 pairs pass the genetic and date cutoffs but have unresolved geography, because the location is
  state-level only and is never used for kilometre screening.
- **Both pairs compare the same segment:** Segment 1 with Segment 1, and Segment 2 with Segment 2.
  Segments are never combined.
- This matches the earlier validated conclusion (state-only location: strict and candidate results
  kept separate).
- All 11 CSV downloads equal the local baseline in content. The only difference is the Python version
  in `provenance.csv`.
- MAFFT alignment and p-distance, dates, the geography precision table, graphs and downloads all ran.

## 6. Sprint 3 + 5 hosted results

**Synthetic integrated example** (all six bundled files):
- Run 71.2 s.
- Readiness: 15 available, 1 partially available, 2 unavailable, 1 experimental.
- 72 downloads; **43 of 44 CSVs equal the local baseline**, and the 44th differs only in the Python
  version.
- The same example also ran from a real browser click ("Run Integrated Example") on the staging
  homepage.

**Fontenele `mmc2.xlsx` + NC_045512.2 FASTA + GenBank:**
- Review 13.1 s; Update 6.4 s; **Run 178.8 s**.
- **33 of 34 CSVs equal the local baseline**, with Python version the only difference.
- Verified in the staging downloads:

| Check | Hosted result |
|---|---|
| Variant observations and samples | 7,973 observations, 52 samples; QC paper preset keeps 7,961 |
| Depth / VAF | Depth from `TOTAL_DP`; alternative counts labelled `Reconstructed` |
| Wilson intervals | 7,970 of 7,970 counted observations, from reconstructed counts (stated) |
| Beta-binomial | **Empirical calibration skipped**: "No documented technical replicates … overdispersion cannot be estimated". Only the hypothetical ρ = 0.002/0.02/0.05 scenarios are shown, labelled "Not calibrated (sensitivity scenario)". **No fabricated empirical estimate** |
| Missing coverage | "Coverage-aware non-detection: Skipped. Unreported variants stay 'coverage unknown' (never 0%)" |
| Temporal | 52 of 52 samples dated; 474 mutation-site series |
| Yue–Clayton / PCoA | Full-count Yue–Clayton unavailable (no A/C/G/T counts). Exploratory Jaccard PCoA only, labelled as such |
| Reference verification | NC_045512.2 Verified compatible, 7,973 of 7,973 bases |
| Annotation | 5,677 of 5,677 SNVs annotated (12 of 12 CDS). Consequences: Missense 5,983, Synonymous 1,787, Stop gained 334, Noncoding 64, Start lost 14, Stop lost 4, Stop retained 2 (identical to local) |
| Mutation grouping | Mode `profiles` (default, no groups): "groups not shown by default because the method did not meet its validation criteria" |

**Missing optional metadata** (`mmc2.xlsx` alone, with no FASTA or annotation):
- Run 163.9 s.
- Readiness: 4 available, 5 partially available, 5 unavailable, 1 experimental, each with its reason.

**The compatibility report** explains every analysis that cannot run. Examples:
- "Needs individual sequences…"
- "Genome coverage known for 0 of 52…"
- "No per-position depth…"

## 7. Sprint 4 free literature search (hosted)

The ten-record H5N1 set was uploaded with enrichment on, then `/enrich/start` (7.0 s). The page states:

> "10 accession records searched. One publication found. No verified collection locations could be
> assigned to these samples. No metadata was changed. Ten accessions represent two distinct isolates…
> The publication is from 2022; the samples were collected in September 2026. … This article was
> published before these samples were collected and does not identify them. Its collection locations
> cannot be applied to your samples. Claude was not used…"

- **Older article rejected:** the **2022 publication (PMID 35302933) is not used** for the 2026
  specimens, as the goal requires.
- **Clear no-change result:** the answer appears in the first lines of the step, before the table.
- **Grouping:** accessions are grouped as two distinct isolates.
- **Back to the analysis:** "Continue with the original metadata" (`/enrich/run`) ran the genomic
  analysis in 14.2 s. Its main findings include "Sprint 4: … No metadata was changed".
- **Paid AI off:** shown as a clear status ("Claude was not used…"; system status "Disabled … No paid
  requests are possible"), not as an error.

**Not tested on staging:**
- Controlled HTTP 429 handling. It is covered only by the existing unit tests; deliberately
  rate-limiting NCBI was not attempted.
- Live Claude inference. **This is not validated.**

## 8. Local versus hosted comparisons

Same script, same inputs. Local: Flask on 127.0.0.1, GeoNames unset, limits as in the Dockerfile.

| Scenario | CSVs equal in content | Differences |
|---|---|---|
| Sprint 2 synthetic | 4 of 11 | GeoNames coordinates (staging has GeoNames; local does not) and the Python version. No change to strict pairs or clusters |
| Sprint 2 NCBI ten-record | 10 of 11 | Python version only |
| Integrated example | 43 of 44 | Python version only |
| Fontenele mmc2 | 33 of 34 | Python version only |
| Formula-safety table | 29 of 30 | Python version only |

- **Byte-level differences** (md5) also occur in PNGs, XLSX and PDF, because they carry rendering
  details and creation dates. The XLSX and PDF contents were checked separately (§9).
- **MAFFT** is 7.505 on staging and 7.526 locally; no output differed.
- **Main findings** are identical except the GeoNames-dependent Sprint 2 synthetic line.

## 9. Download validation (hosted)

| Check | Result |
|---|---|
| CSV / PNG / XLSX / PDF / JSON | All links returned 200. Fontenele: 56 files; integrated: 72 |
| XLSX equals CSV | 17 of 17 sheets equal the staging CSVs cell for cell (Fontenele) |
| PDF | 10 pages, 0.73 MB. All 5 main findings present, plus "overdispersion could not be estimated", "Experimental" and "coverage unknown" |
| Mapping-profile JSON | Downloaded; identical to the local one (ignoring creation time); contains no records |
| Formula safety | A `=HYPERLINK(...)` sample ID gives **0 formula cells** in the XLSX and **0 CSVs** with an unescaped formula |
| On-demand XLSX / PDF time | Fontenele XLSX **68.6 s**, PDF 9.4 s; integrated XLSX 2.1 s, PDF 8.6 s. No crash |
| Isolation and lifecycle | A 15-minute-old results link: **404** (10-minute lifetime). An unknown results ID: 404. Path traversal `/examples/..%2fapp.py`: 404. Links use 128-bit random IDs; they are not tied to a browser (audit P3) |

## 10. Runtime and memory observations

| Run on staging | Time | Local (same code, full CPU) | Local imitation of Render (512 MB, 10% CPU) |
|---|---|---|---|
| Fontenele + FASTA + GenBank, Run | **178.8 s** | 12.5 s | 123.5 s |
| Fontenele alone, Run | **163.9 s** | 11.6 s | — |
| Fontenele, first of two runs at once | 170.3 s | 12.0 s | — |
| Fontenele review (upload) | 13.1 s | 1.0 s | — |
| Integrated example, Run | 71.2 s | 14.6 s | — |
| Fontenele on-demand XLSX | 68.6 s | 4.5 s | 39.8 s |

**Interpretation:**
- **Slower than the imitation.** The staging CPU is about 40% slower than the local "10% of one core"
  imitation used in the audits.
- **The 180 s figure is not a request limit (inferred).** It is gunicorn's `--timeout`, and with the
  `gthread` worker this is a heartbeat check, not a per-request limit: gunicorn's documentation says
  that for non-sync workers it "is not tied to the length of time required to handle a single
  request".
  - So the 178.8 s run is not at risk of being killed by gunicorn, and runs of 164–179 s did complete
    (observed).
  - Whether Render's proxy allows a longer request was not tested.
  - The earlier reports' "under the 180 s timeout" wording overstated this as a hard limit.
- **Usability concern.** A visitor waits **about 3 minutes** for the Fontenele analysis with only the
  browser's loading indicator. There is no progress display.
- **Memory on staging was not observed;** the Render Metrics graph is **Requested**. Under the local
  512 MB imitation, the same run peaked at 302 MB.

**Schumann Berlin35 (local imitation only; not run on staging, as instructed):**
- Review peak **256 MB**.
- Run disabled, with the message "8,862 variants x 109 samples is 965,958 combinations, and this
  server can process 700,000. No samples or variants are removed automatically. Run the analysis on a
  local computer…".
- Upload plus confirmation took 189 s at 10% CPU, so on staging the refusal would take several
  minutes to appear.

## 11. Error recovery (hosted)

| Test | Result |
|---|---|
| FASTQ upload | 200, "Raw sequencing reads cannot be analyzed as processed results." |
| 1.6 M-cell CSV (3.5 MB) | 200 in 0.2 s: "larger than this server can process safely (about 1,600,000 cells; the limit is 1,500,000)". Refused before parsing |
| 11 MB upload | 413: "Upload too large. Files in one request must total less than 10 MB." |
| After the refused uploads | `/health` 200, and the next analyses ran normally |
| Two runs at once | The second was answered in 11.1 s: "Another analysis is running… Your files are still here". The first completed (170.3 s) |
| Worker restart, cold start, behaviour after a restart | **No unplanned restart occurred:** one "Booting worker" line across 06:03–06:32. A deliberate restart and a cold start were **not tested**. Results are in memory by design, so a restart loses download links (as documented) |

## 12. Browser usability findings

**Scripted only.** These come from the hosted pages in the built-in browser at 1024×768, plus element
positions read from the page. **No human usability test was done.**

| # | Finding | Severity |
|---|---|---|
| H1 | Homepage: the **file picker is 1,778 px down** (2.3 screens) and the Review button is at 1,939 px. "How this works" (expanded) and the example cards come first. "Which files do I need?" (2,010 px) appears **after** the Review button | P2 (UX) |
| H2 | No progress display during 1–3 minute analyses; only the browser's loading state | P2 (UX) |
| H3 | Results page for the integrated example: **16,097 px (about 21 screens)**. The download bar (920 px) and main findings (1,433 px) are near the top, which is good. The Sprint 5 section starts at 11,662 px, reachable through the "Jump to" links | P3 |
| H4 | The Fontenele results HTML is 3.4 MB (embedded figures); slow on weak connections | P3 |
| H5 | The Sprint 4 findings page states the answer first ("No verified collection locations… No metadata was changed") with a short status strip. Clear | OK |
| H6 | A large-dataset refusal (Berlin35) appears only after a long review on the free CPU | P3 |

**Short checklist for a later human test** (5–8 people, about 30 minutes each):
1. Find where to upload and say which files are needed.
2. Run `mmc2.xlsx` and state what the page shows while waiting.
3. Find the samples kept after QC and one mutation's interval.
4. Download the workbook.
5. Explain "coverage unknown" and "Experimental".
6. Run the H5N1 Sprint 4 example and say whether any location was changed.

The full checklist is in audit §9.

## 13. Remaining audit issues (not resolved by this task)

| Issue | Affects the staging demonstration? | Misleading-science risk | Recommendation |
|---|---|---|---|
| P2-1 "UTR" gene-track band over ORF3a–N in the genome-position figure | **Yes, visible in the Fontenele demo figure** | Moderate: the figure suggests an untranslated region where genes lie; the tables are correct | Targeted fix before showing that figure: draw the track from the GenBank CDS |
| P2-2 Influenza gene "NA" blank in the XLSX | No (SARS-CoV-2 demo); yes for influenza tables | Low–moderate (blank gene name; the CSV is correct) | `keep_default_na=False` in `StoredTables.get` |
| P2-3 Missing `transl_table` defaults to code 1 without warning | No (SARS-CoV-2 and influenza use code 1) | Moderate for non-standard-code organisms | Add a warning |
| P2-5 Mapping-profile entries not validated (`header_row: -5` accepted) | No (profiles are opt-in) | Low–moderate | Validate as `read_file_options` does |
| P2-7 Download cache: at most 4 runs, evicted by other visitors | Possibly, with several simultaneous visitors | None (availability only) | Age-only eviction or a larger limit |
| Schumann Berlin35 memory | Refused with a clear message (verified locally) | None | Optional categorical storage (533 → 50 MB, identical CSV) |

**None of these is a security blocker,** and none changes a numerical result. P2-1 is the only one that
shows up in the planned Fontenele demonstration, so it should be fixed, or the figure explained,
before showing it.

## 14. Deployment-only fixes and evidence still requested

**No code or configuration was changed for the deployment.** The only additions are:
- `validation/hosted_smoke_test.py` (the test client);
- this report;
- the progress record.

None of these is in the deployed branch, and none has been committed.

**Evidence from the Render dashboard:**
- **Received from the user:** the build and deploy log (commit, MAFFT, dependencies, gunicorn startup)
  and the runtime log for 06:03–06:32 UTC.
- **Still requested:** the Metrics → Memory peak during the Fontenele runs (around 06:13–06:27 UTC).

## 15. Final recommendations

1. **Before a class demonstration:**
   - Fix P2-1 (gene track), or explain the figure.
   - Add a visible "this takes about 3 minutes" note or progress message for pooled analyses (H2).
   - Move the upload controls higher on the homepage (H1).
   - These are targeted UI changes, not a redesign.
2. **Collect the four dashboard items in §14** to close the commit and memory evidence.
3. **Keep production on `sprint4-integration` (`6e122ff`)** until the user decides to promote. When
   promoting, push `5a20fb0` to the production branch deliberately, and re-run
   `validation/hosted_smoke_test.py` against production.
4. **Later:** fix P2-2, P2-3, P2-5 and P2-7, then run a human usability test.

---

**HOSTED VALIDATION PASSED WITH LIMITATIONS**

**Passed (hosted, verified):**
- The staging service runs `5a20fb0`, confirmed by the Render build log.
- The Docker build succeeds, MAFFT 7.505 is installed, and gunicorn `gthread` starts on Render's port.
  The runtime log shows no errors and no worker restart.
- `/health` and the homepage work.
- Sprint 2 (synthetic and the ten-record H5N1 set) and Sprint 3 + 5 (integrated example, Fontenele
  with and without optional files) complete.
- Core outputs equal the local baseline in content: 119 of 130 CSVs are identical, and every
  difference is the Python version or GeoNames being configured.
- All download types work, and the XLSX equals the CSVs.
- The safeguards hold: no fabricated overdispersion, grouping off, coverage unknown never 0%,
  segments kept apart, the 2022 article rejected, no formula execution.
- Unsupported and oversized uploads get controlled messages.
- No paid Claude call was made.
- Production was not modified.

**Limitations:**
- The Render memory graph has not been seen.
- The local Docker test could not run.
- Restarts, cold start, Render's long-request behaviour and HTTP 429 were not tested on staging.
- Hosted runs take 2.5–3 minutes, with no progress display.
- The P2 issues remain, including the misleading gene-track label.
- No human usability test was done.
