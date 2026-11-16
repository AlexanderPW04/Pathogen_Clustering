# Homepage + Sprint 4 UI — report

**Branch:** `feature/homepage-sprint4-ui-20261116` (local only).

**Based on:** `145dc4c`, whose code equals the staging-validated `8690561`.

**What was not done:** nothing was pushed, merged, or deployed; no paid API call was made, and
`AI_ENRICHMENT_ENABLED=false` was used throughout.

**Day-by-day log:** [HOMEPAGE_SPRINT4_UI_PROGRESS.md](HOMEPAGE_SPRINT4_UI_PROGRESS.md).

## 1. Homepage improvements

**Header**
- Title plus the requested sentence: "Analyze individual pathogen sequences, pooled sequencing data,
  and optional publication-based metadata."
- Then the screening-signal label and "No files yet? Try a synthetic example".

**Three analysis cards directly below**

| Card | Workflow strip |
|---|---|
| **Individual Sequences** · Sprint 2 | Upload → Align → Compare → Results |
| **Pooled Sequencing** · Sprint 3 + 5 | Upload → Validate → Analyze mutations → Results |
| **Combined Analysis** · Sprint 2 + 3 + 5 | Process both sources → Compare compatible variants → Results |

- Each card shows its one-line purpose.
- The whole card selects its radio button: the label is stretched over the card in CSS, so there is
  still one native radio group, keyboard arrows work, and the summary is linked with
  `aria-describedby`.
- "Not sure? Detect from my files" stays as a smaller fourth choice and remains the default, so
  automatic, content-based file detection is unchanged.

**Upload fields**
- Only the chosen analysis's fields are shown. This reuses the existing `MODES`, `UPLOAD_FIELDS`, and
  switching script; there is no new backend selection.
- In Combined mode the fields are grouped under *Individual sequences*, *Pooled sequencing results*,
  and *Reference genome and gene annotation*.
- Multiple files, the content checks, field detection, and the user-confirmed mappings on the review
  page are unchanged.

**"What is…?" help**
- Every field has a native `<details>` question with the requested one-line answer, for example "What
  is a processed SNV table?".
- The earlier detailed text sits under a nested "Technical details".

**Kept, but moved below the upload**
- **Sprint 4:** remains an optional step ("3 · Optional: publication-based metadata"), not a fourth
  analysis.
- **Collapsed or moved down:** the long file guide (one collapsed section), the synthetic examples,
  and "How this works".

**Measured positions** (built-in browser, with the 32 px preview banner included)

| Viewport | Cards start | "Add your files" | Horizontal scroll |
|---|---|---|---|
| 1280×800 | 309 px | 563 px | none |
| 1024×768 | 309 px | 587 px | none |
| 375×812 | 433 px (85 px of it is the preview banner) | — | none |

The step pills were then made smaller.

## 2. Sprint 4 workflow improvements

**Step 1 — "Find and Check Publications"**
- One prominent button, on the file review and on the Sprint 4 page.
- It runs the existing free NCBI/PubTator discovery, retrieval, and deterministic evidence check.
- The existing banner appears at once, with an honest message, an elapsed clock, and no percentage,
  and blocks repeated clicks.
- The search is one request and the backend reports no stages, so no stage is invented.

**Step 2 — conclusion first**
- The box shows, in order:
  - the scope ("10 accession records searched.");
  - the headline;
  - one sentence about Claude;
  - the next action;
  - a short "why" list.
- After a search, the page scrolls the conclusion into view and focuses it.
- **H5N1-shaped case (the 2022 article and the September 2026 accessions):**
  - headline: "One publication was found and reviewed. It predates the uploaded samples and does not
    provide verifiable collection locations for them. No metadata was changed."
  - Claude sentence: "Claude was not called because no eligible sample-specific evidence was found."
- **When "reviewed" is used:** only when every publication was retrieved with text and assessed. A
  failed retrieval reads "…could not be retrieved … The evidence review is incomplete. This does not
  mean the literature has no evidence."
- **Claude next step:** when sample-specific evidence exists, a "Next: check the eligible passages
  with Claude (the cost is shown first)" link jumps to the existing Claude step. Cost approval,
  acknowledgement, and access-code checks are unchanged.

**One compact card per publication**
- Title, year, and identifier.
- Five answers:

| Question | Possible answers |
|---|---|
| Geographic information found? | Yes / No / Not checked |
| Linked to uploaded samples? | Yes / No / Uncertain |
| Metadata improvement supported? | Yes / No / Needs Claude review / Not checked |
| Claude analysis performed? | Yes / No / No (offline demonstration) |
| Changes applied? | Yes / No |

- A place named in the article is marked "general mention only" when it is not tied to a sample.
- Publication details (plain explanation, places, dates, retrieval, link counts) are expandable.
- Accession mappings, study relationships, and provenance stay in "Technical details".

**Step 3 — sticky bar, always visible**
- "← Back to file review" and **"Continue to Genomic Analysis"**.
- With approved proposals, Continue applies them to a copy of the inputs; otherwise it runs with the
  original metadata. The hint states which.
- "Skip enrichment (use original metadata)" is offered when proposals exist.

**Less scrolling and fewer clicks**
- The per-sample table and the publication details are collapsed.
- Going back from the file review reopens the same findings without a second search (existing
  behaviour, still tested).
- No additional NCBI or Claude requests were added.

**Safeguards untouched**
- Sample-to-publication matching, publication-date validation, geographic verification, exact
  quotation checks, human approval, GeoNames precision, the paid-call authorization, NCBI rate limits,
  and caching.
- Only `literature/review.py` (summaries) and templates changed. `literature/eligibility`,
  `verification`, `discovery`, `workflow`, `apply`, and `http_client` are unchanged.

## 3. Changed files

| File | Change |
|---|---|
| `app.py` | Mode card text (`title`, `sprint`, `summary`, `flow`), `UPLOAD_GROUPS`, a `group`/`question`/`answer` for each upload field, mode labels naming Sprint 5, scope prefix of the results-page note, `applied_facts` kept after applying |
| `literature/review.py` | "Not checked"; `changes` and `text_checked` for each publication; headline rules ("reviewed" only when read; requested H5N1 wording); `claude_line`; `searched`; the offline demonstration is never shown as a Claude analysis |
| `templates/_upload.html` | Mode cards with workflow strips, grouped fields, "What is…?" help, collapsed file guide, examples row |
| `templates/_literature.html` | Conclusion box, publication cards, collapsed per-sample table, sticky Back/Continue bar, renamed buttons |
| `templates/index.html` | Header text, CSS for cards, strips, groups, and publication cards, step pills on phones, "How this works" moved below the upload on the homepage, loading messages, scroll to the conclusion |
| `templates/_toolbar.html` | `id="back-to-review"` so the sticky bar can use the same Back form |
| `validation/local_ui_preview.py` | **New.** Local-only preview server (see §5) |
| `tests/test_homepage_sprint4_ui.py` | **New.** 6 tests (§6) |
| `tests/test_final_polish_navigation.py` and six other test files | Assertions updated to the intended wording and order |
| `docs/HOMEPAGE_SPRINT4_UI_PROGRESS.md`, `docs/HOMEPAGE_SPRINT4_UI_REPORT.md` | **New** |
| `docs/ui_preview_screenshots/` | **New.** 11 screenshots (§4) |

No new package or executable dependency was added. The screenshots used the locally installed Chrome
and Pillow from the existing environment; neither is needed by the application.

## 4. Screenshots

All screenshots are in [`docs/ui_preview_screenshots/`](ui_preview_screenshots/) and come from the
**synthetic** local preview, which shows a purple banner.

**Homepage**
- `home_1024x768.png`: the default homepage.
- `home_individual_1024x768.png`, `home_pooled_1280x800.png`, `home_combined_1280x800.png`: each mode
  chosen.
- `home_combined_mobile_375.jpg`: the phone layout, from the built-in browser's mobile emulation.

**Sprint 4** (full section, 1024 wide)
- `sprint4_before_1024x768.png`: the file review with the Find and Check Publications button.
- `sprint4_one_publication_1024.png`: the H5N1-shaped case.
- `sprint4_general_mention_1024.png`
- `sprint4_no_publication_1024.png`
- `sprint4_eligible_1024.png`: the Claude step offered.
- `sprint4_retrieval_error_1024.png`

**Not saved as images**
- **Proposals awaiting approval:** checked live (§5), not saved as an image.
- **Sprint 4 at phone width:** checked by measurement. There was no horizontal scroll, the conclusion
  scrolled to the top and received focus, and the Back/Continue bar stayed at the bottom (700–804 of
  812 px). The browser pane could not draw a screenshot while hidden.
- **Headless "375 px" images:** discarded, because headless Chrome keeps a window at least about
  500 px wide and clipped them.

## 5. Opening the local preview

```bash
python validation/local_ui_preview.py
```

Then open <http://127.0.0.1:5058/preview>. You can set `PREVIEW_PORT` to use another port.

**Links on the preview page**
- **Homepage:** the four modes, each opened through the real Back-to-upload route with that mode
  chosen.
- **Before searching:** the real file review with Sprint 4 on.
- **One publication (H5N1-shaped):** predates the samples, nothing changes.
- **General mention:** places are named, but not for these samples.
- **No publication found.**
- **Eligible evidence:** press *Analyze Literature with Claude AI* to see the cost-approval screen,
  then *Run offline demonstration* to see the proposal awaiting your approval.
- **Retrieval error.**

**How it works**
- It stages `tests/fixtures/ncbi_virus_structure_*` (ten invented influenza accessions collected in
  September 2026) through the real upload code.
- It replaces only the publication source with invented in-memory articles from the test suite.
- The session is marked synthetic, so the offline mock provider answers. Claude, NCBI, PubTator, and
  GeoNames are never contacted.
- "Continue to Genomic Analysis" runs the real analysis (MAFFT must be on `PATH`). This was checked
  end to end: approving the synthetic proposal produced "1 of 1 proposed improvement(s) approved;
  1 value(s) used in this analysis".

## 6. Test results

| Run | Result |
|---|---|
| Clean baseline, separate worktree of `145dc4c` | Ran 714 tests — **OK (skipped=16)**. The earlier record of 4 skips came from runs with the optional real-data variables set (`FONTENELE_MMC2`, `CRITS_CHRISTOPH_DIR`, `SCHUMANN_DIR`, `ANNOTATION_MPOX_GENBANK`, and the like) and MAFFT on `PATH`. None of these were set here, so the same 16 tests skip on the baseline and on this branch. |
| Targeted runs during the work | `test_sprint4_interface`, `test_literature_review`, `test_literature_eligibility`, `test_upload_enrichment`, `test_literature_web`, `test_literature_enrichment`, `test_ui_navigation`, `test_mode_upload`, `test_final_polish_navigation`, `test_final_polish_feedback`. All 19 failures after the change were wording or order assertions, and they were updated. |
| New `tests/test_homepage_sprint4_ui.py` | 6 tests, OK. They check: "Changes applied?" is Yes only after a value was applied; an offline demonstration is never reported as a Claude analysis; "reviewed" is used only when the text was read; a chosen mode shows only its field groups; each field has its question and technical details; every preview state opens and is labelled. |
| **Full suite, once at the end** (main checkout, same environment as the baseline) | Ran **720** tests in 631.6 s — **OK (skipped=16)**. That is 714 + the 6 new tests. No failures, no errors, and the same skip set as the baseline. |

## 7. Scientific calculations unchanged

- **Code:** no file under `genomics/`, `wastewater/`, or `offline_workflow/` changed, and neither did
  `compare_genetic_distances.py`, `align_sequences.py`, `detect_potential_clusters.py`,
  `add_temporal_geographic.py`, `geocode_locations.py`, or the literature verification, eligibility,
  discovery, apply, and workflow modules.
  - Not touched: MAFFT, VAF, beta-binomial, overdispersion, clustering, grouping, annotation,
    thresholds, and metadata rules.
- **Outputs:** the integrated example and the Sprint 2 example were run in the baseline worktree and
  on this branch. **All 52 CSV downloads are byte-for-byte identical** (SHA-256).
- **What did change:** the wording of Sprint 4 headlines and statuses. Every statement describes the
  same checks as before. Two are stricter than before: "reviewed" now requires the text to have been
  read, and the offline demonstration is no longer reported as a Claude analysis.

## 8. Remaining usability issues

1. **Default choice:** "Detect from my files" is still the default, so no card is pre-selected. This
   keeps the existing behaviour and tests; a first-time user may still need to notice that the cards
   are choices.
2. **"How this works" placement:** it moved below the upload only on the homepage. On the review,
   Sprint 4, and results pages it still sits above the content (about 70 px).
3. **Two synthetic notices:** synthetic Sprint 4 pages show two SYNTHETIC notices, the page-level one
   and the Sprint 4 one (existing behaviour).
4. **Step-pill numbering:** "3 · Review available analyses" is never highlighted, and "Run · Metadata
   enrichment" is step 4 while results are step 5. This existing numbering may still confuse.
5. **Search progress:** the search reports only an elapsed clock. The backend does not expose stages,
   and no job queue was added (by design).
6. **Claude step placement:** the Claude step card sits below the publication cards; the "Next" link
   jumps to it.
7. **Combined workflow strip:** it wraps onto two lines at 1024–1280 px, and the arrow then starts the
   second line.
8. **Card focus outline:** it uses `:has()`. In a browser without it, only the native radio focus ring
   shows.
9. **Missing image:** there is no saved image of the "proposals awaiting approval" screen or of
   Sprint 4 at phone width (both were checked live or by measurement).
10. **Not on staging:** nothing was deployed, so this UI has no hosted validation yet.
11. **Usability and live Claude:** no human usability test was done, and live Claude NLP remains
    unvalidated. No paid request was made.

## Verdict

**HOMEPAGE + SPRINT 4 UI READY FOR REVIEW**

It is local only, on `feature/homepage-sprint4-ui-20261116`. It still needs your review in the local preview, a human usability check, and, if approved, a separate staging push and hosted revalidation.
