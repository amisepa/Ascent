# 2026-09-25 — Insert regenerated figures into the revision (v13 → v14)

## Goal (user's words)
"continue from TODO.md, the revised figures script is done"

i.e. TODO "NEXT": compose the make_revision_figures.m outputs
(manuscript/figures/regen, finished 14:05) into Figs 4, 5, 6, 7, 10, 11, swap them
into the v13 clean docx → v14, check the text against the new figures, render
to PDF and look at every figure page.

## Outputs (manuscript/, git-ignored)
- Cannard and Delorme 2026 v14 Entropy.docx (clean, Google-Docs-safe) + .pdf (41 pp).
- figures/regen/Fig{4,5,6,7,10,11}_composed.png (3000 px wide, inserted at 5.0 in).
- Rebuild: `python manuscript/tools/compose_insert.py manuscript/figures/regen
  "<v13 clean>.docx" <tmp>.docx` then `python manuscript/tools/v14_text.py <tmp>.docx
  "<v14>.docx"`; PDF via Word COM with a guard that aborts if New-Object attaches
  to an existing WINWORD (session scratchpad render_pdf.ps1; two unrelated WINWORD
  instances, pids 9944/40820, were left untouched).

## Checked against the new figures (numbers from regen CSVs / panels)
- Fig 4: r = 0.81; exponent and offset posterior-dominant; raw alpha posterior,
  corrected alpha weak frontally — text matches. Sliding exp/offset unchanged.
- Fig 5: r = 0.84; exponent max right posterior, offset max right frontotemporal
  (Discussion "right-lateralized" holds); IC14 = highest exponent; IC14 alpha
  maps = posterior max + central minimum (plus a weaker frontal pole).
- Fig 6: HDA panel now shows 63/64 channels significant (one masked near Pz/POz),
  consistent with the text (Fp1 50 ch, P10 13 ch). Word comment 66 deleted.
- Fig 7: every rho in 3.3.1 matches fig7A/B_values.csv. Fixed one claim:
  ExSEnt(HD)–HigFracDim (−0.79) is the strongest *negative* spatial relationship
  (SampEn–FuzzEn is 0.88).
- Figs 10–11: t, g, channel counts already matched; frequencies now from the CSVs
  (Results 3.3.3 and Discussion 4.3): raw 1–30.25 Hz, PO4 8.75 Hz; AF7 26.5 Hz;
  AF7 17.75 Hz; corrected 4.5–23.75 Hz, PO3 9.25 Hz.
- Response letter v4 cites none of the changed values: no edit needed.

## Other edits in v14 (manuscript/tools/v14_text.py)
- Figure refs: "Fig. 1, A" → "Fig. 1A"; "Fig 8. A." etc. → "Fig. 8A"…; caption
  labels "Fig 10./11." → "Fig. 10./11.".
- Fig 6 caption: "channels outside significant clusters are masked (white)".
- Layout: Figs 1 and 6 paragraphs had firstLine = 2608 twips, pushing the 5-in
  image off the right page edge; set to 0. Supplementary Fig S1 5.5 → 5.0 in.
  All 12 figures now sit where v13 had them (checked image bboxes in the PDF).

## Tooling changes (manuscript/tools, ignored)
- compose_insert.py: PSD cluster cells at one scale, lone last cell centred (was
  stretched to double size); Fig 6 row gap collapsed; Fig 6b paragraph removed;
  positive figure-paragraph indents zeroed. Old version: compose_insert_v13.py.
- v14_text.py: the text/layout pass above.

## Found, not fixed in the manuscript
- PSD cluster curves (Figs 10–11 lower panels) plot t ± 1.96 at the peak channel
  under a "Power (dB)" label: plot_clusters gets `tvals` as its data map (GLM
  branch). Same in Figs 8–9 (labels "MSE"/"MFE"/"RCMFE"; ascent_group_analysis.m
  passes tvals too). make_revision_figures.m (tracked) now has a `sections` gate
  and passes A − B, so the curves show the EC−EO difference in dB (mean ± 95% CI)
  — not yet run; needs the MATLAB desktop (`sections = {'psd'}; make_revision_figures`).
- Figs 4–5 time course: right y-label "Corrected (log10 ratio)" touches the "70"
  tick of the panel above (ascent_plot_aperiodic.m:441); cosmetic.
- Blue text in v14 (pre-existing): 0000FF in 2.7 and Fig 1 caption C, 4A86E8 in the
  Intro toolkit paragraph — revision highlighting or leftovers? Ask the user.

## Round 2 (user, ~15:15): "yes phrasing review no uploading online for AI check ...
match my writing style" (Google Scholar papers); "I ran that command:
sections = {'psd'}; make_revision_figures its done. now lets do the other ones,
including multiscale run it if you need stats"; blue text = outdated tracking of
changes since v6 SUBMITTED (G:/My Drive/Ascent/old) → mark all changes vs v6 in red
(word level if feasible, else paragraphs, else explain to reviewers); check
Results/Discussion against the figures.

- Figs 10–11 recomposed with the new curves (15:12 run); captions now say what the
  heatmap, maps and curves show (v14_text.py). v14 docx rebuilt.
- Multiscale stats rerun headless (scratchpad ms_check.m, rng(1), serial stubs):
  current pull_clusters merges adjacent negative/positive regions into one cluster
  labelled by its peak sign (MFE-mean 2 clusters instead of 3; RCMFE-mean fine
  negative merged into the positive one). compute_mcc forms clusters on |t|, so
  mixed-sign significant regions exist. Old published CSVs
  (manuscript/figures/group_results) came from a column-wise sign split.
- Decision: report clusters separately by sign (cell-wise split of the significant
  mask, then pull_clusters per sign) — `pull_clusters_by_sign` local function in
  make_revision_figures.m and ascent_group_analysis.m (the group script also now
  plots EC−EO difference curves). Peaks/t/g identical to the text; extents change:
  RCMFE-sd 2–9 / 8–30; MFE-mean neg 2–4, pos 4–11 with 61 ch; RCMFE-mean pos 4–13,
  neg 12–30 with 48 ch (was 53). Raw PSD gains a 4th small cluster: Iz 2.5–4 Hz,
  2 ch, t −3.35, g −0.52 (scratchpad psd_split.m). Corrected PSD unchanged.
- Text for these: manuscript/tools/v15_text.py (Methods 2.8 sentence, Results
  3.3.2/3.3.3, Discussion 4.2/4.3, Fig 8 caption). Verify against regen CSVs after
  the user's run of `sections = {'aperiodic','multiscale','psd'}` (aperiodic = the
  Figs 4–5 right-axis label, now 'Corrected (log_{10})' in ascent_plot_aperiodic.m).
- Red markup vs v6: manuscript/tools/redline_vs_v6.py (word level; a word is kept
  black if it sits in a ≥3-token run found in one of the 3 most similar v6
  paragraphs; citations compared as placeholders; equations by symbol sequence;
  old blue removed). Trial on v14: 78% of words red, 84/100 equations red (v6
  equations had "lnln", "expexp", spelled-out "tau"; the rebuilt ones differ).
  Paragraph mode (--para) would be 90% red. Word level reads well → use it.
- Phrasing review: style profile agent (scratchpad/style_profile.md); section
  exports in scratchpad/review/*.txt with % new vs v6 per paragraph
  (manuscript/tools/export_review.py). Edits → manuscript/tools/phrasing_edits.json,
  applied by v15_text.py. Only revision-added text is rephrased (v6 text is the
  authors' own and would turn red).
- Letter: v4 has 17/14 tracked ins/del (earlier session). Needs a reply to the
  editor's "highlight revisions" item (not answered anywhere).

- Style profile (agent, 7 papers + v6; scratchpad style_profile.md): his 2020–22
  voice = mean 28-word sentences with high within-paragraph variation, hedges,
  "(i.e., …)", Hence/Thus/Note that, no em-dashes, no clause-joining semicolons,
  colons only before lists. AI tells in the revision (and in v6 Abstract, Intro
  1.1–1.3, early Discussion, Conclusions): even rhythm, ", enabling/yielding"
  tails, colon punchlines, triads, "rather than" (0 in his papers), abstract
  "across", demonstrate/reveal/support, evaluative closers, no hedges.
- Five review agents (hit the usage limit once, resumed): edits A 39, B 71, C1 44,
  C2 69, letter 43; all pass review/check_edits.py (numbers, citations and
  Fig/Section/Eq references identical in old/new). Merged into
  manuscript/tools/phrasing_edits.json (ids = paragraph index of the export) and
  letter_edits.json; my fixes in manual_edits.json + letter_edits.json tail:
  Discussion 4.2 said EC<EO meant "reduced entropy during EO" (fixed to EC);
  "multiscale entropy (mMSE)" → modified (Intro 1.3 and letter); Hedges' g
  parenthesis closed; P0233 "substantially greater demand" (84→105 s) softened;
  typography; plain Discussion headings (4.1–4.5 were claim-style); letter:
  Kamali 2026, MFE location (Section 2.2.2), statistics are external (L0022,
  L0095), τ wording (L0035); new reply to the editor's highlighting item.
- Not changed, for the user: many one-sentence Discussion paragraphs (profile
  prefers fewer, longer ones); ICA-isolates-sources claim repeated in 4.3/4.4.
- One-shot build: `bash manuscript/tools/build_v15.sh <scratch> [outdir]` →
  v15 clean, v15 changes-in-red (vs v6), v15 edits since v14 (tracked), letter v5
  clean + tracked. Dry run OK (250 edits, letter 50/50, 79% red, 16,400 words).

## Round 3 (user, ~19:05): "MATLAB run is done, fix the discussion and add anything
that seems like it's missing, address the claim that 'ICA isolates individual
sources' appears in both 4.3 and 4.4, check all references are formatted
properly, and build the final files, push commits." (+ MATLAB warnings, looked hung)

- The run was not hung: 'progress' was off and the time-resolved aperiodic fit
  is silent for ~10 min per domain; 3 HBN ICA batch jobs from
  eeglab_source_localization_2026 (started 12:58, ~12 cores) slowed it. Outputs
  19:10–19:23. make_revision_figures.m now prints timestamped progress, turns
  ascent_compute progress on, silences name-conflict warnings while adding EEGLAB
  and drops plugin compat/legacy folders from the path.
- Regen CSVs match the v15 text exactly (all 8 summaries; PSD raw 4 clusters incl.
  Iz 2.5–4 Hz, t −3.35, g −0.52). Figs 4–5 label fixed; Figs 8–9 recomposed in the
  old A/B/C layout with EC−EO curves.
- Discussion (v15_text.py + manual_edits.json): 4.3 stray repeat paragraph folded
  into the group-PSD paragraph (posterior alpha peak PO3, "misleading inferences"
  sentence); 4.4 ICA statement rewritten once, hedged (dipolar maps only; complement,
  not replace); duplicate in 4.5 deleted; 10 one-sentence paragraphs merged;
  added limitations (short recordings → coarse-scale entropy less reliable [18];
  each measure tested separately, correction within not between measures →
  illustration, not confirmatory) and a future-directions note on removing the
  evoked response before post-stimulus fits (checked in Gyurkovics 2022 PDF);
  heading 4.3 "Aperiodic parameterization".
- References: 64 DOIs checked against Crossref (scratchpad crossref.json). Fixed:
  [7] title case, [11] "Circuits Syst. I", [16] article no. imag-2-00054, [17]
  journal abbreviation, [18] article 138, [22] Waschke (Crossref has the preprint
  typo "Washcke"), [23] now 10, 1431–1440, [29] NeurIPS 8 chapter pp. 145–151
  (last PDF page headed 151), [33] H2039–H2049, [37] Phys. D, [42]/[46] location,
  dates, DOIs, [43] Phys. A, [55] den Nijs, PLoS ONE/PLoS Comput. Biol. unified.
  In-text citations: all 68 cited in first-appearance order, no duplicates.
- Red markup: references are marked per entry (36 of 68 new since v6), body text
  word level: 83% of words red; equations 84/100 red.
- Final build 19:3x: manuscript/ v15 clean/red/tracked docx + clean/red PDFs
  (41 pp; every figure page checked, nothing overflows); letter v5 clean/tracked
  (manuscript 16,500 words in the editor reply).
- Committed and pushed: code (make_revision_figures.m, ascent_group_analysis.m,
  ascent_plot_aperiodic.m) and notes. Not committed: figures/ deletion (README
  embeds figures/figure2.png, fig3.png, figure4.png), papers/ renames, IDEA.md.
- TODO "AI-generated check" (not started). Quick counts, v14 body 17.3k words:
  em-dash 25; "highlight" 8, "additionally" 3, "moreover" 2, "notably",
  "underscore", "comprehensive", "compelling", "strikingly" 1 each ("robust" 16 is
  mostly technical). Letter: 5.2k words, no em-dashes.
