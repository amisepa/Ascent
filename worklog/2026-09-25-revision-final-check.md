# 2026-09-25 — Final pre-submission check of Entropy revision (v11 → v12)

## Goal (user's words)
"I want to submit the revisions today. Go through the whole paper and reviewers
comments and see if we're missing anything. Pay special attention to equations
(both referenced and in the text) are in right format and correct. Check the
code again if needed. Check that figures are latest, and that the text
describes them properly. Then that the methods/results are described properly
and corresponding interpretations are correct. Finally update the conclusions,
abstract, and intro if needed. Keep the journal template formatting."

## Inputs
- Manuscript: ~/Downloads/Cannard and Delorme 2026 v11 Entropy - Revisions(4).docx
  (newer than the manuscript/ copy of the same name)
- Reviewer comments + responses: ~/Downloads/Revision 1 - Reviewer comments v2(1).docx

## How the edits were made
- All edits applied through Word COM as tracked changes (author
  "Claude (pre-submission check)") on a copy, so the MDPI template is untouched.
  Scripts (session scratchpad): word_ops.ps1 (runner), build_ops.py + edits*.py
  (edit list), renumber.py (citations/refs), build_rev.py (response letter).
- Citation renumbering and the reference-list rebuild were applied UNtracked
  (mechanical); they assume all tracked changes are accepted.
- Equations: all 34 rebuilt from UnicodeMath, number right-aligned "(n)"
  instead of the "Equation (n):" prefix.

## Verified against code / data (numbers and where they came from)
- Aperiodic sample-subject numbers: headless rerun of current code
  (scratchpad ap_full.m) reproduces the text exactly: channel r = 0.901,
  exp 1.045±0.540, offset 1.137±0.430, alpha 5.76 vs 4.34 µV²/Hz (24.6 %
  aperiodic); sliding exp 1.352±0.478, offset 1.824±0.462 (SD over channels AND
  time bins); ICA top IC14, r = 0.884, sliding 0.760±0.411 / −0.549±0.682;
  ICA alpha 49 % aperiodic (mean across ICs).
- Group stats: serial rerun of eeg_robust_statistics on the cached outputs
  (~/Documents/biosemi_data/ascent_outputs_*_new.mat) reproduces every reported
  single-scale cluster and MSE/MFE/RCMFE cluster with CLUSTER-BASED correction
  (mcc_type 2). TFCE was never used -> Methods/captions/response fixed.
- ExSEnt HDA: ascent_group_analysis.m line 549 (compute_mcc for ExSEnt3) is
  commented out, so HDA reused the HA mask. Correct result: 2 positive
  clusters, 63 channels (Fp1 50 ch t = 8.83 g = 1.37; P10 13 ch t = 4.60 g = 0.71).
  Fig 6 HDA panel still shows the wrong mask -> user must regenerate.
- Group data (computed Apr 29 / May 18) used: MSE/MFE with per-scale
  z-scoring (old compute_SampEn/FuzzEn behaviour), RCMFE fixed r, RCMFE-mean
  includes scale 1, PSD 2-s windows, P_max 6. Current defaults differ (zNorm = 0,
  P_max = 3, robust fit) -> documented in 2.8 and Discussion 4.2.
- Recording length per condition 62–80 s (mean 68 s), not 2 min.
- EC/EO preprocessed .set spectra both show a ~50 Hz low-pass roll-off (no
  filtering confound); Methods lists only a notch -> Word comment left for user.
- specparam validation (.mat recomputed by agent): ASCENT exp recovery 0.946 vs
  specparam 0.938 (text had them swapped); offset recovery 0.52 vs 0.97
  (now disclosed); specparam version 2.0.0rc6.

## Environment notes (this machine)
- MATLAB graphics (drawnow) hang under `matlab -batch` in this session, and a
  second parfor call hangs; headless numeric runs work with serial stubs
  (scratchpad/stubs). Figures must be regenerated from the MATLAB desktop.
- Hung processes left running (hook blocks killing): several MATLAB.exe from
  08:12–08:40 and a WINWORD /Automation from 08:29.

## Outputs (manuscript/)
- Cannard and Delorme 2026 v12 Entropy - Revisions (tracked check).docx
- Cannard and Delorme 2026 v12 Entropy - all changes accepted.pdf
- Revision 1 - Reviewer comments v3 (tracked check).docx
- figures/Fig1_patched.png, Fig4_new.png, Fig5_new.png
- Word COM renumbering pass hung (>1.5 h); renumbering, reference rebuild,
  caption bold fixes and response-letter edits were done by direct XML edit
  (scratchpad phase2_xml.py, rev_xml.py, docx_xml.py). Refs 1–68 verified cited
  in first-appearance order; Garrett 2014 dropped (was only mis-cited).

## Unfinished / for the user
- Regenerate Fig 6 (HDA panel) after restoring line 549.
- Optionally regenerate Figs 4–5 in the MATLAB desktop (inserted panels are the
  Sept 3 rerun; time-course right-axis label overlaps; raw legend missing).
- Fig 1A screenshot predates the ICA-domain selector (optional re-shoot).
- Confirm low-pass/notch wording in 2.8; λ default (0.001 via ascent_compute vs
  0.01 in compute_ExSEnt).
- Abstract ~230 words (MDPI guideline ~200).
- Plugin code changes described in the paper (P_max 3, robust fit, pruning,
  knee correction, zNorm) are uncommitted in the working tree.

## Round 2 (afternoon)
- Google Docs import: tracked edits inside equations and Word "#(n)" numbering
  don't survive -> deliver a clean copy (v13) with equations followed by a
  tab + "(n)" in plain text, math style normalized (Word had stored bold-italic).
- specparam port fixed by agent (compute_AperiodicFit.m): robust initial fit
  now specparam's; ASCENT vs Python r > 0.999 (was offset r 0.45 synthetic);
  ground-truth offset recovery 0.52 -> 0.97. Committed + pushed 5002f81.
- Subject-level rerun (fixed code): channel r(exp,off) 0.901 -> 0.812, offset
  1.055±0.327, alpha 22% aperiodic; ICA IC14, r 0.884 -> 0.844, 47%; sliding
  unchanged to 2 dp.
- Group aperiodic recomputed with current code/defaults (scratchpad
  group_aperiodic_rerun.m; written into the cached biosemi_data .mat files,
  originals kept as *_backup_20260925.mat): exponent 52 ch (AF7 t 5.63),
  offset 61 ch (P8 t 5.89); raw PSD 1-30 Hz +, 20.75-40 Hz − (AF7), 17.5-18 −;
  corrected PSD now 2 clusters: 4.5-24 Hz + (PO3 9 Hz), 34.75-40 Hz − (POz).
  Fig 7 values recomputed (old recomputation reproduces published values).
- Filtering: data show a 50 Hz low-pass in both conditions -> Methods updated.
- GLM aperiodic project (~/Documents/MATLAB/eeg_glm_aperiodic_covariate)
  cross-check: EO/EC fitted-aperiodic changes are leakage/artefact-prone;
  spatial r ~ 0 is not evidence of dissociation; ratio correction = multiplicative.
  Discussion 4.1/4.3/4.5/4.6 and Methods 2.4 caveats added. Kałamała now
  Psychophysiology 2026.
- Figures 4, 5, 6, 7, 10, 11 need the MATLAB desktop: make_revision_figures.m ->
  manuscript/figures/regen; scratchpad compose_insert.py composes and swaps them.
