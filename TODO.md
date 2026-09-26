# ASCENT — TODO

## Status 2026-09-25 (supersedes older notes below where they conflict)
- DONE: knee-mode aperiodic correction (full Lorentzian removed) — item 2.
- DONE: compute_AperiodicFit now ports specparam 2.0 exactly (robust initial
  fit, peak search in Hz with edge/overlap drops, plain final refit). Validation
  vs Python specparam (validation_specparam/matlab_results_v2.mat): r > 0.999 for
  exponent and offset, real and synthetic; identical peak counts. Old 2c numbers
  are obsolete.
- DONE: ascent_group_analysis.m ExSEnt HDA now gets its own compute_mcc (was
  reusing the HA mask).
- DONE 2026-09-25: v14 docx + PDF (manuscript/) with the regenerated Figs 4–7,
  10, 11, text checked against them, Fig 6 comment removed, Figs 1/6 no longer
  overflow the page. Details: worklog/2026-09-25-insert-regenerated-figures.md.
- DONE 2026-09-25 evening: v15 (manuscript/): clean, changes since v6 SUBMITTED
  in red, edits since v14 (tracked), PDFs; letter v5 clean + tracked. Figures
  4–11 regenerated, phrasing review, Discussion restructured, references checked
  against Crossref. Rebuild: `bash manuscript/tools/build_v15.sh <scratch>`.
  Details: worklog/2026-09-25-insert-regenerated-figures.md.
- NEXT (user): import v15 clean into Google Docs as the new master; submit the
  red version + letter v5; delete the old blue-tracked Google Doc copy.
- OPEN (eeg_robust_statistics repo): pull_clusters' 2D path labels clusters on
  the unsigned mask, so adjacent positive/negative regions merge; ASCENT scripts
  now use a local pull_clusters_by_sign. Consider a sign-split option there.
- OPEN: knee mode still differs from specparam (offset MAE ~0.5) — starting
  guesses/bounds of the knee fit.
- OPEN: ascent_compute passes no tau to compute_SampEn (tau always 1);
  'slidPeakWidthLimits' name-value never matches (lowercased key);
  sd/var coarse-graining drops scale 2 twice in ascent_compute:503-506 plots.
- OPEN: cluster overlay default (item 1).

# Older notes (Sept 2026 revision session)

## 1. Cluster overlay on the main multiscale heatmap (ascent_plot_multiscale.m)
- User is unsure about the recent automatic cluster detection that decides what
  gets plotted in the scalp topo + curve plots (manuscript Fig. 2).
- Last time it required manual threshold adjustment to look right; user thinks
  there should probably be only ONE cluster shown, not up to 3.
- IDEAS PROPOSED (do not implement without asking):
  a) Show only the top-1 cluster by default.
  b) Add an interactive SLIDER on the figure so the user can move the
     percentile threshold live; the topo and curve panels update dynamically
     as the slider moves. (uicontrol 'style','slider' on the figure, callback
     re-runs find_entropy_clusters + redraws the topo/curve subplots.)
- Related files: ascent_plot_multiscale.m (find_entropy_clusters,
  clusterThresh opt, default 0.75), ascent_plot.m ('ClusterThresh' pass-through).

## 2. Knee mode + aperiodic subtraction (audit findings, awaiting decision)
- Current code disables aperiodic correction in knee mode (ascent_compute.m,
  compute_AperiodicFit_sliding.m, GUI callback greys out the checkbox).
- AUDIT RESULT (Sept 2026): upstream specparam/FOOOF has NO such restriction.
  The correction "log10 P - (b - chi*log10 f)" as currently written is only
  valid in fixed mode, BUT the correct knee-mode analogue exists and is used
  upstream: subtract the full aperiodic MODEL curve (Lorentzian, evaluated
  with the fitted offset/knee/exponent), i.e.
      psd_corrected = psd ./ 10.^(offset - log10(knee + f^chi))
  FOOOF computes the flattened (aperiodic-removed) spectrum internally on
  EVERY fit regardless of mode (_spectrum_flat = power_spectrum - _ap_fit,
  specparam algorithm.py ~line 186; same in FOOOF 1.x fit.py).
- RECOMMENDATION: enable correction in knee mode by subtracting the full
  knee-model curve (compute_AperiodicBandPower needs a knee branch), rather
  than the current fixed-only linear formula. Knee-mode alpha_osc split needs
  the same treatment (currently NaN in knee mode).
- DECISION PENDING with Cedric before touching code.

## 2b-IMPLEMENTED (Sept 2026): P_max default 6 -> 3; SPRiNT peak pruning
- P_max lowered to 3 per Kałamała et al. 2025 (reliability decreases with
  allowed peaks; their best-reliability caps were 0-3 peaks; SPRiNT's
  empirical EEG settings use max 3). Changed in: ascent_get_params.m
  (default + doc), ascent_compute.m (doc), ascent_compute_gui.m (GUI field),
  compute_AperiodicFit.m (default + doc).
- SPRiNT outlier peak pruning implemented in compute_AperiodicFit_sliding.m
  as Step 5 (default ON): peaks with <3 similar neighbours (within 2.5 Hz,
  6 time bins) are removed; aperiodic model refit at pruned bins via OLS on
  the peak-removed spectrum; corrected PSD regenerated; peaks_t updated.
  Toggles: 'PrunePeaks' (true), 'PruneMinPeaks' (3), 'PruneFreqTol' (2.5),
  'PruneTimeBins' (6) - SPRiNT's recommended values.
- ALSO FIXED latent bug: peaks_t was collected from info_win.peaks, a field
  that never existed in compute_AperiodicFit (correct name: peak_params) -
  peaks_t had silently been empty. Now collects info_win.peak_params.
- ascent_plot_aperiodic.m: removed stale src/srcLabel vars (dead since the
  yyaxis fix); corrected traces now read psd_corr_t directly; fixed
  try/catch alignment lint warning.
- NOT YET RERUN - batch with other pending changes per user instruction.

## 2b. Literature cross-check of aperiodic pipeline (Sept 2026, all READ)
- Kałamała et al. 2025 (bioRxiv 2025.11.10.687541): Welch > FFT for smooth
  spectra (ASCENT default Welch: OK); reliability DECREASES with more allowed
  peaks (fooof0 > fooof1 > fooof3) and more positive-slope outliers; alpha
  confound: low peak caps miss alpha differentially across conditions (their
  EC/EO slope effect sign-flips between fooof0 and fooof1/3). Flagged: v11
  calls P_max=6 "conservative" - defensible only vs specparam's unlimited
  default; consider one clause of nuance. ASCENT saves fit-quality flags
  (info.flags) consistent with their screening recommendation - mention as
  a strength.
- Gyurkovics et al. 2022 (J Neurosci 42:7144): aperiodic activity is
  non-stationary, ERPs must be removed before post-event fitting (v11 cites
  correctly). Their settings: peak_width_limits=[2,8], max_n_peaks=1, fixed
  mode, 2-25 Hz fitting range - ASCENT sliding [2,8] matches exactly.
- Wilson et al. 2022 SPRiNT (eLife 11:e77348): ASCENT sliding implementation
  matches SPRiNT core: 1s (fn default) Hann windows, 50% overlap, nAvg=5
  consecutive windows averaged per time bin, time bin = centre of middle
  window, per-bin specparam, fixed mode 1-40 Hz. SPRiNT sim settings: peak
  width (0.5,6), max 3 peaks; LEMON EEG settings: (1.5,6), max 3 peaks,
  min peak amplitude 0.5, fixed mode. ASCENT uses winSec=2 in GUI (0.5 Hz
  res, documented in Methods) and peak widths [2,8] - defensible.
- GAP: SPRiNT's optional-but-recommended 4th step (outlier peak pruning:
  remove peaks with <3 similar neighbours within 2.5 Hz / 6 time bins, then
  REFIT aperiodic at pruned bins) is NOT implemented in
  compute_AperiodicFit_sliding.m. Would mainly suppress transient spurious
  peaks in the sliding spectrograms. Not implemented without user decision.

## 2c. specparam MATLAB port validation vs Python reference (DONE Sept 2026)
- Scripts: validate_specparam_port.m (MATLAB side), validation_specparam/
  fit_specparam_python.py (specparam 2.0.0rc7 reference),
  compare_specparam.py (stats + figure), diagnose_diff.py (root-cause).
- Benchmark: 200 synthetic spectra (exp 0.5-2.5, off 0.5-1.5, 0-3 log-space
  Gaussians 8-14 Hz, pink multiplicative noise) + 64 real Welch PSDs from the
  sample dataset. Matched settings: fixed mode, 1-40 Hz, widths [1,12],
  max 3 peaks, h_min 0.05, threshold 2.0.
- RESULTS: real EEG exponent r=0.995 (MAE 0.052, bias +0.008, LoA +/-0.148);
  real offset r=0.962 (MAE 0.134, bias +0.083, LoA +/-0.287). Synthetic
  exponent r=0.985; offset r=0.449 driven ENTIRELY by peak-count differences
  (specparam fits more, partly spurious, peaks: 2.28/spectrum vs ASCENT 0.78
  vs 1.5 true). On peak-free spectra (n=43) the two agree to noise level
  (dOff +0.009 +/- 0.013, dExp +0.001 +/- 0.010) - aperiodic core numerically
  equivalent; deviations come from peak fitting, symmetric (both leak
  unmodeled peaks into the aperiodic fit). Disagreement concentrates in
  low-exponent tail (<0.3). Ground-truth recovery M exp r=0.946 vs P 0.938.
- GOTCHA for reruns: eeglab/plugins/roiconnect/libs/mvgc_v1.0/utils/legacy
  ships a legacy randi.m that shadows MATLAB's built-in (randi [0 3] fails
  with "Too many input arguments"); rerun scripts must rmpath it.
- Figure: validation_specparam/specparam_validation_figure.png (8 panels:
  scatter + Bland-Altman, synthetic + real). Manuscript sentence delivered
  in chat for 2.6 Diagnostics + Supplementary note.

## 3. MATLAB rerun of aperiodic figures (DONE Sept 3 2026, pending user check)
- rerun_aperiodic_figs.m ran clean; outputs in figures/subject_level/rerun/:
  aperiodic_channel_static / _timecourse, aperiodic_ica_static / _timecourse
  (.png 300dpi + .fig).
- Sept 3 (later): band-panel fix in ascent_plot_aperiodic.m (raw vs corrected
  were different quantities sharing one axis; now dual y-axes with correct
  labels). Time-course figures regenerated via rerun_aperiodic_timecourse_only.m.
- USER PROCESS NOTE: before ANY future figure-regeneration run, first collect
  and agree on ALL pending code improvements so everything is batched into a
  single rerun (avoid repeated partial reruns).
- New numbers for manuscript (FINAL, Sept 10 2026 rerun with P_max=3 + robust
  fit + SPRiNT pruning; sample dataset, 64ch, ~75s @256Hz) — user has pasted
  these into the Google Doc §3.2.2 already:
  - Channel static: r(exp,off)=0.901; exp mean 1.045 sd 0.540; offset mean
    1.137 sd 0.430; alpha raw 5.76 uV2/Hz, above-fit 4.34 (~25% aperiodic).
  - Channel sliding (70 bins): exp 1.352 sd 0.478; off 1.824 sd 0.462.
    Pruning removed 176 outlier peaks (ch) / 212 (ICA).
  - ICA static: top IC = IC14, r=0.884.
  - ICA sliding: exp 0.760 sd 0.411; off -0.549 sd 0.682.
- Known env quirk: EEGLAB's Biosig/FieldTrip plugin folders put broken
  compat stubs (pwelch.m, mean.m, strncmpi.m) on the path and break
  parpool/Signal toolbox. rerun script removes them; keep in mind for any
  future -batch runs. Parallel was disabled for these runs (small dataset).

## 4. Manuscript v11 remaining cleanups (from the earlier v8 pass, spot-check v11)
- Verify the small typos found in v8 were fixed in v11: Fig.1 caption
  ("tracjing", "paramters", "estinating", "the the", "computing..",
  "andn = 2"), "iscomputationally" (2.3), "zzz-normalized" (RCmvMFE),
  ExSent->ExSEnt, "but but" (Fig 3 caption). v11 renumbered figures (1-11)
  and most spots were rewritten; do a final grep-style proofread.
- 4.1 has a stray fragment "alpha oscillations modulate... Lombardi et al.
  (2023)" (orphan note to self, duplicated paragraph lead) - needs merging.
- "[ADD REF]" placeholder still in 4.1 (EC inhibitory dominance sentence).
- Equation (33) in Methods shows "kappa + f^chi" with 'hi' glyph corruption
  in the PDF text layer - check the source doc renders correctly.
- Discussion 4.2 mentions Kosciessa mMSE "bandpass (the original ... defaults
  to lowpass)" - double check compute_mMSE default filter_mode is what the
  group analysis actually used (ascent_group_analysis used defaults =>
  narrowband; manuscript says bandpass - consistent, but verify wording).