# ASCENT: TODO

## Status 2026-10-30 (supersedes everything below where they conflict)
Revision 1 is with the editor. The working tree is v2.2, not committed, and differs from the
ASCENT2.1.zip release.
- DONE in code: mMSE low-pass filt-skip default ('FilterType' fir or butter, 'StartPoints'
  pooled or average), band-pass kept as the check; spectral covariates saved next to every
  measure except Aperiodic (EEG.ascent.<measure>.covariates); scale frequency labels
  ('scalefreq'); tau reaches SampEn, n the fuzzy measures, Kernel MFE, fuzzy_mode FuzzEn, MFE and
  CMFE, IncludeScale1 MFE and CMFE; zNorm 1 scales r by std^n for the exponential fuzzy kernel
  (std for the Gaussian one), so it equals z-scoring each series with the population SD; RCmvMFE
  scale 1; scale-2 drop in the SD/variance plots; slidPeakWidthLimits; Parallel Computing
  Toolbox optional for real; eegh command complete and valid; 'save','auto' compares the data;
  warnings for epochs and boundary events; channel selection keeps the user's order and errors
  on a misspelled label; group script and benchmark script edited, outputs of the group script
  in their own files (outputs_tag), so the cached outputs behind the figures are not touched.
- CHECKED on one subject (EC, 6 channels) against the cached outputs behind the figures, with the
  group script's settings: SampEn, FuzzEn, HigFracDim, ExSEnt, exponent and RCMFE (SD and mean)
  match exactly. MSE and MFE with 'zNorm', 1 and mMSE band-pass/average differ by at most .012
  (median < .001): the group run z-scored each series with the sample SD inside
  compute_SampEn/compute_FuzzEn (each starting point for mMSE), and doing that reproduces them to
  1e-15 (Oz); zNorm_scale uses the population SD.
- GROUP RESULT 2026-10-30 (EC vs EO, 40 subjects, 2000 permutations, cluster correction, as in
  make_revision_figures): low-pass mMSE has one negative cluster (EC < EO), scales 2-12, all 64
  channels, peak F2 at scale 2, t = -10.17, g = -1.58. The cached band-pass mMSE gives no cluster,
  as in the paper. Averaged over channels, EC-EO at scales 2-4 correlates with the EC-EO alpha
  difference across subjects (r = -.50 to -.55), so the effect is at least partly spectral; only the
  surrogate test can say more. Arrays: biosemi_data/ascent_mmse_lowpass_v22.mat.
- NOT CHECKED: the group script end to end (its compute part ran on 1 subject, 4 channels, 12 s,
  both conditions and both coarse-grainings; statistics and figures need the desktop); benchmark
  sections after MSE; the new tick labels on a real figure; the case without the Parallel Computing
  Toolbox (it is installed here).
- NEXT (code): the figures need the desktop: run the group script's stats and
  make_revision_figures there and look at the new tick labels; rerun
  assess_speed_and_outputs.m on an idle machine for serial timings; rebuild the release zip
  (version 2.2) with the revision; EEGLAB upload.
- EO ICA used 'maxiter', 500 as EC (confirmed 2026-11-02); the group script now says so.
- FOR REVISION 2 (text, with the group result):
  - mMSE: low-pass default and pooled match counts (Wu et al. 2014, a new reference) in 2.2.5, the group Methods,
    4.2 and the reply on default bandpass filtering; band-pass as Kosciessa et al.'s check. Remove
    "not implemented". State the filter: FIR, closer to an ideal low-pass than the 6th-order
    Butterworth. Say the v2.1 results were band-pass.
  - Band-pass caveat: the top of each band (12% at scale 2, 38% at scale 30) folds back on
    decimation, so "a defined frequency band" is approximate.
  - Spectral covariates and scale frequency labels: Methods, and the limitation on power control.
  - Surrogate test (see below): wanted in the paper.
  - zNorm 1 equals z-scoring each series with the population SD, for MSE and the fuzzy measures
    (the group run used the sample SD).
  - 3.1 benchmark (sample data, 64 channels, code of this version, max |difference|): SampEn 1.1e-4
    vs Azami; FuzzEn 5.8e-5 local and 7.7e-5 global vs Azami; ExSEnt 3.1e-4 vs the authors' code; MSE
    1.1e-4 to 6.5e-4 over scales 1-5 vs EntropyHub with the same tolerance rule, as the reference in
    place of our port of Costa. No fixed-tolerance or 50-scale rerun: drop the stated cause of the
    old 5- and 50-scale differences, or run them. Serial timings on an idle machine still to come:
    on the loaded machine FuzzEn global took 334 s vs 270 s for the reference, so check
    "reduced for all benchmarked measures".
  - The 69 GiB crash claim belongs to RCmvMFE, not mvFuzzEn (text and reply).
  - "exceed available RAM" -> "exceed a set size (2 GB by default)".
  - 2.7 and Data Availability: name the commit once the group script has run end to end; say the
    published figures come from the cached outputs, and that the current code gives slightly
    different MSE and MFE (up to .012 on the one subject compared; sample vs population SD).
  - "ASCENT does not introduce new entropy definitions" sits badly with the mMSE differences and
    the mvFuzzEn averaging.
  - Reply to R1, now "Not applicable.": the title and scope were revised (abstract, 1.3:
    computation and export, statistics done externally), and the authors will gladly address the
    details given to the editor.
- PLANNED FOR REVISION 2: surrogate test, optional for every entropy/fractal measure: IAAFT or
  Fourier phase randomization, n surrogates; output the measure, the surrogate mean and original
  minus surrogate per channel and scale (Kosciessa et al. 2020; Courtiol et al. 2016). In the
  plugin first, with a check on the sample data; group level if time allows.
- PLANNED after the paper (the v2.2 items below move to the next version): Kosciessa's high-pass
  check; MRMSSD (Courtiol et al. 2016); alpha band-stop check; simulation mode varying one
  parameter at a time (Donoghue et al. 2024 design; fuzzy and multiscale entropy were not in it);
  eBOSC-style rhythm detection; corrected outputs labeled by rule and exports for the oneoverf
  work, once that paper is out; per-epoch computation with pooled counts instead of joining epochs;
  skip boundary segments; cap the band-pass edge at the post-skip Nyquist, if band-pass stays.
- KEPT AS IS (decided): compute_SampEn counts m and m+1 matches over different template sets (as
  the published values); same-template counts would change them slightly, and matter only for
  very short windows or tau > 1.
- DONE 2026-11-02: one fill-and-z-score step for MSE, MFE, CMFE, RCMFE and RCmvMFE (dead channels
  give NaN, NaN gaps are filled, so the zNorm fallback no longer jumps); RCmvMFE leaves out dead
  channels; ExSEnt warns that it has no tau; knee sign as in the fit in the band power; datasets
  without electrode locations (ECG, RR) are computed without figures; no IC maps for the
  multivariate measures; a warning when the sampling rate leaves few bins for the aperiodic fit;
  the tutorial saves the dipole fit to a copy; single-scale console tables removed from the group
  script; RCmvMFE help.

## Status 2026-10-22 (supersedes the 10-14 "NEXT" below)
Plugin v2.1 is on GitHub (release v2.1, ASCENT2.1.zip). Manuscript v16 and letter v6 are built
in manuscript/ (no PDFs), with all figures current; next is the final read and submission.
- NEXT (paper, to submit the revision):
  1. DONE 2026-10-22 14:48: aperiodic figures regenerated, v16 and letter v6 rebuilt
     (manuscript/), README Figs 4-5 updated.
  2. Final read of v16 red + letter v6, then submit (see 3).
  3. Submit v16 red + letter v6. IRB protocol code: none available,
     v16 cites the IONS IRB and the dataset paper.
- NEXT (plugin):
  - Upload ASCENT2.1.zip through the EEGLAB plugin form once it works again (down on
    2026-10-22, a server-side PHP issue according to Arnaud Delorme): name ASCENT, version 2.1,
    docs https://github.com/amisepa/Ascent.
  - DONE 2026-10-22: old copies deleted; eeglab/plugins/ASCENT is the clean install from
    ASCENT2.1.zip (52 files, no .git).
  - Optional: ask GitHub Support to purge unreferenced objects after the history rewrite.
  - Optional: GitHub Actions smoke test (matlab-actions + EEGLAB, README commands on a 20-s
    excerpt of the sample).
- DONE 2026-10-23: IC outputs with ICLabel classes, AAL labels (ascent_ic_info), dipole fit
  quality and coregistration check; Fig 3, v16 and README regenerated; reference 46 (AAL)
  added to the manuscript. Release zip rebuilt. Dipole figure for every IC measure except the
  multivariate ones (ascent_plot_dipoles); guards for non-MNI fits, unreadable head models,
  DIPFIT models that no longer match the ICs, and plotting errors (results are kept).
- DONE 2026-10-23: Revision 1 submitted to Entropy (v16, letter v6, Figure S1, cover letter).
- NEXT: wait for the editor's decision; EEGLAB plugin form upload when their server is fixed.
- PLANNED v2.2: brain-only option (ICLabel) for the IC figures and cluster detection.
- PLANNED v2.2: source-space entropy through ROIconnect (LCMV/eLORETA ROI time series), the
  route for any claim about where entropy differs.
- PLANNED v2.2 (after the paper; not in the paper): eLORETA option for IC maps (template BEM
  leadfield, ft_inverse_eloreta as in the MMAGIC-validation code), shown next to the dipoles;
  report a lobe-level label, keep the RV < 15% dipole view as the default.
- PLANNED v2.2 (after the paper; not in the paper): attractor measure
  - 3-D time-delay embedding of 20-s windows, x(t), x(t+tau), x(t+2*tau), tau about 4-5 ms
    (1 sample at 250 Hz), data 1-50 Hz; minimum-volume enclosing ellipsoid (Khachiyan
    algorithm); radius ratio = shortest / longest axis (1 = sphere, awake; lower = flatter,
    loss of consciousness). Sliding windows (5-s steps) for the time course. Optional
    correlation dimension (Grassberger-Procaccia).
  - Outlier-robust option (trim extreme points or principal axes), since the enclosing
    ellipsoid is stretched by single artifacts. Note in the help: at short delays the ratio
    tracks fast vs slow activity (rho 0.71-0.91 with gamma power in Eagleman et al. 2019).
  - Figures: topography of the ratio, 3-D point clouds with the fitted ellipsoid for chosen
    channels or ICs, ratio over time. GUI entry, README tutorial step, references.
  - References: Eagleman, Chander, Reynolds, Ouellette & MacIver (2019) PLoS ONE 14,
    e0223921; Eagleman, Drover, Drover, Ouellette & MacIver (2018) Front Hum Neurosci 12,
    173; Eagleman, Vaughn, et al. (2018) Front Neurosci 12, 645; MacIver & Bland (2014)
    Front Syst Neurosci 8, 203; Khachiyan (1980); Takens (1981).

## Status 2026-10-14 (supersedes the 09-25 "NEXT": do NOT submit v15)
Review of v15 done; v16 is one build away.
- FIXED in the working tree, NOT committed: compute_psd.m (detrend across channels),
  compute_AperiodicFit_sliding.m (pruning column), ascent_plot_multiscale.m (Fig 3 dipole
  panel shows all dipolar ICs), ascent_compute_gui.m (headers), make_revision_figures.m
  ('rcmfe' section, single-scale tests with scalp adjacency), Fig S1.
- DONE (desktop): cached group outputs updated (backups *_before_detrend_fix.mat); Figs 2,
  4-7, 10-11 regenerated; regen CSVs match the headless checks.
- NEXT:
  1. DONE 2026-10-14 11:09: Figs 2-3 regenerated (exact cluster outlines, 13-pt fonts,
     Fig 2 panels stacked, MRI without the stray slab) and checked composed.
  2. DONE 2026-10-22: Fig 1 screenshots in manuscript/figures/fig1_v16/ (same files as
     img/gui_main.png, gui_multiscale.png, gui_aperiodic.png); caption A updated.
  3. IRB: statements fixed and Informed Consent added (v16_edits_stage1.py). No protocol
     code in [46] or on OSF; parked 2026-10-14: v16 goes without one.
  4. `bash manuscript/tools/build_v16.sh <scratch> <fig1 dir>`; PDFs with
     `pwsh -File manuscript/tools/render_pdf.ps1 <docx> ...`; check every figure page and
     the numbers against manuscript/figures/regen/*.csv; then submit v16 red + letter v6.
  5. DONE 2026-10-22: code committed and pushed (history rewritten and slimmed to 57 MB);
     GitHub release v2.1 with ASCENT2.1.zip. v16 manuscript + letter v6 built (no PDFs).
  6. Upload ASCENT2.1.zip to the EEGLAB extension manager (plugin_uploader form).
- Decided: no disclosure statement added; Fig 2 uses the default 75th-percentile threshold.
- OPEN (code, not fixed): ascent_compute.m:505-509 drops scale 2 in SD/var plots;
  (fixed 2026-10-14: 'save' option ignored); tau never passed to SampEn, n never
  to MFE/RCmvMFE; ascent_group_analysis.m no longer reproduces the paper (needs
  'zNorm' 1 for MSE, 2 for MFE, IncludeScale1 for RCMFE-mean; EO preprocessing
  block differs from EC); assess_speed_and_outputs.m cannot be rerun as shipped;
  single-scale tests in ascent_group_analysis.m:463-468 pass [chan x subj], which
  the statistics code clusters along the channel list (pass [chan x 1 x subj] as
  make_revision_figures.m now does); compute_RCmvMFE scale 1 is always 0
  (compute_RCmvMFE.m:150, by reading, not run).
  (fixed 2026-10-22: ascent_sample_data.set no longer carries stale EEG.ascent outputs.)

## Status 2026-09-25 (supersedes older notes below where they conflict)
- DONE: knee-mode aperiodic correction (full Lorentzian removed) - item 2.
- DONE: compute_AperiodicFit now ports specparam 2.0 exactly (robust initial
  fit, peak search in Hz with edge/overlap drops, plain final refit). Validation
  vs Python specparam (validation_specparam/matlab_results_v2.mat): r > 0.999 for
  exponent and offset, real and synthetic; identical peak counts. Old 2c numbers
  are obsolete.
- DONE: ascent_group_analysis.m ExSEnt HDA now gets its own compute_mcc (was
  reusing the HA mask).
- DONE 2026-09-25: v14 docx + PDF (manuscript/) with the regenerated Figs 4–7,
  10, 11, text checked against them, Fig 6 comment removed, Figs 1/6 no longer
  overflow the page.
- DONE 2026-09-25 evening: v15 (manuscript/): clean, changes since v6 SUBMITTED
  in red, edits since v14 (tracked), PDFs; letter v5 clean + tracked. Figures
  4–11 regenerated, phrasing review, Discussion restructured, references checked
  against Crossref. Rebuild: `bash manuscript/tools/build_v15.sh <scratch>`.
- NEXT: import v15 clean into Google Docs as the new master; submit the
  red version + letter v5; delete the old blue-tracked Google Doc copy.
- OPEN (eeg_robust_statistics repo): pull_clusters' 2D path labels clusters on
  the unsigned mask, so adjacent positive/negative regions merge; ASCENT scripts
  now use a local pull_clusters_by_sign. Consider a sign-split option there.
- OPEN: knee mode still differs from specparam (offset MAE ~0.5) - starting
  guesses/bounds of the knee fit.
- OPEN: ascent_compute passes no tau to compute_SampEn (tau always 1);
  'slidPeakWidthLimits' name-value never matches (lowercased key);
  sd/var coarse-graining drops scale 2 twice in ascent_compute:503-506 plots.
- OPEN: cluster overlay default (item 1).

# Older notes (Sept 2026 revision session)

## 1. Cluster overlay on the main multiscale heatmap (ascent_plot_multiscale.m)
- Unsure about the recent automatic cluster detection that decides what
  gets plotted in the scalp topo + curve plots (manuscript Fig. 2).
- Last time it required manual threshold adjustment to look right; there
  should probably be only ONE cluster shown, not up to 3.
- IDEAS (not implemented):
  a) Show only the top-1 cluster by default.
  b) Add an interactive SLIDER on the figure so that one can move the
     percentile threshold live; the topo and curve panels update dynamically
     as the slider moves. (uicontrol 'style','slider' on the figure, callback
     re-runs find_entropy_clusters + redraws the topo/curve subplots.)
- Related files: ascent_plot_multiscale.m (find_entropy_clusters,
  clusterThresh opt, default 0.75), ascent_plot.m ('ClusterThresh' pass-through).

## 2. Knee mode + aperiodic subtraction (done, see the 09-25 status)
- Current code disables aperiodic correction in knee mode (ascent_compute.m,
  compute_AperiodicFit_sliding.m, GUI callback greys out the checkbox).
- CHECK (Sept 2026): upstream specparam/FOOOF has NO such restriction.
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
- NOT YET RERUN - batch with other pending changes.

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
  peaks in the sliding spectrograms. Not implemented.

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
  scatter + Bland-Altman, synthetic + real). Manuscript sentence for
  2.6 Diagnostics + Supplementary note.

## 3. MATLAB rerun of aperiodic figures (DONE Sept 3 2026)
- rerun_aperiodic_figs.m ran clean; outputs in figures/subject_level/rerun/:
  aperiodic_channel_static / _timecourse, aperiodic_ica_static / _timecourse
  (.png 300dpi + .fig).
- Sept 3 (later): band-panel fix in ascent_plot_aperiodic.m (raw vs corrected
  were different quantities sharing one axis; now dual y-axes with correct
  labels). Time-course figures regenerated via rerun_aperiodic_timecourse_only.m.
- PROCESS NOTE: before ANY future figure-regeneration run, first collect
  and agree on ALL pending code improvements so everything is batched into a
  single rerun (avoid repeated partial reruns).
- New numbers for manuscript (FINAL, Sept 10 2026 rerun with P_max=3 + robust
  fit + SPRiNT pruning; sample dataset, 64ch, ~75s @256Hz) - pasted
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