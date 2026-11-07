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
  - Citations to add in the next revision round, before publication (user, 2026-11-06; abstracts
    read only, read the full texts first): Ghassemkhani & Dotta 2026, Brain Topogr
    39:36, doi:10.1007/s10548-026-01176-9, cites the 2022 preprint; MDD > controls in MSE and HFD
    alongside higher beta power and a flatter aperiodic exponent, not separated: an example for
    the spectral covariates and the surrogate test. Racz et al. 2025, arXiv:2505.18117 (MMSPM):
    two aperiodic slopes, steeper below ~2 Hz; at most one clause where the fit range is
    discussed ('knee' mode exists; the 1-Hz high-pass and 1-40 Hz default fit avoid that range).
  - zNorm 1 equals z-scoring each series with the population SD, for MSE and the fuzzy measures
    (the group run used the sample SD).
  - 3.1 benchmark (sample data, 64 channels, code of this version, max |difference|): SampEn 1.1e-4
    vs Azami; FuzzEn 5.8e-5 local and 7.7e-5 global vs Azami; ExSEnt 3.1e-4 vs the authors' code; MSE
    1.1e-4 to 6.5e-4 over scales 1-5 vs EntropyHub with the same tolerance rule, as the reference in
    place of our port of Costa. No fixed-tolerance or 50-scale rerun: drop the stated cause of the
    old 5- and 50-scale differences, or run them. Serial timings 2026-11-02 (reference vs ASCENT,
    64 channels, 75 s): SampEn 70.9 vs 30.1 s, ExSEnt 1178 vs 25.2 s, MSE 5 scales 1415 vs 108 s.
    The machine was not idle (other MATLAB jobs), so rerun on a checked idle machine. FuzzEn:
    controlled run (1 thread, interleaved, 3 repeats, 4 channels): ASCENT = Azami's code in both
    modes (local 37.6 vs 38.4 s, global 44.1 vs 39.6 s); both compute all pairs with pdist, so
    serially there is no speedup to claim for FuzzEn: its gain is the parallel pool. Say that in
    3.1, with pool size and hardware.
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
- DONE 2026-11-02 (v2.2): surrogate test, 'surrogates', N for every measure except Aperiodic
  (compute_surrogates, surrogate_segment): the same call on N surrogates, saved as original, mean,
  SD, diff, z and rank p per channel and scale. Default 'gaussphase': on the segment whose ends
  match (at most 10% cut), each signal is rank-transformed to a Gaussian, then phase-randomized
  with the same phases for all signals; the original is the measure of that transformed signal.
  Chosen by simulation (SampEn, 100 signals of 3000 samples per case, 19 surrogates, % with
  p <= .05, with the exact rank test of 2026-11-03: the distance from the mean of all N+1 values;
  the first version centred on the surrogates alone and gave 6.2% under an exact null): linear
  Gaussian 5 (phase 7, IAAFT 20); linear then exp() 6 (phase 100, IAAFT 100); noisy logistic map
  100 for all three. z is now diff/(sd*sqrt(1+1/N)). IAAFT, which Kosciessa et al. used, is biased
  (mean z +1.2 on linear Gaussian signals). Plain phase rejects for a non-Gaussian amplitude distribution
  alone: on the sample dataset (eyes closed) SampEn was below its phase surrogates in 59 of 64
  channels (z -4.9) but not below its IAAFT surrogates (median z +0.3; excess kurtosis 0.31).
  In the paper: the method, these sizes and power, and what the null is.
- RUN 2026-11-03: surrogate test on the group data (EC and EO, low-pass mMSE, 30 scales, 8 channels
  F3 Fz F4 Cz P3 Pz P4 Oz, 19 'gaussphase' surrogates per recording, seed = recording index; 4.6 h).
  EC vs EO tested on the raw value, the original tested, the surrogate mean (the spectral part) and
  their difference (beyond the spectrum): sign-flip max-|t| over scales on the channel mean, and the
  cluster test. RESULT: in both conditions the data are below their surrogates at every scale
  (FWER p <= .0002); the EC < EO effect (scales 2-6, t up to -10.1) is the same in the surrogate
  mean (t -9.4 at scale 2) and absent from the difference (max |t| 1.85, FWER p >= .52): the
  low-pass mMSE effect is spectral (alpha +4.3 dB). For revision 2 (with the caveat that phase
  surrogates are stationary).
- RUN 2026-11-04/05: the same surrogate test for SampEn, FuzzEn, MSE, MFE and RCMFE (mean and SD
  coarse-graining), same 80 recordings, channels and seeds (laptop; RCMFE mean alone 13.7 h).
  RESULT: the same in all eight. EC vs EO raw and surrogate mean |t| 7.6-10; data minus
  surrogates max channel-mean |t| 0.7-1.9 (FWER p .27-.87), no cluster. RCMFE mean crosses over (EC < EO at scales 1-3 and
  20-30, EC > EO at 5-10) and the surrogates reproduce the crossover. One marginal exception:
  RCMFE SD scale 5, t +2.78, FWER p .042, no cluster (chance level for 8 tests). Arrays:
  group_surr_stats_<measure>.mat (not in the repo).
- RUN 2026-11-06: 60-Hz line noise in the group data (8 channels, every paper measure, 80 recordings,
  6.3 h). A 60-Hz peak (> 3 dB above the flanks, median channel) in 52 of 80 recordings; amplitude in
  the 8 channels median 0.029 uV, 95th percentile 0.091 uV. Notch: pop_eegfiltnew 58-62 Hz band-stop
  (-47.5 dB at 60 Hz, flat below 55 Hz). Effects in units of the between-subject SD, channel mean.
  RESULT: the notch moves no measure by more than .1 SD (TOST, 90% CI, both conditions, every
  scale) in 15 of 16 outputs; worst means .02 SD (MSEmean), worst single channel .08 SD. The EC-EO
  effect is unchanged: same t to 2 decimals, same cluster cells (MSEsd 73 -> 76, all 73 kept); the
  change in the effect reaches FWER p < .05 for MSEsd (.015), RCMFEsd (.027) and RCMFEmean (.008)
  but is at most .013 SD (every recording shifts by nearly the same tiny amount). Exception: the
  band-pass mMSE check at scale 2 (42.7-64 Hz, the only band holding 60 Hz): -.14 SD EC, -.08 EO
  (other scales <= .003). Positive control (28 recordings without a peak): adding the 95th-percentile
  amplitude shifts measures by <= .02 SD (band-pass mMSE .17); 10 times it, .2-1.6 SD for the entropy
  measures (band-pass mMSE 3.7), .06 HFD, 0 Aperiodic; the notch then removes it (<= .008 SD). Lengths:
  EC-EO -0.9 s (SD 4.8, p .24), no crop; r(length difference, EC-EO difference) up to .38
  (RCMFEmean scale 3, p .015 uncorrected over scales and measures). For revision 2: one sentence
  that line noise does not drive the results, and the band-pass caveat.
- NEXT: sampling rate (the paper lists it as a limitation). Published: SampEn falls as fs rises
  at a fixed delay (gait: McCamley et al. 2018, Raffalt et al. 2019; the only EEG study,
  Espinosa et al. 2020, upsampled 160 Hz data); a scale maps to fs/(2s), so the same scale list
  covers different bands at 256 and 512 Hz (Kosciessa et al. 2020; Zheng et al. 2023 propose
  scales in absolute time); Higuchi's k is a decimation factor, so kmax spans kmax/fs s. Nothing
  measured for FuzzEn, MFE, CMFE, RCMFE, mMSE, Katz, ExSEnt (durations in samples) or the
  multivariate measures. STUDY BUILT 2026-11-03 (branch timescales, not merged):
  assess_rate_and_length.m with simulate_eeg.m (1/f + alpha, 1/f^1 and ^2, Rossler + 1/f;
  band-passed 1-45 Hz or broadband with a white floor and an amplifier low-pass), 128-1024 Hz
  incl. 250/500/1000, four settings ('samples', 'time' = 'timebase' 256, 'resample' to 256 Hz,
  fixed number of samples), criteria fixed in advance in units of the real between-subject SD.
  New options: 'timescales' (seconds), 'timebase' (count tau, Higuchi intervals and scales at a
  reference rate), 'kmax'. RESULTS 2026-11-04 (1958 files; scratchpad sim_out/summary.csv,
  units: real between-subject SD): same settings at another rate shift every entropy and
  fractal measure by a median 1-5 SD at 512-1024 Hz (Aperiodic unaffected). Resampling to
  256 Hz is exact (0.00 SD) for band-passed data, every measure; on broadband data (amplifier
  low-pass at fs/5) SampEn, FuzzEn, ExSEnt and box-counting FD still differ by 1-4 SD, since
  tau = 1 sees the band between the two anti-alias corners. 'timebase' is close at integer
  ratios for MFE, CMFE, RCMFE, RCmvMFE, mvFuzzEn, HFD, mMSE at 512 Hz; approximate at 250/500/
  1000 Hz; no help for MSE, RCMFEsd, ExSEnt, box-counting FD. Rule: low-pass every dataset
  below the lowest rate's anti-alias corner, then resample to one rate. TO DECIDE: the default
  (warning when the rate is not 256 Hz, or an internal resampling option).
- NEXT: minimum length per measure (also a stated limitation). Published: SampEn 10^m to 30^m
  points (100-900 for m = 2; Pincus & Goldberger 1994); on 1/f noise it is biased upward below
  ~10^4 points, and pooling segments works (Grandy et al. 2016); MSE needs N/s_max >= 10^m, better
  20^m-30^m (Azami et al. 2017); fuzzy variants stay defined on short data, MFE accurate at 400
  points on HRV (Borin et al. 2021); Higuchi N >= 100, kmax(N) rule for fBm only (Wanliss &
  Wanliss 2022); ExSEnt stable at 500-3000 points on the logistic map (its authors); specparam
  segments >= 0.5 s. No study on EEG-like signals for any measure, none for filter-and-skip
  mMSE. STUDY BUILT 2026-11-03 (same script): segments of 2-60 s at 256 Hz and 2-30 s at 512 Hz
  ('time') against the process value; real EEG split-half ICC from 2 to 30 s (40 EC recordings,
  8 channels). One minimum-length rule now in all multiscale measures: StableMinBins (100 = 10^m
  for m = 2) on templates per series, or pooled over offsets where the estimator pools (CMFE on
  data under 101*S-1 samples, i.e. 11.8 s at 256 Hz, is now NaN at its coarsest scales; v2.2 gave
  values there). RESULTS 2026-11-04: real split-half ICC(A,1) lower bound >= .75 at 20 s for
  SampEn, FuzzEn, HFD, mvFuzzEn, Aperiodic and scale 2 of every multiscale measure; MSE scale 5
  needs 30 s; MFE/CMFE/RCMFE coarse scales reach ~.8 (lower ~.7) at 30 s; mMSE scales >= 10
  stay below .5 at 30 s (scale 30: .09); box-counting FD .3-.4 at every length and biased by
  length itself (-2.75 SD at 60 s vs 180 s); ExSEnt HA/HDA under .7. Simulation RMSE < .5 SD:
  MSE/mMSE need ~900-1500 samples per coarse series (30^m), fuzzy variants far fewer. The bias
  equivalence test (16 signals) is underpowered at coarse scales, so "never" there means not
  shown, not biased. TO DO: per-measure length warning (20 s; MSE/mMSE N/s < 900); a note
  that box-counting FD needs equal lengths; RCmvMFE allows a 2-GB pdist temporary per worker,
  which ran out of memory once at 8 workers.
- NEXT (rigor of the real-data part): the group recordings are ~65 s, so halves stop at 30 s
  (measures still below .75 there are censored, not settled), adjacent halves of one minute
  mix estimator noise with drift and leave out between-session variance, and the yardstick SD
  comes from 60-s values. Planned: (1) Wang et al. 2022 test-retest dataset (OpenNeuro
  ds004148, 60 adults, 5 min EO then 5 min EC per session, sessions 90 min and 1 month apart,
  Brain Products 63/64 ch, 500 Hz; 33 GB for all five states, EC/EO only is a subset): split
  halves up to 150 s, within-session non-adjacent segments, same-day and 1-month ICC, yardstick
  from 5-min values. (2) MPI LEMON raw EEG (Babayan et al. 2019; 62 ch, 2500 Hz, hardware band
  0.015 Hz-1 kHz, 16 alternating 60-s blocks, 8 EC; ~317 MB per subject, ~220 subjects): the
  rate study on real broadband EEG (resample to 1024/512/500/256/250/128) and reliability across
  non-adjacent EC blocks; a subset of ~50 subjects. (3) Simulations: references of 300 s and
  segments up to 120-300 s for the coarse scales. Downloads not started.
- RUNNING 2026-11-05 (Mac, from 11:14): ds004148 EC preprocessing (180 recordings), then analysis
  A trimmed (no CMFE, no FuzzEnGlobal, no 295-s value in session 2, REST arm deferred), 6 workers,
  runners in .claude/mac-runners/, results ~/ascent/results/A (one file per recording x measure x L).
  Statistics written and checked on simulations; they run on their own when A ends (expected
  2026-11-06 evening; ~170 worker-h in all). QUEUED on the Mac, in this order (user, 2026-11-06:
  all of these on the Mac): (1) read A: minimum length, ICCs, drowsiness sensitivity; (2) the
  paper's 3.1 timings on the idle Mac; (3) LEMON (B), the 2500-Hz raw (~25 GB; the LEMON copy in
  the challenge folder is the 250-Hz preprocessed release): the Mac has 85 GB free since the user
  set duplicate Proton folders and DATA/Ascent/raw to "Online Only" (cloud copies kept), guard
  12 GB; (4) surrogates that keep non-stationarity
  (PLANNED below), built and validated on simulations, then the group data.
- DONE 2026-11-05 (branch timescales, commit ab93b01): exact speedups, outputs bit-identical to the
  previous code. functions/sampen_exact.m: SampEn match counts with templates sorted by their first
  coordinate (a lag sweep over sorted neighbours within r), O(N) memory; all pairs with pdist for
  short series (<= 5e5 pairs) and when > .4/maxNumCompThreads of pairs are within r on the first
  coordinate. compute_SampEn, compute_MSE and compute_mMSE use it (MSE/mMSE without the
  compute_SampEn front end). Blocked fuzzy sum: across-block tiles from pdist2. Validated: new
  tests/test_exact_counts.m (320 settings x 3 routes vs pdist counts, short and 25000-sample series,
  pdist2 tile = max-abs loop); old vs new on real EEG (sample data on the laptop, ds004148 sub-01 on
  the Mac): every output isequaln for SampEn, MSE (mean, sd, zNorm 1-4), mMSE (6 variants incl.
  time-resolved), ExSEnt, FuzzEn/MFE/mvFuzzEn/RCmvMFE on the blocked path. Mac, one thread, 1 ch:
  SampEn 3.3x (60 s) and 5.5x (290 s), MSE 4.2x/5.1x, mMSE 2.2x/3.6x, ExSEnt 3.4x, FuzzEn 120 s
  1.03x, RCmvMFE 120 s scale 1 1.28x. The paper's 3.1 timings (SampEn, MSE vs reference) need
  rerunning with this code. Not done: sharing CMFE/RCMFE/MFE sums; a MEX.
- SPEED (audit 2026-11-04, scratch audit/). All entropy measures share one
  path: pdist (Chebyshev) into a <= 2 GB vector, then sum(exp(-d.^2/r)); above 2 GB a blocked
  loop. Bit-identical wins: pdist2 for the blocked across-block tile (1.49x at dim 17, so ~1.25x
  for RCmvMFE); sorted-template SampEn (prototype 3.6x, 60/60 cases identical) used directly by
  MSE/mMSE (mMSE loses ~37% to compute_SampEn call overhead); share the (pm, pm1) sums of
  CMFE/RCMFE/MFE when settings match; parfor over row blocks for RCmvMFE scales 1-5 (wall time).
  Round-off-level: global FuzzEn by lag blocks (1.49x, |dFuzzEn| <= 4.7e-15). Only a MEX
  (distance + kernel fused, SIMD exp) gives 10-30x; needs compilers for Windows/macOS arm64 and a
  test tolerance ~1e-12. Not reachable: GEMM (max-norm), truncating the fuzzy kernel (too wide).
- PLANNED: surrogates that keep non-stationarity, to split each measure into a spectral, an
  amplitude-modulation and a nonlinear part. Full texts read 2026-11-04 (Lancaster 2018,
  Breakspear 2003, Xiong 2017, Keylock 2019): no published benchmark on linear non-stationary
  data, so validate ours on such simulations; one randomization shared by all channels; keep
  slow trends or high-pass data and surrogates alike; pinned-wavelet results over a range of
  pinned energy fractions, not one. Phase surrogates are stationary, and for SampEn on
  linear AR processes amplitude modulation alone gives > 80% rejections with FT surrogates
  (Faes et al. 2009, doi:10.1109/TBME.2008.2009358). First choice: band-envelope surrogates
  (filter bank, ~4 bands per octave; envelopes smoothed over max(0.5 s, 3 cycles); carriers
  phase-randomized with the same phases for all channels, then remodulated; IAAFT-style
  iterations for the global spectrum, envelopes and ranks), built on Chavez & Cazelles 2019
  (doi:10.1038/s41598-019-43571-2) and Keylock's IAAWT; it reduces to the current surrogate
  for infinite smoothing. Checks: pinned wavelet IAAFT (Keylock 2007, doi:10.1016/
  j.physd.2006.10.012) at several pinning fractions; time-varying AR (Faes 2009) on a few
  channels. Validate before real data: size on linear non-stationary processes (bursting
  alpha, frequency drift, spectral switching, slow drift, static distortion, 64-channel
  mixtures), power on nonlinear ones (phase-locked harmonics, phase-amplitude coupling,
  asymmetric waveforms, Stuart-Landau, neural mass incl. bistable alpha, chaos by SNR), a
  sweep of the smoothing and band settings fixed before real data. Still to read in full:
  Lancaster et al. 2018 (doi:10.1016/j.physrep.2018.06.001), Breakspear et al. 2003
  (doi:10.1016/S0167-2789(03)00136-2).
- ascent_group_analysis.m (branch timescales): optional surrogate section added (run_surrogates,
  off by default; results in data_path/ascent_surrogates_<tag>/). Seen, not fixed: the NaN
  check after each compute loop tests MFE(:,2:end,:) only and its disp comes after error; the
  MSE section has no 'any(mask)' guard before pull_clusters_by_sign/plot_results; output names
  say 'tfce'/'uncorrected' where mcc_type 2 is cluster correction; the correlation section runs
  2000 permutations for t-values only; 'Prosessing' in progress titles; pull_clusters_by_sign
  fails on single-scale maps (no Start/End).
- CHECK on the desktop (figures hang headless): multiscale plots with 'timescales'/'timebase'
  (axes in ms) and mMSE plots (now numeric scale factors).
- NOTE: specparam with the default settings fits about 3 Gaussians to pure 1/f noise and raises
  the exponent by about .01 (MaxPeaks 0: exact). Defaults unchanged; mention if relevant.
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
  - EEGLAB upload: version 2.2 only, not 2.1 (user, 2026-11-06), from the zip rebuilt with the
    revision. Not on the EEGLAB plugin list as of 2026-11-06. The form failed on 2026-10-22 (a
    server-side PHP issue according to Arnaud Delorme): name ASCENT, docs
    https://github.com/amisepa/Ascent.
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
- NEXT: wait for the editor's decision; EEGLAB plugin form upload with version 2.2.
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
  6. Upload to the EEGLAB extension manager (plugin_uploader form): version 2.2, not 2.1.
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