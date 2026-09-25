function [exponent_t, offset_t, times, freqs, psd_t, psd_corrected_t, fit_error_t, peaks_t] = ...
        compute_AperiodicFit_sliding(data, fs, varargin)
% compute_AperiodicFit_sliding  Time-resolved aperiodic estimation following
%   the SPRiNT architecture (Wilson et al., 2022, eLife).
%
%   Short overlapping FFT windows are computed across the full signal.
%   At each time bin, the power spectra of nAvg consecutive windows are
%   averaged to produce a smoothed PSD estimate, which is then
%   parameterized with specparam.  This separates frequency resolution
%   (controlled by winSec) from spectral SNR (controlled by nAvg),
%   yielding finer temporal resolution than a single long Welch window.
%
% INPUTS
%   data   - [nChan x nSamples], µV
%   fs     - sampling rate (Hz)
%
% OPTIONAL NAME-VALUE
%   'winSec'          - short FFT window length in s (default: 1)
%                       Determines frequency resolution (1/winSec Hz).
%   'winOverlap'      - fractional overlap between consecutive windows 0-1
%                       (default: 0.75). Controls step size and thus the
%                       temporal resolution of output time bins.
%   'nAvg'            - number of consecutive windows to average per time bin
%                       (default: 5). Trades temporal resolution for spectral
%                       SNR: larger nAvg = smoother PSDs but coarser time axis.
%   'freqRange'       - [fMin fMax] Hz for both PSD and fit (default: [1 40])
%   'aperiodicMode'   - 'fixed' | 'knee' (default: 'fixed')
%   'maxPeaks'        - max Gaussian peaks per fit (default: 6)
%   'minPeakHeight'   - minimum peak height above aperiodic, log10 units (default: 0.05)
%   'peakThreshold'   - peak detection threshold in SD of flattened spectrum (default: 2.0)
%   'peakWidthLimits' - [min max] peak width in Hz (default: [2 8])
%   'correctAperiodic'- subtract aperiodic model from PSD (default: true;
%                       in knee mode the full Lorentzian model is removed)
%   'PrunePeaks'      - SPRiNT-style outlier peak pruning (default: true).
%                       After fitting, peaks with fewer than PruneMinPeaks
%                       similar peaks (within PruneFreqTol Hz) in nearby time
%                       bins (PruneTimeBins) are removed, and the aperiodic
%                       model is refit at the affected bins on the
%                       peak-removed spectrum. Reduces transient spurious
%                       peaks in the time-resolved spectrograms.
%   'PruneMinPeaks'   - minimum similar neighbouring peaks to keep a peak
%                       (default: 3; SPRiNT recommendation)
%   'PruneFreqTol'    - centre-frequency tolerance for "similar" peaks, Hz
%                       (default: 2.5; SPRiNT recommendation)
%   'PruneTimeBins'   - neighbouring time-bin window for the neighbour count
%                       (default: 6; SPRiNT recommendation)
%   'Parallel'        - parfor over channels inside each time bin (default: true)
%   'Progress'        - print progress (default: true)
%
% OUTPUTS
%   exponent_t      - [nChan x nTimes]  aperiodic exponent per time bin
%   offset_t        - [nChan x nTimes]  aperiodic offset per time bin (log10 units)
%   times           - [1 x nTimes]      centre time of each averaged window group (s)
%   freqs           - [1 x nFreqs]      frequency vector (Hz)
%   psd_t           - [nChan x nFreqs x nTimes]  averaged PSD per time bin (µV²/Hz)
%   psd_corrected_t - [nChan x nFreqs x nTimes]  aperiodic-corrected PSD;
%                     empty ([]) if correctAperiodic = false
%   fit_error_t     - [nChan x nTimes]  fit MAE per time bin
%   peaks_t         - {nChan x nTimes}  peak params [CF power BW] per bin
%                     (outlier peaks removed when PrunePeaks = true)
%
% NOTES
%   Architecture follows SPRiNT (Wilson et al., 2022):
%     Step 1 - Compute Hann-windowed FFT power for every short window.
%     Step 2 - At each time bin, average the power of nAvg consecutive
%              windows to obtain a smoothed single-sided PSD.
%     Step 3 - Fit the specparam aperiodic model to that averaged PSD.
%     Step 4 - (optional) Subtract the aperiodic model.
%     Step 5 - (optional, SPRiNT-recommended) Prune outlier peaks that do
%              not recur across neighbouring time bins, and refit the
%              aperiodic model at the pruned bins.
%
%   Temporal resolution = winSec * (1 - winOverlap) seconds per step.
%   Frequency resolution = 1/winSec Hz.
%   Total output time bins = nWins - nAvg + 1, where
%     nWins = floor((nSamples - winSamples) / stepSamples) + 1.
%
%   Aperiodic correction removes the fitted aperiodic model in log space:
%   fixed: log10(P_corr) = log10(P) - (offset - exponent * log10(f))
%   knee:  log10(P_corr) = log10(P) - (offset - log10(knee + f^exponent))
%
% REFERENCES
%   Wilson LE, da Silva Castanheira J, Baillet S (2022). Time-resolved
%     parameterization of aperiodic and periodic brain activity.
%     eLife, 11, e77348. https://doi.org/10.7554/eLife.77348
%   Donoghue T, et al. (2020). Parameterizing neural power spectra into
%     periodic and aperiodic components. Nature Neuroscience, 23, 1655-1665.
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Ascent Plugin - Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

%% --- Parse inputs ---
p = inputParser;
addRequired(p,  'data',            @(x) isnumeric(x) && ismatrix(x));
addRequired(p,  'fs',              @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'winSec',          1,       @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'winOverlap',      0.75,     @(x) isnumeric(x) && x >= 0 && x < 1);
addParameter(p, 'nAvg',            5,       @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'freqRange',       [1 40],  @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'aperiodicMode',   'fixed', @(s) ischar(s) || isstring(s));
addParameter(p, 'maxPeaks',        6,       @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'minPeakHeight',   0.05,    @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'peakThreshold',   2.0,     @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'peakWidthLimits', [2 8],  @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'correctAperiodic',true,    @islogical);
addParameter(p, 'PrunePeaks',      true,    @islogical);
addParameter(p, 'PruneMinPeaks',   3,       @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'PruneFreqTol',    2.5,     @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'PruneTimeBins',   6,       @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'Parallel',        true,    @islogical);
addParameter(p, 'Progress',        true,    @islogical);
parse(p, data, fs, varargin{:});
opt = p.Results;

opt.aperiodicMode = lower(char(opt.aperiodicMode));
if ~ismember(opt.aperiodicMode, {'fixed','knee'})
    error('compute_AperiodicFit_sliding:badMode', ...
        'aperiodicMode must be ''fixed'' or ''knee''.');
end

%% --- Window setup ---
[nChan, nSamp] = size(data);
winSamp  = round(opt.winSec * fs);
stepSamp = round(opt.winSec * (1 - opt.winOverlap) * fs);

starts = 1 : stepSamp : (nSamp - winSamp + 1);
nWins  = numel(starts);

if nWins < opt.nAvg
    error('compute_AperiodicFit_sliding:tooFewWindows', ...
        ['Only %d windows available but nAvg=%d. ' ...
         'Reduce winSec, winOverlap, or nAvg.'], nWins, opt.nAvg);
end

nTimes = nWins - opt.nAvg + 1;   % one output bin per valid averaged group

freqRes = 1 / opt.winSec;
minPeakWidth = 2 * freqRes;
if opt.peakWidthLimits(1) < minPeakWidth
    fprintf('Note: peakWidthLimits lower bound raised from %.2f to %.2f Hz (2x freq resolution for winSec=%.1fs)\n', ...
        opt.peakWidthLimits(1), minPeakWidth, opt.winSec);
    opt.peakWidthLimits(1) = minPeakWidth;
end

%% --- Frequency vector from FFT resolution ---
% Single-sided frequencies for a winSamp-point FFT at rate fs
freqs_fft = (0 : winSamp-1) * (fs / winSamp);
freq_mask  = freqs_fft >= opt.freqRange(1) & freqs_fft <= opt.freqRange(2);
freqs      = freqs_fft(freq_mask);
nFreqs     = numel(freqs);

if nFreqs == 0
    error('compute_AperiodicFit_sliding:noFreqs', ...
        ['No FFT frequencies fall within freqRange [%g %g] Hz. ' ...
         'Increase winSec or widen freqRange.'], opt.freqRange(1), opt.freqRange(2));
end

%% --- Centre times ---
% Time bin i corresponds to averaged windows i : i+nAvg-1.
% Centre = centre of the middle window of that group.
mid_offset = floor(opt.nAvg / 2) * stepSamp;   % samples from first window start
times = ((starts(1:nTimes) - 1) + winSamp/2 + mid_offset) / fs;

%% --- Step 1: FFT power for every short window ---
% Hann taper and single-sided µV²/Hz normalization.
hann_win    = hann(winSamp)';               % 1 x winSamp row vector
norm_factor = 2 / (fs * sum(hann_win.^2)); % factor for µV²/Hz, single-sided

if opt.Progress
    fprintf(['SPRiNT-style sliding aperiodic fit\n' ...
             '  winSec=%.1fs  overlap=%.0f%%  step=%.2fs  nAvg=%d  ->  %d time bins\n'], ...
        opt.winSec, opt.winOverlap*100, stepSamp/fs, opt.nAvg, nTimes);
    fprintf('  Step 1: computing FFT power for %d windows...\n', nWins);
end

psd_all = nan(nChan, nFreqs, nWins);
for iWin = 1:nWins
    seg     = data(:, starts(iWin) : starts(iWin)+winSamp-1);  % nChan x winSamp
    seg_win = seg .* hann_win;                                  % apply taper
    fft_pwr = abs(fft(seg_win, [], 2)).^2;                      % nChan x winSamp
    psd_all(:, :, iWin) = norm_factor * fft_pwr(:, freq_mask);
end

%% --- Pre-allocate outputs ---
exponent_t = nan(nChan, nTimes);
offset_t   = nan(nChan, nTimes);
psd_t      = nan(nChan, nFreqs, nTimes);
if opt.correctAperiodic
    psd_corrected_t = nan(nChan, nFreqs, nTimes);
else
    psd_corrected_t = [];
end

%% --- Steps 2-4: average windows, fit specparam, optional correction ---
if opt.Progress
    fprintf('  Step 2: fitting aperiodic model to %d time bins...\n', nTimes);
end

fit_error_t = nan(nChan, nTimes);
peaks_t     = cell(nChan, nTimes);
knee_t      = cell(nChan, nTimes);   % full aperiodic param vector per bin (knee = param 2)

for iTime = 1:nTimes

    % Average power across nAvg consecutive windows
    psd_avg = mean(psd_all(:, :, iTime : iTime+opt.nAvg-1), 3);   % nChan x nFreqs
    psd_t(:, :, iTime) = psd_avg;

    % Fit aperiodic model to smoothed PSD
    [exp_win, off_win, info_win] = compute_AperiodicFit(freqs, psd_avg, ...
        'FreqRange',       opt.freqRange,       ...
        'AperiodicMode',   opt.aperiodicMode,   ...
        'MaxPeaks',        opt.maxPeaks,        ...
        'MinPeakHeight',   opt.minPeakHeight,   ...
        'PeakThreshold',   opt.peakThreshold,   ...
        'PeakWidthLimits', opt.peakWidthLimits, ...
        'Parallel',        opt.Parallel,        ...
        'Progress',        false);

    exponent_t(:, iTime) = exp_win;
    offset_t(:, iTime)   = off_win;

    if isstruct(info_win)
    if isfield(info_win, 'error'),       fit_error_t(:, iTime) = info_win.error;       end
    if isfield(info_win, 'peak_params'), peaks_t(:, iTime)     = info_win.peak_params; end
    if isfield(info_win, 'knee'),        knee_t(:, iTime)      = info_win.knee;        end
    end

    % Optional: subtract aperiodic model (both modes).
    % Ratio to the fit; the model itself comes from compute_AperiodicBandPower
    % so this stays identical to the static path in ascent_compute. Only
    % ap_model is wanted here, so the band spans freqs to keep the helper's
    % coverage warning from firing once per time bin. Knee mode: the knee
    % parameter is ap_params(2) of each channel's fit (info.knee{ch} = full
    % param vector), passed so the full Lorentzian model is used.
    if opt.correctAperiodic
        if strcmp(opt.aperiodicMode, 'knee')
            knee_vec = cellfun(@(k) k(2), info_win.knee);
            [~, ~, ap_model] = compute_AperiodicBandPower(freqs, psd_avg, ...
                                    exp_win, off_win, [freqs(1) freqs(end)], knee_vec);
        else
            [~, ~, ap_model] = compute_AperiodicBandPower(freqs, psd_avg, ...
                                    exp_win, off_win, [freqs(1) freqs(end)]);
        end
        psd_corrected_t(:, :, iTime) = psd_avg ./ ap_model;
    end

    if opt.Progress
        fprintf('  time bin %3d/%3d  (t = %.2f s)  done\n', iTime, nTimes, times(iTime));
    end
end

%% --- Step 5 (optional): SPRiNT outlier peak pruning + aperiodic refit
% Faithful port of SPRiNT's remove_outliers helper (Wilson et al., 2022,
% SPRiNT.m remove_outliers):
%   * a peak is flagged when fewer than PruneMinPeaks similar peaks
%     (|dCF| <= PruneFreqTol Hz) occur within +/-PruneTimeBins time bins;
%     SPRiNT's timeRange = maxtime * window-step, so bins here are steps.
%   * flagging is ITERATIVE (SPRiNT: "while any(remove)"): after each removal
%     pass, neighbour counts are recomputed on the survivors until stable.
%   * at every time bin whose surviving peak count changed, the aperiodic
%     model is refit on the log spectrum minus the KEPT peaks' model
%     (pruned peaks are treated as noise and not subtracted), using the
%     original exponent as initial guess for the nonlinear knee refit.
%   * fit error and the corrected PSD are regenerated at refit bins.
if opt.PrunePeaks
    if opt.Progress
        fprintf('  Step 5: pruning outlier peaks (min %d neighbours within %.1f Hz across %d bins)...\n', ...
            opt.PruneMinPeaks, opt.PruneFreqTol, opt.PruneTimeBins);
    end
    nPruned = 0;
    for iChan = 1:nChan
        % Flatten per-channel peaks: [CF(Hz) power(log10) BW(Hz) bin]
        pkList = zeros(0, 4);
        origCount = zeros(1, nTimes);
        for iTime = 1:nTimes
            pp = peaks_t{iChan, iTime};
            if isempty(pp), continue; end
            origCount(iTime) = size(pp, 1);
            pkList = [pkList; [pp, iTime * ones(size(pp, 1), 1)]]; %#ok<AGROW>
        end
        if isempty(pkList), continue; end

        % --- Iterative neighbour-count pruning (SPRiNT: while any(remove)) ---
        while true
            remove = false(size(pkList, 1), 1);
            for iP = 1:size(pkList, 1)
                % SPRiNT counts the current peak as its own neighbour
                nNb = sum(abs(pkList(:, 2) - pkList(iP, 2)) <= opt.PruneTimeBins & ...
                          abs(pkList(:, 1) - pkList(iP, 1)) <= opt.PruneFreqTol) - 1;
                if nNb < opt.PruneMinPeaks
                    remove(iP) = true;
                end
            end
            if ~any(remove), break; end
            nPruned = nPruned + sum(remove);
            pkList  = pkList(~remove, :);
            if isempty(pkList), break; end
        end

        % --- Refit aperiodic at bins whose surviving peak count changed ---
        for iTime = 1:nTimes
            if origCount(iTime) == 0
                continue;                       % never had peaks (SPRiNT skips)
            end
            keptNow = pkList(pkList(:, 2) == iTime, 1:3);
            if size(keptNow, 1) == origCount(iTime)
                continue;                       % unchanged (SPRiNT: continue)
            end

            logf = log10(freqs(:));
            logp = log10(max(psd_t(iChan, :, iTime), eps));
            logp = logp(:);

            % Peak model of the KEPT peaks only (SPRiNT subtracts survivors).
            % ASCENT Gaussians live in linear frequency (as specparam):
            % height = power (log10), sigma = BW/2 (Hz).
            peak_model = zeros(size(logf));
            for iP = 1:size(keptNow, 1)
                peak_model = peak_model + keptNow(iP, 2) .* ...
                    exp(-0.5 .* ((freqs(:) - keptNow(iP, 1)) ./ max(keptNow(iP, 3) / 2, eps)).^2);
            end
            logp_rm = logp - peak_model;

            % Aperiodic refit; previous exponent feeds the nonlinear guess
            % (SPRiNT simple_ap_fit gets aperiodic_params(end) of the last fit)
            [exp_new, off_new, ap_vec, ap_params_new] = refit_aperiodic_sliding(logf, logp_rm, ...
                opt.aperiodicMode, exponent_t(iChan, iTime));
            exponent_t(iChan, iTime) = exp_new;
            offset_t(iChan, iTime)   = off_new;
            knee_t{iChan, iTime}     = ap_params_new;   % full param vector

            % Recompute fit error (MAE, matching compute_AperiodicFit's metric)
            fit_error_t(iChan, iTime) = mean(abs(logp - (ap_vec + peak_model)));

            % Regenerate the corrected PSD at this bin from the new fit
            if opt.correctAperiodic
                if strcmp(opt.aperiodicMode, 'knee') && numel(ap_params_new) == 3
                    [~, ~, ap_model] = compute_AperiodicBandPower(freqs, ...
                        psd_t(iChan, :, iTime), exponent_t(iChan, iTime), ...
                        offset_t(iChan, iTime), [freqs(1) freqs(end)], ap_params_new(2));
                else
                    [~, ~, ap_model] = compute_AperiodicBandPower(freqs, ...
                        psd_t(iChan, :, iTime), exponent_t(iChan, iTime), ...
                        offset_t(iChan, iTime), [freqs(1) freqs(end)]);
                end
                psd_corrected_t(iChan, :, iTime) = psd_t(iChan, :, iTime) ./ ap_model;
            end

            % Store surviving peaks for this bin
            peaks_t{iChan, iTime} = keptNow;
        end
    end
    if opt.Progress
        fprintf('  pruned %d outlier peak(s).\n', nPruned);
    end
end

if opt.Progress
    fprintf('Sliding aperiodic fit complete.\n');
end
end

%% ===== Local helpers (Step 5) =====
function [exponent, offset, ap_vec, ap_params] = refit_aperiodic_sliding(logf, logp, mode, exp_guess)
% Refit the aperiodic model on a peak-removed log spectrum.
% Mirrors SPRiNT's simple_ap_fit: OLS for fixed mode; nonlinear least squares
% for knee mode, warm-started from the previous exponent (exp_guess) and the
% previous behaviour of using the first power value as the offset guess.
% Returns the parameter evaluation ap_vec on logf (log10 units) and the full
% parameter vector ap_params ([offset exponent] fixed; [offset knee exponent]
% knee), matching compute_AperiodicFit's param ordering.
logf = logf(:);
logp = logp(:);
if strcmp(mode, 'fixed')
    X      = [ones(size(logf)), -logf];
    beta   = X \ logp;
    offset   = beta(1);
    exponent = beta(2);
    ap_params = [offset; exponent];
    ap_vec   = offset - exponent .* logf;
else
    % Knee model: log10 P = offset - log10(knee + f^exponent)
    % Initial guess: fixed-mode OLS for offset/exponent, knee from midpoint
    X     = [ones(size(logf)), -logf];
    beta0 = X \ logp;
    off0  = beta0(1);
    exp0  = max(0, beta0(2));
    if isfinite(exp_guess) && exp_guess > 0
        exp0 = exp_guess;                 % SPRiNT: warm start from last fit
    end
    knee0 = max(10^median(logf) / 2, eps);
    model_fn = @(b) sum((logp - (b(1) - log10(abs(b(2)) + 10.^(b(3) .* logf)))).^2);
    opts_ls  = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 2000);
    try
        if ~isempty(which('fminsearch'))
            b = fminsearch(model_fn, [off0, knee0, exp0], opts_ls);
        else
            b = [off0, knee0, exp0];      % degenerate fallback: keep guesses
        end
    catch
        b = [off0, knee0, exp0];
    end
    offset   = b(1);
    exponent = max(0, b(3));
    ap_params = [offset; abs(b(2)); exponent];
    ap_vec   = offset - log10(abs(b(2)) + 10.^(b(3) .* logf));
end
end