function [exponent_t, offset_t, times, freqs, psd_t, psd_corrected_t] = ...
        compute_AperiodicFit_sliding(data, fs, varargin)
% compute_AperiodicFit_sliding  Time-varying aperiodic estimation via a
%   sliding-window Welch PSD + specparam fit.  At each time step, a Welch
%   PSD is estimated from the local data segment and the aperiodic model is
%   fitted to that PSD, yielding continuous exponent(t) and offset(t) traces.
%   Optionally the aperiodic component is subtracted to produce a corrected
%   PSD(t) that isolates periodic activity.
%
% INPUTS
%   data   - [nChan x nSamples], µV
%   fs     - sampling rate (Hz)
%
% OPTIONAL NAME-VALUE
%   'slidWinSec'      - outer sliding window length in s (default: 4)
%   'slidStepSec'     - step between successive windows in s (default: 1)
%   'freqRange'       - [fMin fMax] Hz for both PSD and fit (default: [1 40])
%   'psdWinSec'       - Welch segment length within each window in s
%                       (default: slidWinSec / 2; must be <= slidWinSec)
%   'psdOverlap'      - Welch fractional overlap 0-1 (default: 0.5)
%   'aperiodicMode'   - 'fixed' | 'knee' (default: 'fixed')
%   'maxPeaks'        - max Gaussian peaks per fit (default: 6)
%   'minPeakHeight'   - minimum peak height above aperiodic, log10 units (default: 0.05)
%   'peakThreshold'   - peak detection threshold in SD of flattened spectrum (default: 2.0)
%   'peakWidthLimits' - [min max] peak width in Hz (default: [1 12])
%   'correctAperiodic'- subtract aperiodic model from PSD (default: true;
%                       fixed mode only — disabled automatically for knee mode)
%   'Parallel'        - parfor over channels inside each window (default: true)
%   'Progress'        - print per-window progress (default: true)
%
% OUTPUTS
%   exponent_t      - [nChan x nTimes]  aperiodic exponent per window
%   offset_t        - [nChan x nTimes]  aperiodic offset per window (log10 units)
%   times           - [1 x nTimes]      centre time of each window (s)
%   freqs           - [1 x nFreqs]      frequency vector (Hz)
%   psd_t           - [nChan x nFreqs x nTimes]  linear PSD per window (µV²/Hz)
%   psd_corrected_t - [nChan x nFreqs x nTimes]  aperiodic-corrected PSD;
%                     empty ([]) if correctAperiodic = false
%
% NOTES
%   The outer window loop is serial; parallelism is applied inside each window
%   across channels via compute_psd and compute_AperiodicFit.  Nesting parfor
%   loops is not supported in MATLAB and is avoided here.
%
%   Aperiodic correction is implemented for fixed mode only (same restriction
%   as in the static pipeline).  Model: log10(P_corr) = log10(P) - (offset - exponent * log10(f))
%   converted back to linear scale.
%
%   Recommended minimum slidWinSec to obtain a reliable Welch estimate:
%   with psdWinSec = slidWinSec/2 and psdOverlap = 0.5, each window yields
%   approximately 3 Welch segments.  Longer windows improve spectral
%   smoothness at the cost of temporal resolution.
%
% REFERENCES
%   Donoghue T, et al. (2020). Parameterizing neural power spectra into
%     periodic and aperiodic components. Nature Neuroscience, 23, 1655-1665.
%   Kałamała P, et al. (2026). How to improve the reliability of aperiodic
%     parameter estimates in M/EEG. [preprint]
%   Wilson LE, et al. (2022). Time-resolved parameterization of aperiodic
%     and periodic brain activity. eLife, 11, e77348.
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Ascent Plugin - Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

%% --- Parse inputs ---
p = inputParser;
addRequired(p,  'data', @(x) isnumeric(x) && ismatrix(x));
addRequired(p,  'fs',   @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'slidWinSec',      4,       @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'slidStepSec',     1,       @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'freqRange',       [1 40],  @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'psdWinSec',       [],      @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(p, 'psdOverlap',      0.5,     @(x) isnumeric(x) && x >= 0 && x < 1);
addParameter(p, 'aperiodicMode',   'fixed', @(s) ischar(s) || isstring(s));
addParameter(p, 'maxPeaks',        6,       @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'minPeakHeight',   0.05,    @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'peakThreshold',   2.0,     @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'peakWidthLimits', [1 12],  @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'correctAperiodic',true,    @islogical);
addParameter(p, 'Parallel',        true,    @islogical);
addParameter(p, 'Progress',        true,    @islogical);
parse(p, data, fs, varargin{:});
opt = p.Results;

opt.aperiodicMode = lower(char(opt.aperiodicMode));
if ~ismember(opt.aperiodicMode, {'fixed','knee'})
    error('compute_AperiodicFit_sliding:badMode', ...
        'aperiodicMode must be ''fixed'' or ''knee''.');
end

% Default psdWinSec = half the sliding window
if isempty(opt.psdWinSec)
    opt.psdWinSec = opt.slidWinSec / 2;
end
if opt.psdWinSec > opt.slidWinSec
    warning('compute_AperiodicFit_sliding:psdWinTooLong', ...
        'psdWinSec (%.1f s) > slidWinSec (%.1f s); clipping to %.1f s.', ...
        opt.psdWinSec, opt.slidWinSec, opt.slidWinSec);
    opt.psdWinSec = opt.slidWinSec;
end

% Correction not available for knee mode
if strcmp(opt.aperiodicMode, 'knee') && opt.correctAperiodic
    warning('compute_AperiodicFit_sliding:kneeCorrection', ...
        'Aperiodic correction not implemented for knee mode. Disabling.');
    opt.correctAperiodic = false;
end

%% --- Window indices ---
[nChan, nSamp] = size(data);
slidWinSamp  = round(opt.slidWinSec  * fs);
slidStepSamp = round(opt.slidStepSec * fs);

starts = 1 : slidStepSamp : (nSamp - slidWinSamp + 1);
nTimes = numel(starts);
if nTimes == 0
    error('compute_AperiodicFit_sliding:dataTooShort', ...
        'Data (%.2f s) is shorter than the sliding window (%.1f s).', ...
        nSamp/fs, opt.slidWinSec);
end

% Centre time of each window (s)
times = ((starts - 1) + slidWinSamp/2) / fs;

%% --- Frequency vector: determined from a single dummy PSD call ---
[~, freqs] = compute_psd(data(1,:), fs, ...
    'freqRange', opt.freqRange,  ...
    'winSec',    opt.psdWinSec,  ...
    'overlap',   opt.psdOverlap, ...
    'detrend',   true,           ...
    'Progress',  false,          ...
    'Parallel',  false);
nFreqs = numel(freqs);

%% --- Pre-allocate outputs ---
exponent_t = nan(nChan, nTimes);
offset_t   = nan(nChan, nTimes);
psd_t      = nan(nChan, nFreqs, nTimes);
if opt.correctAperiodic
    psd_corrected_t = nan(nChan, nFreqs, nTimes);
else
    psd_corrected_t = [];
end

if opt.Progress
    fprintf(['Sliding aperiodic fit: %d windows x %d channels\n' ...
             '  slidWin=%.1fs  step=%.1fs  psdWin=%.1fs  overlap=%.0f%%\n'], ...
        nTimes, nChan, opt.slidWinSec, opt.slidStepSec, ...
        opt.psdWinSec, opt.psdOverlap*100);
end

%% --- Outer window loop (serial; channel parallelism handled inside) ---
log_freqs = log10(freqs);   % pre-compute for correction step

for iWin = 1:nTimes
    i1  = starts(iWin);
    i2  = i1 + slidWinSamp - 1;
    seg = data(:, i1:i2);   % [nChan x slidWinSamp]

    %% Step 1: Welch PSD for this window
    [psd_win, ~] = compute_psd(seg, fs, ...
        'freqRange', opt.freqRange,   ...
        'winSec',    opt.psdWinSec,   ...
        'overlap',   opt.psdOverlap,  ...
        'detrend',   true,            ...
        'Progress',  false,           ...
        'Parallel',  opt.Parallel);

    %% Step 2: Aperiodic fit
    [exp_win, off_win, ~] = compute_AperiodicFit(freqs, psd_win, ...
        'FreqRange',       opt.freqRange,       ...
        'AperiodicMode',   opt.aperiodicMode,   ...
        'MaxPeaks',        opt.maxPeaks,        ...
        'MinPeakHeight',   opt.minPeakHeight,   ...
        'PeakThreshold',   opt.peakThreshold,   ...
        'PeakWidthLimits', opt.peakWidthLimits, ...
        'Parallel',        opt.Parallel,        ...
        'Progress',        false);

    exponent_t(:, iWin) = exp_win;
    offset_t(:, iWin)   = off_win;
    psd_t(:, :, iWin)   = psd_win;

    %% Step 3 (optional): aperiodic correction — fixed mode only
    % Model: log10(P_corr) = log10(P) - (offset - exponent * log10(f))
    if opt.correctAperiodic
        psd_corr = nan(nChan, nFreqs);
        for iChan = 1:nChan
            if isnan(exp_win(iChan)) || isnan(off_win(iChan)), continue; end
            ap_model = off_win(iChan) - exp_win(iChan) .* log_freqs;
            psd_corr(iChan, :) = 10.^(log10(psd_win(iChan,:)) - ap_model);
        end
        psd_corrected_t(:, :, iWin) = psd_corr;
    end

    if opt.Progress
        fprintf('  window %3d/%3d  [%6.2f - %6.2f s]  done\n', ...
            iWin, nTimes, (i1-1)/fs, i2/fs);
    end
end

if opt.Progress
    fprintf('Sliding aperiodic fit complete.\n');
end
end
