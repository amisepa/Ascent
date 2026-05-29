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
%                       fixed mode only, disabled automatically for knee mode)
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
%
% NOTES
%   Architecture follows SPRiNT (Wilson et al., 2022):
%     Step 1 - Compute Hann-windowed FFT power for every short window.
%     Step 2 - At each time bin, average the power of nAvg consecutive
%              windows to obtain a smoothed single-sided PSD.
%     Step 3 - Fit the specparam aperiodic model to that averaged PSD.
%     Step 4 - (optional) Subtract the aperiodic model.
%
%   Temporal resolution = winSec * (1 - winOverlap) seconds per step.
%   Frequency resolution = 1/winSec Hz.
%   Total output time bins = nWins - nAvg + 1, where
%     nWins = floor((nSamples - winSamples) / stepSamples) + 1.
%
%   Aperiodic correction is implemented for fixed mode only.
%   Model: log10(P_corr) = log10(P) - (offset - exponent * log10(f))
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
addParameter(p, 'Parallel',        true,    @islogical);
addParameter(p, 'Progress',        true,    @islogical);
parse(p, data, fs, varargin{:});
opt = p.Results;

opt.aperiodicMode = lower(char(opt.aperiodicMode));
if ~ismember(opt.aperiodicMode, {'fixed','knee'})
    error('compute_AperiodicFit_sliding:badMode', ...
        'aperiodicMode must be ''fixed'' or ''knee''.');
end
if strcmp(opt.aperiodicMode, 'knee') && opt.correctAperiodic
    warning('compute_AperiodicFit_sliding:kneeCorrection', ...
        'Aperiodic correction not implemented for knee mode. Disabling.');
    opt.correctAperiodic = false;
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

log_freqs = log10(freqs);   % pre-compute for correction step

%% --- Steps 2-4: average windows, fit specparam, optional correction ---
if opt.Progress
    fprintf('  Step 2: fitting aperiodic model to %d time bins...\n', nTimes);
end

fit_error_t = nan(nChan, nTimes);
peaks_t     = cell(nChan, nTimes);

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
    if isfield(info_win, 'error'),  fit_error_t(:, iTime) = info_win.error;  end
    if isfield(info_win, 'peaks'),  peaks_t(:, iTime)     = info_win.peaks;  end
    end

    % Optional: subtract aperiodic model (fixed mode only)
    % Model: log10(P_corr) = log10(P) - (offset - exponent * log10(f))
    if opt.correctAperiodic
        psd_corr = nan(nChan, nFreqs);
        for iChan = 1:nChan
            if isnan(exp_win(iChan)) || isnan(off_win(iChan)), continue; end
            ap_model = off_win(iChan) - exp_win(iChan) .* log_freqs;
            psd_corr(iChan, :) = 10.^(log10(psd_avg(iChan,:)) - ap_model);
        end
        psd_corrected_t(:, :, iTime) = psd_corr;
    end

    if opt.Progress
        fprintf('  time bin %3d/%3d  (t = %.2f s)  done\n', iTime, nTimes, times(iTime));
    end
end

if opt.Progress
    fprintf('Sliding aperiodic fit complete.\n');
end
end