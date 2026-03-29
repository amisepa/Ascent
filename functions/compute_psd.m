function [psd, freqs] = compute_psd(data, fs, varargin)
% compute_psd - Compute Welch PSD (linear) for use with specparam/FOOOF.
%
% INPUTS:
%   data        - [channels x samples] or [1 x samples]
%   fs          - sampling rate (Hz)
%
% OPTIONAL NAME-VALUE:
%   'freqRange'  - [fMin fMax] Hz to retain (default: [1 40])
%   'winSec'     - segment length in seconds (default: 4; gives 0.25 Hz resolution)
%   'overlap'    - fractional overlap 0-1 (default: 0.5)
%   'window'     - window type string accepted by window() (default: 'hann')
%   'nfft'       - DFT points; [] = next power of 2 >= seg length (default: [])
%   'detrend'    - linear detrend entire signal before Welch (default: true)
%   'Progress'   - print per-channel progress to console (default: true)
%   'Parallel'   - use parfor across channels (default: true)
%
% OUTPUTS:
%   psd         - [channels x freqs] linear power spectral density (uV^2/Hz)
%   freqs       - [1 x freqs] frequency vector (Hz)
%
% NOTES:
%   - Returns one-sided linear PSD, ready for specparam/FOOOF (no dB conversion).
%   - winSec >= 2 s is recommended so that freq resolution <= 0.5 Hz, satisfying
%     specparam's requirement that peak_width_limits(1) >= 2 * freq_resolution
%     (Donoghue et al. 2020, Nature Neuroscience).
%
% USAGE:
%   [psd, freqs] = compute_psd(EEG.data, EEG.srate, 'freqRange', [1 40], ...
%       'winSec', 4, 'overlap', 0.5, 'window', 'hann', 'nfft', [], 'detrend', true);

% --- parse inputs ---
p = inputParser;
addRequired(p,  'data');
addRequired(p,  'fs',         @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'freqRange',  [1 40],    @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'winSec',     4,         @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'overlap',    0.5,       @(x) isnumeric(x) && x >= 0 && x < 1);
addParameter(p, 'window',     'hann',    @ischar);
addParameter(p, 'nfft',       [],        @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'detrend',    true,      @islogical);
addParameter(p, 'Progress',   true,     @islogical);
addParameter(p, 'Parallel',   true,     @islogical);
parse(p, data, fs, varargin{:});
opt = p.Results;

% --- build window vector ---
segLen   = round(opt.winSec * fs);              % samples per segment
noverlap = round(opt.overlap * segLen);          % overlap in samples
win      = window(str2func(opt.window), segLen); % e.g. hann(segLen)

% --- nfft: zero-pad to next power of 2 for smoother curve ---
if isempty(opt.nfft)
    nfft = 2^nextpow2(segLen);
else
    nfft = opt.nfft;
end

% validate: nfft must be >= segLen
if nfft < segLen
    warning('nfft (%d) < segLen (%d); setting nfft = segLen.', nfft, segLen);
    nfft = segLen;
end

% --- warn if freq resolution is too coarse for specparam ---
freqRes = fs / segLen;
if freqRes > 0.5
    warning(['Frequency resolution is %.3f Hz (winSec = %.1f s). ' ...
             'specparam recommends freq_res <= 0.5 Hz (winSec >= 2 s) ' ...
             'to avoid overfitting noise as peaks.'], freqRes, opt.winSec);
end

fprintf(['PSD settings: window = %s, segLen = %d samples (%.2f s), ' ...
         'overlap = %.0f%%, nfft = %d, freq_res = %.4f Hz\n'], ...
        opt.window, segLen, segLen/fs, opt.overlap*100, nfft, freqRes);

% --- optionally linear-detrend entire signal before Welch ---
if opt.detrend
    data = detrend(data, 'linear');  % operates row-wise: [channels x samples]
end

% --- compute Welch PSD per channel ---
nChan = size(data, 1);

% first call to get frequency vector
[pxx_tmp, freqs_full] = pwelch(data(1,:)', win, noverlap, nfft, fs);

% frequency mask
fMask = freqs_full >= opt.freqRange(1) & freqs_full <= opt.freqRange(2);
freqs = freqs_full(fMask)';
psd   = nan(nChan, sum(fMask));

if opt.Parallel

    % parfor: progress tracking via per-channel fprintf (no in-place overwrite,
    % parfor workers write asynchronously so order is not guaranteed)
    if opt.Progress
        fprintf('Computing PSD in parallel (%d channels)...\n', nChan);
    end
    parfor iChan = 1:nChan
        pxx = pwelch(data(iChan,:)', win, noverlap, nfft, fs);  %#ok<PFBNS>
        psd(iChan,:) = pxx(fMask)';                             %#ok<PFBNS>
        if opt.Progress
            fprintf('  channel %d/%d done\n', iChan, nChan);
        end
    end
    if opt.Progress, fprintf('PSD done.\n'); end

else

    if opt.Progress, progressbar('Computing PSD'); end
    for iChan = 1:nChan
        pxx = pwelch(data(iChan,:)', win, noverlap, nfft, fs);
        psd(iChan,:) = pxx(fMask)';
        if opt.Progress, progressbar(iChan/nChan); end
    end

end