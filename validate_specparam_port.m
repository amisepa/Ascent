% validate_specparam_port.m — ASCENT vs Python specparam validation (TODO 2c)
%
% Benchmarks ASCENT's compute_AperiodicFit (MATLAB) against the Python
% specparam 2.0 reference implementation on:
%   (A) 200 synthetic spectra: 1/f (exp 0.5-2.5, off 0.5-1.5) + 0-3 log-space
%       Gaussians (8-14 Hz, ~alpha-like) + pink multiplicative noise.
%   (B) 64 real channel PSDs from the sample dataset (Welch, 1-40 Hz).
%
% Matched settings: fixed mode, freqs 1-40 Hz, peak_width_limits [1 12],
% max_n_peaks 3, min_peak_height 0.05, peak_threshold 2.0 — ASCENT defaults.
%
% Outputs (written to validation_specparam/):
%   matlab_results.mat  — exponent/offset/peaks/r2/MAE per spectrum
%   (Python side writes python_results.csv; compare_specparam.m merges both,
%    computes agreement stats, and renders the supplementary figure.)

root    = 'C:\Users\ccann\Documents\MATLAB\Ascent';
eegRoot = 'C:\Users\ccann\Documents\MATLAB\eeglab';
outdir  = fullfile(root, 'validation_specparam');
if ~exist(outdir, 'dir'), mkdir(outdir); end

addpath(genpath(eegRoot));
addpath(root);
addpath(fullfile(root, 'functions'));
warning('off', 'all');
bad = {fullfile(eegRoot, 'plugins', 'Biosig3.8.5'), ...
       fullfile(eegRoot, 'plugins', 'Fieldtrip-lite250523'), ...
       fullfile(eegRoot, 'plugins', 'roiconnect', 'libs', 'mvgc_v1.0', 'utils', 'legacy')};
for k = 1:numel(bad)
    if exist(bad{k}, 'dir'), rmpath(genpath(bad{k})); end
end
% randi must resolve to the built-in (roiconnect's legacy utils shadow it)
assert(isempty(regexpi(which('randi'), 'eeglab', 'once')), 'randi shadowed by: %s', which('randi'));

%% ---------- Part A: synthetic benchmark ----------
rng(2026);
nSyn   = 200;
f_ax   = (1:0.25:40).';
nF     = numel(f_ax);
logf   = log10(f_ax);

expTrue = 0.5 + 2.0*rand(nSyn,1);
offTrue = 0.5 + 1.0*rand(nSyn,1);
nPkTrue = randi([0 3], nSyn,1);
cfTrue  = 8 + 6*rand(nSyn,3);      % peaks between 8 and 14 Hz
ampTrue = 0.2 + 0.4*rand(nSyn,3);  % log10 height above 1/f
sigTrue = 0.05 + 0.10*rand(nSyn,3);% sigma_log10 (~0.5-1.5 Hz at 10 Hz)

psdSyn = zeros(nSyn, nF);
for k = 1:nSyn
    base = 10.^(offTrue(k) - expTrue(k).*logf.');
    pk = zeros(1, nF);
    for p = 1:nPkTrue(k)
        pk = pk + 10.^ampTrue(k,p) .* exp(-0.5.*((logf.' - log10(cfTrue(k,p)))./sigTrue(k,p)).^2);
    end
    psdSyn(k,:) = base + pk;
end
% pink-ish multiplicative noise: log-normal, correlated across freq
noise = 10.^((randn(nSyn, nF) * 0.02) + bsxfun(@times, randn(nSyn,1)*0.01, ones(1,nF)));
psdSyn = psdSyn .* noise;

%% ---------- Fit with ASCENT (MATLAB) ----------
nTest   = nSyn;
expMT   = nan(nTest,1); offMT = nan(nTest,1);
r2M     = nan(nTest,1); maeM = nan(nTest,1);
npkM    = nan(nTest,1); peaksM = cell(nTest,1);

for k = 1:nTest
    [e, o, info] = compute_AperiodicFit(f_ax.', psdSyn(k,:), ...
        'AperiodicMode', 'fixed', 'FreqRange', [1 40], ...
        'MaxPeaks', 3, 'MinPeakHeight', 0.05, 'PeakThreshold', 2.0, ...
        'PeakWidthLimits', [1 12]);
    expMT(k) = e; offMT(k) = o;
    if isstruct(info)
        r2M(k)  = info.r_squared;
        maeM(k) = info.error;
        npkM(k) = size(info.peak_params{1}, 1);
        peaksM{k} = info.peak_params{1};
    end
end
fprintf('ASCENT fits done: %d/%d valid\n', sum(isfinite(expMT)), nTest);

%% ---------- Part B: real PSDs from the sample dataset ----------
scratch = fullfile(root, 'scratch_aperiodic');
EEG = pop_loadset('filename', 'ascent_sample_data.set', 'filepath', scratch);
% Welch PSD identical to ascent_compute's static path (4 s hann, 50%, nfft 2^nextpow2)
[psdReal, freqsR] = compute_psd(EEG.data, EEG.srate, ...
    'winSec', 4, 'overlap', 0.5, 'window', 'hann', 'nfft', [], 'detrend', true);
nReal = size(psdReal, 1);
fprintf('real PSDs: %d chan x %d freqs (%.2f-%.2f Hz)\n', nReal, numel(freqsR), freqsR(1), freqsR(end));

expMR = nan(nReal,1); offMR = nan(nReal,1); r2MR = nan(nReal,1); maeMR = nan(nReal,1);
for ch = 1:nReal
    [e, o, info] = compute_AperiodicFit(freqsR, psdReal(ch,:), ...
        'AperiodicMode', 'fixed', 'FreqRange', [1 40], ...
        'MaxPeaks', 3, 'MinPeakHeight', 0.05, 'PeakThreshold', 2.0, ...
        'PeakWidthLimits', [1 12]);
    expMR(ch) = e; offMR(ch) = o;
    if isstruct(info), r2MR(ch) = info.r_squared; maeMR(ch) = info.error; end
end

%% ---------- Save for Python side ----------
save(fullfile(outdir, 'matlab_results.mat'), 'f_ax', 'psdSyn', ...
    'expTrue', 'offTrue', 'nPkTrue', 'cfTrue', 'ampTrue', 'sigTrue', ...
    'expMT', 'offMT', 'r2M', 'maeM', 'npkM', 'peaksM', ...
    'freqsR', 'psdReal', 'expMR', 'offMR', 'r2MR', 'maeMR', '-v7');
fprintf('saved %s\n', fullfile(outdir, 'matlab_results.mat'));
fprintf('==== MATLAB SIDE DONE ====\n');