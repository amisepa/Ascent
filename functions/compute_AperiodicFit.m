function [exponent, offset, info] = compute_AperiodicFit(freqs, psd, varargin)
% Aperiodic exponent and offset via spectral parameterization (FOOOF/specparam).
% Pure MATLAB reimplementation — no Python, no external toolboxes required.
%
%   [exponent, offset, info] = compute_AperiodicFit(freqs, psd, ...
%       'FreqRange', [1 40], 'AperiodicMode', 'fixed', ...
%       'MaxPeaks', 3, 'MinPeakHeight', 0.05, 'PeakThreshold', 2.0, ...
%       'PeakWidthLimits', [1 12], ...
%       'Parallel', true, 'Progress', true)
%
% Inputs
%   freqs   : frequency vector [1 x nFreqs], linear spacing, Hz
%   psd     : power spectra [nChan x nFreqs], linear scale (not log)
%             Rows are channels/components; matches EEGLAB convention.
%   Name-Value pairs:
%     'FreqRange'       : [fmin fmax] Hz to fit (default [1 40])
%     'AperiodicMode'   : 'fixed' (default) | 'knee'
%                         'fixed'  — log P = offset - exponent * log F
%                         'knee'   — log P = offset - log(knee + F^exponent)
%     'MaxPeaks'        : maximum number of Gaussian peaks to fit (default 3)
%     'MinPeakHeight'   : minimum peak height above aperiodic, in log10 power
%                         units (default 0.05)
%     'PeakThreshold'   : peak detection threshold in units of SD of the
%                         flattened spectrum (default 2.0; FOOOF default = 2.0)
%     'PeakWidthLimits' : [min max] peak width in Hz; width = 2*sigma
%                         (default [1 12])
%     'Parallel'        : true|false, parfor over channels (default true)
%     'Progress'        : true|false, print per-channel output (default true)
%
% Outputs
%   exponent : aperiodic exponent per channel [nChan x 1]
%              steeper 1/f → larger exponent
%   offset   : aperiodic offset per channel [nChan x 1] (log10 units)
%   info     : per-channel diagnostics struct:
%                .knee        – knee parameter (NaN if AperiodicMode='fixed')
%                .ap_fit      – fitted aperiodic spectrum [nFit x 1], log10
%                .flat_spec   – flattened (aperiodic-removed) spectrum
%                .peak_params – [nPeaks x 3]: [CF, power, BW] per peak
%                               CF = center freq (Hz), power = height above
%                               aperiodic (log10), BW = 2*sigma (Hz)
%                .gauss_params– [nPeaks x 3]: underlying Gaussian [CF (Hz),
%                               height (log10), sigma (Hz)], linear-frequency
%                               Gaussians as in specparam
%                .r_squared   – R² of full model vs. log10 PSD
%                .error       – MAE of full model vs. log10 PSD
%                .freqs_used  – frequency vector used for fitting
%                .flags       – struct: fitFailed, highError, noConverge
%              plus .params   – copy of options used
%
% Algorithm (Donoghue et al., 2020, Nat Neurosci; mirrors specparam 2.0
% SpectralFitAlgorithm._fit step by step)
%   1. Trim spectrum to FreqRange; log10 power.
%   2. Robust initial aperiodic fit: simple fit, flatten, clip negatives to
%      0, refit on points <= 0.025th percentile of the clipped residual
%      (i.e. points at/below the first fit), as specparam _robust_ap_fit.
%   3. Subtract initial aperiodic → flattened spectrum.
%   4. Greedy peak search on the flattened spectrum (Gaussians in LINEAR Hz,
%      guess width from FWHM, guess subtracted), drop edge / overlapping
%      guesses, then one joint bounded fit of all Gaussians.
%   5. Simple (non-robust) aperiodic re-fit on the peak-removed spectrum.
%   6. Final full model = aperiodic + sum of Gaussians; compute R² and MAE.
%
% The aperiodic model in log10 space:
%   fixed: log10_P(f) = offset - exponent * log10(f)
%   knee:  log10_P(f) = offset - log10(knee + f^exponent)
%
% References
%   Donoghue T, Haller M, Peterson EJ, et al. (2020).
%     Parameterizing neural power spectra into periodic and aperiodic
%     components. Nature Neuroscience, 23, 1655-1665.
%     https://doi.org/10.1038/s41593-020-00744-x
%
% Note on comparability
%   Results should be numerically close to Python specparam/FOOOF 1.x
%   with equivalent settings, as this follows the same algorithm. Minor
%   differences may arise from optimizer tolerances and Gaussian init.
%   The Brainstorm and FieldTrip toolboxes use the same MATLAB reimplementation
%   approach; cite Donoghue et al. (2020) regardless of which port is used.
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Ascent Plugin - Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

% Parse inputs
p = inputParser;
p.addRequired('freqs', @(x) isnumeric(x) && isvector(x));
p.addRequired('psd',   @(x) isnumeric(x) && ismatrix(x));
p.addParameter('FreqRange',       [1 40],      @(x) isnumeric(x) && numel(x)==2);
p.addParameter('AperiodicMode',   'fixed',     @(s) ischar(s) || isstring(s));
p.addParameter('MaxPeaks',        3,           @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('MinPeakHeight',   0.05,        @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('PeakThreshold',   2.0,         @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('PeakWidthLimits', [1 12],     @(x) isnumeric(x) && numel(x)==2);
p.addParameter('Parallel',        true,        @(x) islogical(x) && isscalar(x));
p.addParameter('Progress',        true,        @(x) islogical(x) && isscalar(x));
p.parse(freqs, psd, varargin{:});
opts = p.Results;

opts.AperiodicMode = lower(char(opts.AperiodicMode));
if ~ismember(opts.AperiodicMode, {'fixed','knee'})
    error('compute_AperiodicFit:badMode', ...
        'AperiodicMode must be ''fixed'' or ''knee''.');
end

% Shape inputs 
freqs = freqs(:).';              % row
if size(psd,2) ~= numel(freqs)
    if size(psd,1) == numel(freqs)
        psd = psd.';             % -> [nChan x nFreqs]
    else
        error('compute_AperiodicFit:dimMismatch', ...
            'psd dimensions do not match freqs length.');
    end
end
[nchan, ~] = size(psd);

% Frequency trimming
fmask = freqs >= opts.FreqRange(1) & freqs <= opts.FreqRange(2);
if sum(fmask) < 4
    error('compute_AperiodicFit:tooFewFreqs', ...
        'FreqRange [%.1f %.1f] yields fewer than 4 frequency bins.', ...
        opts.FreqRange(1), opts.FreqRange(2));
end
f_use   = freqs(fmask);           % [1 x nFit]
logf    = log10(f_use(:));        % [nFit x 1]

% Warn if peak width lower limit < 2x frequency resolution
freq_res = mean(diff(f_use));
if opts.PeakWidthLimits(1) < 2 * freq_res
    warning('compute_AperiodicFit:narrowPeakWidth', ...
        ['Lower PeakWidthLimits (%.2f Hz) < 2x frequency resolution (%.2f Hz). ' ...
         'Consider widening to avoid fitting noise.'], ...
        opts.PeakWidthLimits(1), 2*freq_res);
end

% Pre-allocate outputs 
exponent = nan(nchan,1);
offset   = nan(nchan,1);
info             = struct();
info.knee        = cell(nchan,1);
info.ap_fit      = cell(nchan,1);
info.flat_spec   = cell(nchan,1);
info.peak_params = cell(nchan,1);
info.gauss_params= cell(nchan,1);
info.r_squared   = nan(nchan,1);
info.error       = nan(nchan,1);
info.freqs_used  = f_use;
info.flags       = cell(nchan,1);
info.params      = opts;

% Progress header
if opts.Progress
    parStr = 'off'; if opts.Parallel, parStr = 'on'; end
    fprintf('AperiodicFit: %d channel(s) | mode=%s | freqs=[%.0f %.0f] Hz | parallel=%s\n', ...
        nchan, opts.AperiodicMode, opts.FreqRange(1), opts.FreqRange(2), parStr);
end

% Iterate channels
if opts.Parallel && ~isempty(ver('parallel'))

    knee_all   = cell(nchan,1);
    apfit_all  = cell(nchan,1);
    flat_all   = cell(nchan,1);
    pp_all     = cell(nchan,1);
    gp_all     = cell(nchan,1);
    rsq_all    = nan(nchan,1);
    err_all    = nan(nchan,1);
    flags_all  = cell(nchan,1);

    parfor ch = 1:nchan
        logp = log10(max(psd(ch, fmask), eps));  % [1 x nFit], guard zero
        [exp_ch, off_ch, chInfo] = fooof_single(f_use(:), logf, logp(:), opts);
        exponent(ch) = exp_ch;
        offset(ch)   = off_ch;
        knee_all{ch}  = chInfo.knee;
        apfit_all{ch} = chInfo.ap_fit;
        flat_all{ch}  = chInfo.flat_spec;
        pp_all{ch}    = chInfo.peak_params;
        gp_all{ch}    = chInfo.gauss_params;
        rsq_all(ch)   = chInfo.r_squared;
        err_all(ch)   = chInfo.error;
        flags_all{ch} = chInfo.flags;
        if opts.Progress
            fprintf('  ch %3d/%3d: exp=%6.4f  offset=%7.4f  R²=%6.4f  nPeaks=%d\n', ...
                ch, nchan, exp_ch, off_ch, chInfo.r_squared, size(chInfo.peak_params,1));
        end
    end

    info.knee         = knee_all;
    info.ap_fit       = apfit_all;
    info.flat_spec    = flat_all;
    info.peak_params  = pp_all;
    info.gauss_params = gp_all;
    info.r_squared    = rsq_all;
    info.error        = err_all;
    info.flags        = flags_all;

else
    useWB = opts.Progress && usejava('desktop');
    hWB = [];
    if useWB
        try, hWB = waitbar(0,'Fitting aperiodic component...','Name','compute_AperiodicFit');
        catch, hWB = []; end
    end

    for ch = 1:nchan
        logp = log10(max(psd(ch, fmask), eps));
        [exponent(ch), offset(ch), chInfo] = fooof_single(f_use(:), logf, logp(:), opts);
        info.knee{ch}         = chInfo.knee;
        info.ap_fit{ch}       = chInfo.ap_fit;
        info.flat_spec{ch}    = chInfo.flat_spec;
        info.peak_params{ch}  = chInfo.peak_params;
        info.gauss_params{ch} = chInfo.gauss_params;
        info.r_squared(ch)    = chInfo.r_squared;
        info.error(ch)        = chInfo.error;
        info.flags{ch}        = chInfo.flags;
        if opts.Progress
            fprintf('  ch %3d/%3d: exp=%6.4f  offset=%7.4f  R²=%6.4f  nPeaks=%d\n', ...
                ch, nchan, exponent(ch), offset(ch), chInfo.r_squared, ...
                size(chInfo.peak_params,1));
            if ~isempty(hWB) && isvalid(hWB)
                try, waitbar(ch/nchan, hWB, sprintf('Aperiodic fit... (%d/%d)', ch, nchan)); end
            end
        end
    end
    if ~isempty(hWB) && isvalid(hWB), try, close(hWB); end, end
end
end


%% Helpers

% Single-spectrum worker
function [exp_out, off_out, out] = fooof_single(f_lin, logf, logp, opts)
% f_lin: [nFit x 1] frequencies (Hz); logf, logp: [nFit x 1] log10 values.
% Step order and rules follow specparam 2.0 SpectralFitAlgorithm._fit.

out = init_out();

% Initial aperiodic fit, robust to peaks (specparam _robust_ap_fit)
ap0 = fit_aperiodic_robust(logf, logp, opts.AperiodicMode);
if any(~isfinite(ap0))
    out.flags.fitFailed = true;
    out.knee = ap0;   % param-vector shape ([NaN NaN] or [NaN NaN NaN])
    exp_out = NaN; off_out = NaN;
    return
end

% Flatten spectrum by subtracting initial aperiodic
ap_spec0  = aperiodic_model(logf, ap0, opts.AperiodicMode);
flat_spec = logp - ap_spec0;

% Peak search + joint Gaussian fit, Gaussians in LINEAR frequency (Hz)
gauss_params = find_peaks(f_lin, flat_spec, opts);  % [nPeaks x 3]: [cf_Hz height std_Hz]
gauss_fit    = eval_gaussians(reshape(gauss_params.', 1, []), f_lin, size(gauss_params,1));

% Final aperiodic fit: simple (non-robust) fit on the peak-removed spectrum,
% also when no peak was found (specparam _simple_ap_fit on _spectrum_peak_rm)
ap_params = fit_aperiodic(logf, logp - gauss_fit, opts.AperiodicMode);
if any(~isfinite(ap_params)), ap_params = ap0; end

% Compute final full model and goodness of fit
ap_fit   = aperiodic_model(logf, ap_params, opts.AperiodicMode);
full_fit = ap_fit + gauss_fit;
ss_res = sum((logp - full_fit).^2);
ss_tot = sum((logp - mean(logp)).^2);
r2  = 1 - ss_res / max(ss_tot, eps);
mae = mean(abs(logp - full_fit));

% Peak params (specparam 'log_sub' / 'full_width' converters):
% CF = Gaussian centre (Hz); power = full model minus aperiodic at the bin
% nearest CF (= periodic model there, log10); BW = 2*std (Hz)
peak_params = zeros(size(gauss_params));
for pk = 1:size(gauss_params,1)
    [~, ind] = min(abs(f_lin - gauss_params(pk,1)));
    peak_params(pk,:) = [gauss_params(pk,1), gauss_fit(ind), 2 * gauss_params(pk,3)];
end

% Pack outputs
exp_out = ap_params(end);    % last param is always the exponent
off_out = ap_params(1);
out.knee         = ap_params;   % full param vector; knee is ap_params(2) if knee mode
out.ap_fit       = ap_fit;
out.flat_spec    = flat_spec;
out.peak_params  = peak_params;
out.gauss_params = gauss_params;
out.r_squared    = r2;
out.error        = mae;
out.flags.fitFailed  = false;
out.flags.highError  = mae > 0.1;
out.flags.noConverge = any(~isfinite(ap_params));
end


% Aperiodic model
function ap = aperiodic_model(logf, params, mode)
% params for 'fixed': [offset, exponent]
% params for 'knee' : [offset, knee, exponent]
% All in log10 space; logf is log10(Hz)
if strcmp(mode,'fixed')
    ap = params(1) - params(2) .* logf;
else
    % knee model: log10_P = offset - log10(knee + 10^(exponent * logf))
    ap = params(1) - log10(abs(params(2)) + 10.^(params(3) .* logf));
end
end


% Robust aperiodic fitting (specparam _robust_ap_fit)
function params = fit_aperiodic_robust(logf, logp, mode)
% 1) simple fit; 2) flatten, clip negative residuals to 0; 3) keep points at
% or below the 0.025th percentile of the clipped residual (numpy 'linear'
% percentile) - in practice every point on/below the first fit, so peaks
% cannot pull the line up; 4) refit on those points, warm-started.
logf = logf(:); logp = logp(:);
p0 = fit_aperiodic(logf, logp, mode);
if any(~isfinite(p0)), params = p0; return; end
flat = logp - aperiodic_model(logf, p0, mode);
flat(flat < 0) = 0;
s    = sort(flat);
pos  = 0.025/100 * (numel(s) - 1);          % numpy.percentile, linear interp
lo   = floor(pos);
thr  = s(lo+1) + (pos - lo) * (s(min(lo+2, numel(s))) - s(lo+1));
mask = flat <= thr;
if nnz(mask) < numel(p0) + 1, params = p0; return; end
params = fit_aperiodic(logf(mask), logp(mask), mode, p0);
if any(~isfinite(params)), params = p0; end
end


% Aperiodic fitting
function params = fit_aperiodic(logf, logp, mode, x0)
% Least-squares fit of the aperiodic model in log10 space (specparam
% _simple_ap_fit). Fixed mode is linear in log-log -> exact OLS. Knee mode
% uses lsqcurvefit if available, otherwise fminsearch; x0 optional guess.

logf = logf(:); logp = logp(:);

if strcmp(mode, 'fixed')
    X      = [ones(size(logf)) -logf];
    beta   = X \ logp;
    params = [beta(1); beta(2)];

else
    % Knee mode: nonlinear — use lsqcurvefit or fminsearch
    % Initial guess: run fixed first for offset/exponent, knee = median power
    X       = [ones(size(logf)) -logf];
    beta0   = X \ logp;
    if nargin < 4 || isempty(x0)
        off0    = beta0(1);
        exp0    = max(0, beta0(2));
        knee0   = 10^median(logf) / 2;
        x0      = [off0, knee0, exp0];
    else
        x0      = [x0(1), max(0, x0(2)), max(0, x0(3))];
    end
    lb      = [-inf, 0,   0   ];
    ub      = [ inf, inf, inf ];
    model_fn = @(b, lf) b(1) - log10(abs(b(2)) + 10.^(b(3) .* lf));
    opts_ls  = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',5000);
    try
        if ~isempty(which('lsqcurvefit'))
            params_v = lsqcurvefit(model_fn, x0, logf, logp, lb, ub, opts_ls);
        else
            cost = @(b) sum((logp - model_fn(b, logf)).^2);
            params_v = fminsearch(cost, x0, opts_ls);
        end
        params = params_v(:);
    catch
        % Fallback: return fixed-mode fit with knee=0
        params = [beta0(1); 0; max(0,beta0(2))];
    end
end
end


% Peak detection (specparam _fit_peaks, Gaussians in linear Hz)
function gauss_params = find_peaks(f, flat_spec, opts)
% Returns [nPeaks x 3]: [CF (Hz), height (log10 power), std (Hz)], by CF.

f         = f(:);
flat_spec = flat_spec(:);
freq_res  = f(2) - f(1);
std_lim   = opts.PeakWidthLimits / 2;           % 2-sided BW -> 1-sided std
f_range   = [f(1) f(end)];

% 1) Greedy search: take the max, stop on relative/absolute height, guess
%    the width from the FWHM, subtract the GUESS Gaussian, repeat.
guess     = zeros(0,3);
flat_iter = flat_spec;
while size(guess,1) < opts.MaxPeaks
    [max_height, max_ind] = max(flat_iter);
    if max_height <= opts.PeakThreshold * std(flat_iter, 1), break; end  % numpy std (ddof=0)
    if ~(max_height > opts.MinPeakHeight), break; end

    fwhm = estimate_fwhm(flat_iter, max_ind, freq_res);
    if isnan(fwhm)
        g_std = mean(opts.PeakWidthLimits);      % specparam fallback
    else
        g_std = fwhm / (2 * sqrt(2 * log(2)));
    end
    if g_std < std_lim(1), g_std = std_lim(1); end
    if g_std > std_lim(2), g_std = std_lim(1); end   % sic: specparam resets to lower limit

    cur = [f(max_ind), max_height, g_std];
    guess(end+1,:) = cur; %#ok<AGROW>
    flat_iter = flat_iter - eval_gaussians(cur, f, 1);
end

% 2) Drop guesses too close to the edges (|CF - edge| <= 1 std) ...
if ~isempty(guess)
    keep  = abs(guess(:,1) - f_range(1)) > guess(:,3) & ...
            abs(guess(:,1) - f_range(2)) > guess(:,3);
    guess = guess(keep,:);
end
% ... and the lower of any adjacent pair overlapping at +/-0.75 std
if size(guess,1) > 1
    guess = sortrows(guess, 1);
    bnds  = [guess(:,1) - 0.75*guess(:,3), guess(:,1) + 0.75*guess(:,3)];
    drop  = false(size(guess,1),1);
    for k = 1:size(guess,1)-1
        if bnds(k,2) > bnds(k+1,1)
            if guess(k,2) <= guess(k+1,2), drop(k) = true; else, drop(k+1) = true; end
        end
    end
    guess = guess(~drop,:);
end

if isempty(guess)
    gauss_params = zeros(0,3);
    return
end

% 3) Joint bounded fit of all guesses on the whole flattened spectrum:
%    CF within +/- 2*1.5*std of guess (clipped to range), height >= 0,
%    std within PeakWidthLimits/2.
nP = size(guess,1);
x0 = reshape(guess.', 1, []);                  % [cf1 h1 s1 cf2 h2 s2 ...]
lb = zeros(1, 3*nP); ub = zeros(1, 3*nP);
for k = 1:nP
    lcf = guess(k,1) - 2*1.5*guess(k,3);
    hcf = guess(k,1) + 2*1.5*guess(k,3);
    lb(3*k-2:3*k) = [max(lcf, f_range(1)), 0,   std_lim(1)];
    ub(3*k-2:3*k) = [min(hcf, f_range(2)), Inf, std_lim(2)];
end
fun     = @(b, x) eval_gaussians(b, x, nP);
opts_ls = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',3000, ...
                   'MaxFunEvals',5000);
try
    if ~isempty(which('lsqcurvefit'))
        x_fit = lsqcurvefit(fun, x0, f, flat_spec, lb, ub, opts_ls);
    else
        cost  = @(b) sum((flat_spec - fun(b, f)).^2);
        x_fit = min(max(fminsearch(cost, x0, opts_ls), lb), ub);
    end
catch
    x_fit = x0;                                % keep guesses
end
gauss_params = sortrows(reshape(x_fit, 3, []).', 1);
end


function fwhm = estimate_fwhm(flat, peak_ind, freq_res)
% specparam estimate_fwhm: half-max crossing on each side, keep the shorter
% side (robust to overlapping neighbours); NaN if neither side crosses.
half = 0.5 * flat(peak_ind);
le = find(flat(2:peak_ind-1) <= half, 1, 'last');       % numpy range stops at index 1
if ~isempty(le), le = le + 1; end
ri = find(flat(peak_ind+1:end) <= half, 1, 'first');
if ~isempty(ri), ri = ri + peak_ind; end
sides = [abs(le - peak_ind), abs(ri - peak_ind)];
if isempty(sides)
    fwhm = NaN;
else
    fwhm = min(sides) * 2 * freq_res;
end
end


function y = eval_gaussians(params, x, n)
% Sum of n Gaussians; params = [mu1 amp1 sig1  mu2 amp2 sig2 ...]
params = reshape(params, 3, n);
y = zeros(size(x));
for k = 1:n
    y = y + params(2,k) .* exp(-0.5 .* ((x - params(1,k)) ./ params(3,k)).^2);
end
end

function out = init_out()
out = struct();
out.knee         = NaN;
out.ap_fit       = [];
out.flat_spec    = [];
out.peak_params  = zeros(0,3);
out.gauss_params = zeros(0,3);
out.r_squared    = NaN;
out.error        = NaN;
out.flags        = struct('fitFailed', false, 'highError', false, 'noConverge', false);
end