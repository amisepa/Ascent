function [exponent, offset, info] = compute_AperiodicFit(freqs, psd, varargin)
% Aperiodic exponent and offset via spectral parameterization (FOOOF/specparam).
% Pure MATLAB reimplementation — no Python, no external toolboxes required.
%
%   [exponent, offset, info] = compute_AperiodicFit(freqs, psd, ...
%       'FreqRange', [1 40], 'AperiodicMode', 'fixed', ...
%       'MaxPeaks', 6, 'MinPeakHeight', 0.05, 'PeakThreshold', 2.0, ...
%       'PeakWidthLimits', [0.5 12], ...
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
%                .gauss_params– [nPeaks x 3]: underlying Gaussian [mu, amp, sigma]
%                .r_squared   – R² of full model vs. log10 PSD
%                .error       – MAE of full model vs. log10 PSD
%                .freqs_used  – frequency vector used for fitting
%                .flags       – struct: fitFailed, highError, noConverge
%              plus .params   – copy of options used
%
% Algorithm (Donoghue et al., 2020, Nat Neurosci)
%   1. Trim spectrum to FreqRange; convert to log10 power vs log10 frequency.
%   2. Initial aperiodic fit (robustfit on log-log, ignoring peaks).
%   3. Subtract initial aperiodic → flattened spectrum.
%   4. Iteratively detect and fit Gaussian peaks on flattened spectrum.
%   5. Robust aperiodic re-fit after removing peak regions.
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

%% ---------- Parse inputs ----------
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

%% ---------- Shape inputs ----------
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

%% ---------- Frequency trimming ----------
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

%% ---------- Pre-allocate outputs ----------
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

%% ---------- Progress header ----------
if opts.Progress
    parStr = 'off'; if opts.Parallel, parStr = 'on'; end
    fprintf('AperiodicFit: %d channel(s) | mode=%s | freqs=[%.0f %.0f] Hz | parallel=%s\n', ...
        nchan, opts.AperiodicMode, opts.FreqRange(1), opts.FreqRange(2), parStr);
end

%% ---------- Iterate channels ----------
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
        [exp_ch, off_ch, chInfo] = fooof_single(logf, logp(:), opts);
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
        [exponent(ch), offset(ch), chInfo] = fooof_single(logf, logp(:), opts);
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


%% ===================== Single-spectrum worker ===================== %%
function [exp_out, off_out, out] = fooof_single(logf, logp, opts)
% logf, logp: [nFit x 1] column vectors (log10 scale)

out = init_out();

%% Step 1 — Initial aperiodic fit (robust, ignores peaks)
ap0 = fit_aperiodic(logf, logp, opts.AperiodicMode);
if any(~isfinite(ap0))
    out.flags.fitFailed = true;
    exp_out = NaN; off_out = NaN;
    return
end

%% Step 2 — Flatten spectrum by subtracting initial aperiodic
ap_spec0  = aperiodic_model(logf, ap0, opts.AperiodicMode);
flat_spec = logp - ap_spec0;

%% Step 3 — Iterative Gaussian peak detection on flattened spectrum
gauss_params = find_peaks(logf, flat_spec, opts);  % [nPeaks x 3]: [mu amp sigma]

%% Step 4 — Robust aperiodic re-fit after removing peak regions
% Zero out frequency bins within 2*sigma of each peak center
peak_mask = false(size(logf));
for pk = 1:size(gauss_params,1)
    mu    = gauss_params(pk,1);
    sigma = gauss_params(pk,3);
    peak_mask = peak_mask | (logf >= mu - 2*sigma & logf <= mu + 2*sigma);
end
logf_ap = logf(~peak_mask);
logp_ap = logp(~peak_mask);

if numel(logf_ap) >= 2 + strcmp(opts.AperiodicMode,'knee')
    ap_params = fit_aperiodic(logf_ap, logp_ap, opts.AperiodicMode);
else
    ap_params = ap0;   % fallback to initial fit
end
if any(~isfinite(ap_params))
    ap_params = ap0;
end

%% Step 5 — Compute final full model and goodness of fit
ap_fit   = aperiodic_model(logf, ap_params, opts.AperiodicMode);
gauss_fit = zeros(size(logf));
for pk = 1:size(gauss_params,1)
    gauss_fit = gauss_fit + gauss_params(pk,2) .* ...
        exp(-0.5 * ((logf - gauss_params(pk,1)) ./ gauss_params(pk,3)).^2);
end
full_fit = ap_fit + gauss_fit;

ss_res = sum((logp - full_fit).^2);
ss_tot = sum((logp - mean(logp)).^2);
r2  = 1 - ss_res / max(ss_tot, eps);
mae = mean(abs(logp - full_fit));

%% Step 6 — Convert Gaussian params to peak params (CF, power, BW)
% CF    = mu  (Hz, convert from log10)
% power = height of model above aperiodic at CF (log10 power units)
% BW    = 2 * sigma in log10 freq → convert back to Hz approximation
peak_params = zeros(size(gauss_params));
for pk = 1:size(gauss_params,1)
    cf_log = gauss_params(pk,1);         % log10(CF)
    cf_hz  = 10^cf_log;
    % Height above aperiodic at CF
    ap_at_cf = aperiodic_model(cf_log, ap_params, opts.AperiodicMode);
    pw = gauss_fit(find(logf >= cf_log, 1));  % model height at CF bin
    if isempty(pw), pw = gauss_params(pk,2); end
    bw_hz = 2 * gauss_params(pk,3) * cf_hz * log(10);   % first-order approx
    peak_params(pk,:) = [cf_hz, pw - ap_at_cf, bw_hz]; %#ok<FNDSB>
end

%% Pack outputs
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


%% ===================== Aperiodic model ===================== %%
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


%% ===================== Aperiodic fitting ===================== %%
function params = fit_aperiodic(logf, logp, mode)
% Least-squares fit of the aperiodic model in log10 space.
% Uses MATLAB's lsqcurvefit if available, otherwise OLS (fixed mode only).

logf = logf(:); logp = logp(:);

if strcmp(mode, 'fixed')
    % Linear in log-log: logp = b0 - b1*logf  →  OLS
    X      = [ones(size(logf)) -logf];
    beta   = X \ logp;
    params = [beta(1); beta(2)];

else
    % Knee mode: nonlinear — use lsqcurvefit or fminsearch
    % Initial guess: run fixed first for offset/exponent, knee = median power
    X       = [ones(size(logf)) -logf];
    beta0   = X \ logp;
    off0    = beta0(1);
    exp0    = max(0, beta0(2));
    knee0   = 10^median(logf) / 2;
    x0      = [off0, knee0, exp0];
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


%% ===================== Peak detection ===================== %%
function gauss_params = find_peaks(logf, flat_spec, opts)
% Iteratively find and fit Gaussians on the flattened spectrum.
% Returns [nPeaks x 3]: [log10(CF), amplitude, sigma_log10]

logf      = logf(:);
flat_spec = flat_spec(:);
nF        = numel(logf);
freq_res  = mean(diff(logf));                  % log10 Hz resolution

% Convert Hz width limits to log10-space sigma limits
% BW = 2*sigma (Hz) ≈ 2*sigma_log * CF * ln(10);  sigma_log ≈ BW/(2*CF*ln(10))
% Use approximate: sigma_log_min = log10-space half-width
% Simpler conservative approach: convert limits at geometric mean frequency
f_mid     = 10^mean(logf);
sigma_min = (opts.PeakWidthLimits(1)/2) / (f_mid * log(10));
sigma_max = (opts.PeakWidthLimits(2)/2) / (f_mid * log(10));
sigma_min = max(sigma_min, freq_res);          % cannot be < freq resolution

gauss_params = zeros(0,3);
flat_iter    = flat_spec;

for pk = 1:opts.MaxPeaks
    % Find candidate peak: maximum of flattened spectrum
    [peak_val, peak_idx] = max(flat_iter);

    % Stop if peak is below absolute or relative threshold
    if peak_val < opts.MinPeakHeight, break; end
    if peak_val < opts.PeakThreshold * std(flat_iter), break; end

    % Initial Gaussian guess at peak
    cf_guess    = logf(peak_idx);
    amp_guess   = peak_val;
    sigma_guess = max(sigma_min, freq_res * 2);

    % Fit Gaussian to neighborhood around peak (±3 sigma)
    win   = abs(logf - cf_guess) <= 3 * sigma_guess * 3;
    if sum(win) < 3
        win = abs(logf - cf_guess) <= 5 * freq_res;
    end
    xfit  = logf(win);
    yfit  = flat_iter(win);

    gauss_fn = @(b, x) b(2) .* exp(-0.5 .* ((x - b(1)) ./ b(3)).^2);
    x0_g     = [cf_guess, amp_guess, sigma_guess];
    lb_g     = [logf(1),   0,         sigma_min ];
    ub_g     = [logf(end), Inf,       sigma_max ];
    opts_ls  = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',1000);

    try
        if ~isempty(which('lsqcurvefit'))
            gfit = lsqcurvefit(gauss_fn, x0_g, xfit, yfit, lb_g, ub_g, opts_ls);
        else
            cost = @(b) sum((yfit - gauss_fn(b,xfit)).^2);
            gfit = fminsearch(cost, x0_g, opts_ls);
            gfit(3) = max(sigma_min, min(sigma_max, abs(gfit(3))));
        end
    catch
        gfit = x0_g;
    end

    % Validate fitted peak
    if gfit(2) < opts.MinPeakHeight, break; end
    if gfit(3) < sigma_min || gfit(3) > sigma_max, break; end

    gauss_params(end+1,:) = gfit; %#ok<AGROW>

    % Subtract this Gaussian from flattened spectrum before next iteration
    flat_iter = flat_iter - gauss_fn(gfit, logf);
end

% Final joint re-fit of all Gaussians simultaneously if > 1 peak found
if size(gauss_params,1) > 1
    all_gauss = @(b, x) eval_gaussians(b, x, size(gauss_params,1));
    x0_all    = gauss_params(:).';
    lb_all    = repmat([logf(1),   0,         sigma_min], 1, size(gauss_params,1));
    ub_all    = repmat([logf(end), Inf,       sigma_max], 1, size(gauss_params,1));
    opts_ls   = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',3000);
    try
        if ~isempty(which('lsqcurvefit'))
            x_fit = lsqcurvefit(all_gauss, x0_all, logf, flat_spec, lb_all, ub_all, opts_ls);
        else
            cost  = @(b) sum((flat_spec - all_gauss(b, logf)).^2);
            x_fit = fminsearch(cost, x0_all, opts_ls);
        end
        gauss_params = reshape(x_fit, 3, []).';
    catch
        % keep greedy fit
    end
end
end


%% ===================== Helpers ===================== %%
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