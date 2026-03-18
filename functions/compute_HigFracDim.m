function [D, SD, info] = compute_HigFracDim(data, varargin)
% Higuchi Fractal Dimension for multichannel EEG data.
%
%   [D, SD, info] = compute_HiguchiFD(data, ...
%       'Kmax', 16, 'RobustFit', 'theilsen', ...
%       'Parallel', true, 'Progress', true, 'MinScales', 4)
%
% Inputs
%   data            : [n_channels x n_samples] (vector → 1 channel)
%                     Data are assumed to be already preprocessed (cleaned,
%                     artifact-free). No additional artifact rejection is
%                     applied.
%   Name-Value pairs:
%     'Kmax'       : maximum subsampling interval (default 16; Doyle et al.
%                    2017 recommend 16–36 for EEG; original Higuchi 1988
%                    used 8)
%     'RobustFit'  : 'theilsen' (default) | 'ols'
%                    Theil-Sen robust regression on the log(L(k)) vs
%                    log(1/k) scaling plot. Falls back to OLS when fewer
%                    than MinScales valid k values are available.
%     'Parallel'   : true|false, parfor over channels (default true)
%     'Progress'   : true|false, print per-channel results (default true)
%                    Parallel+Progress → text only (parfor-safe)
%     'MinScales'  : minimum number of valid k values required for robust
%                    fit (default 4); channels below this threshold return
%                    NaN with info.flags.fewScales = true
%
% Outputs
%   D   : Higuchi FD per channel [n_channels x 1], expected range [1, 2]
%         D ≈ 1 → smooth/regular signal; D ≈ 2 → highly complex/irregular
%   SD  : standard error of the slope estimate [n_channels x 1]
%   info: struct with per-channel diagnostics:
%           .Lk          – mean curve length at each k [Kmax x 1]
%           .logLk       – log(Lk)
%           .logInvk     – log(1/k), i.e. the x-axis of the scaling plot
%           .fitMethod   – 'theilsen' or 'ols'
%           .nUsedScales – number of valid k values used in fit
%           .flags       – struct: fewScales, highSE, outOfRange
%         plus .params   – copy of options used
%
% Algorithm
%   Implements Higuchi (1988) exactly. For each k = 1..Kmax and each
%   starting offset m = 1..k, the normalised mean curve length is:
%
%     L_m(k) = [ sum_{i=1}^{floor((N-m)/k)} |x(m+ik) - x(m+(i-1)k)| ]
%              * (N-1) / ( floor((N-m)/k) * k^2 )
%
%   L(k) = mean_{m=1}^{k} L_m(k)
%
%   HFD = slope of log(L(k)) vs log(1/k) via linear regression.
%
%   Note: normalisation uses floor((N-m)/k) per the original paper,
%   not round() as in some circulating implementations (e.g. Lord, UA).
%
% References
%   Higuchi, T. (1988). Physica D, 31, 277-283.
%   Accardo et al. (1997). Biol Cybern, 77, 339-350.
%   Esteller et al. (2001). IEEE Trans Biomed Eng, 48(2), 177-183.
%   Doyle et al. (2017). Clin Neurophysiol, 128(3), 438-443.
%     [validates Kmax range for EEG; recommends 16-36]
%   Sen, P.K. (1968). JASA, 63, 1379-1389.  [Theil-Sen regression]
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Ascent Plugin - Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

%% ---------- Parse inputs ----------
p = inputParser;
p.addRequired('data',       @(x) isnumeric(x) && ismatrix(x));
p.addParameter('Kmax',      16,         @(x) isnumeric(x) && isscalar(x) && x >= 2);
p.addParameter('RobustFit', 'theilsen', @(s) ischar(s) || isstring(s));
p.addParameter('Parallel',  true,       @(x) islogical(x) && isscalar(x));
p.addParameter('Progress',  true,       @(x) islogical(x) && isscalar(x));
p.addParameter('MinScales', 4,          @(x) isnumeric(x) && isscalar(x) && x >= 2);
p.parse(data, varargin{:});
opts = p.Results;

if strcmpi(opts.RobustFit, 'huber')
    warning('compute_HiguchiFD:HuberNotImplemented', ...
        'Huber fit not available. Using Theil-Sen instead.');
    opts.RobustFit = 'theilsen';
end

%% ---------- Shape ----------
if size(data,1) > size(data,2)
    data = data.';
end
[nchan, ~] = size(data);
D  = nan(nchan,1);
SD = nan(nchan,1);

%% ---------- Pre-allocate info ----------
info             = struct();
info.Lk          = cell(nchan,1);
info.logLk       = cell(nchan,1);
info.logInvk     = cell(nchan,1);
info.fitMethod   = cell(nchan,1);
info.nUsedScales = zeros(nchan,1);
info.flags       = cell(nchan,1);
info.params      = opts;

%% ---------- Progress header ----------
if opts.Progress
    parStr = 'off'; if opts.Parallel, parStr = 'on'; end
    fprintf('HiguchiFD: %d channel(s) | Kmax=%d | robust=%s | parallel=%s\n', ...
        nchan, opts.Kmax, lower(string(opts.RobustFit)), parStr);
end

%% ---------- Iterate channels ----------
if opts.Parallel && ~isempty(ver('parallel'))

    Lk_all     = cell(nchan,1);
    logLk_all  = cell(nchan,1);
    logIk_all  = cell(nchan,1);
    method_all = cell(nchan,1);
    nsc_all    = zeros(nchan,1);
    flags_all  = cell(nchan,1);

    parfor ch = 1:nchan
        [D(ch), SD(ch), chInfo] = higuchi_single(data(ch,:), opts);
        Lk_all{ch}     = chInfo.Lk;
        logLk_all{ch}  = chInfo.logLk;
        logIk_all{ch}  = chInfo.logInvk;
        method_all{ch} = chInfo.fitMethod;
        nsc_all(ch)    = chInfo.nUsedScales;
        flags_all{ch}  = chInfo.flags;
        if opts.Progress
            fprintf('  ch %3d/%3d: HFD=%7.5f  SE=%7.5f\n', ch, nchan, D(ch), SD(ch));
        end
    end

    info.Lk          = Lk_all;
    info.logLk       = logLk_all;
    info.logInvk     = logIk_all;
    info.fitMethod   = method_all;
    info.nUsedScales = nsc_all;
    info.flags       = flags_all;

else
    useWB = opts.Progress && usejava('desktop');
    hWB = [];
    if useWB
        try, hWB = waitbar(0, 'Computing Higuchi FD...', 'Name', 'compute_HiguchiFD');
        catch, hWB = []; end
    end

    for ch = 1:nchan
        [D(ch), SD(ch), chInfo] = higuchi_single(data(ch,:), opts);
        info.Lk{ch}           = chInfo.Lk;
        info.logLk{ch}        = chInfo.logLk;
        info.logInvk{ch}      = chInfo.logInvk;
        info.fitMethod{ch}    = chInfo.fitMethod;
        info.nUsedScales(ch)  = chInfo.nUsedScales;
        info.flags{ch}        = chInfo.flags;
        if opts.Progress
            fprintf('  ch %3d/%3d: HFD=%7.5f  SE=%7.5f\n', ch, nchan, D(ch), SD(ch));
            if ~isempty(hWB) && isvalid(hWB)
                try, waitbar(ch/nchan, hWB, sprintf('Higuchi FD... (%d/%d)', ch, nchan)); end
            end
        end
    end
    if ~isempty(hWB) && isvalid(hWB), try, close(hWB); end, end
end
end


%% ===================== Single-channel worker ===================== %%
function [dimension, standard_dev, out] = higuchi_single(sig, opts)

sig = sig(:).';
N   = numel(sig);

% Remove non-finite samples (should be rare in cleaned EEG)
if any(~isfinite(sig))
    sig = sig(isfinite(sig));
    N   = numel(sig);
end

% Guard: signal too short to support requested Kmax
if N < 2 * opts.Kmax
    dimension = NaN; standard_dev = NaN;
    out = make_empty_out(opts, 'tooShort');
    return
end

% Guard: constant signal
if (max(sig) - min(sig)) < eps
    dimension = NaN; standard_dev = NaN;
    out = make_empty_out(opts, 'constant');
    return
end

%% Core Higuchi algorithm (Higuchi 1988, equations 1-3)
Kmax = opts.Kmax;
Lk   = nan(Kmax, 1);

for k = 1:Kmax
    Lm = zeros(k,1);
    for m = 1:k
        nsteps = floor((N - m) / k);   % floor per Higuchi (1988), not round
        if nsteps < 1, continue; end
        idx        = m : k : m + nsteps*k;
        curve_len  = sum(abs(diff(sig(idx))));
        norm_factor = (N - 1) / (nsteps * k^2);
        Lm(m)      = curve_len * norm_factor;
    end
    valid_m = Lm > 0;
    if any(valid_m)
        Lk(k) = mean(Lm(valid_m));
    end
end

% Valid k values for regression
validK = find(isfinite(Lk) & Lk > 0);
if numel(validK) < 2
    dimension = NaN; standard_dev = NaN;
    out = make_empty_out(opts, 'fewScales');
    return
end

x2 = log(1 ./ validK(:));   % log(1/k)  - x-axis
y2 = log(Lk(validK));       % log(L(k)) - y-axis

%% Regression
forceOLS = strcmpi(opts.RobustFit, 'ols');
smallN   = numel(x2) < opts.MinScales || rank([ones(numel(x2),1) x2]) < 2;

if forceOLS || smallN
    X        = [ones(size(x2)) x2];
    beta     = X \ y2;
    e        = y2 - X*beta;
    s2       = (e.'*e) / max(1, numel(y2)-2);
    se       = sqrt(s2) .* sqrt(diag(pinv(X.'*X)));
    slope    = beta(2);
    slope_se = se(2);
    methodUsed = 'ols';
else
    % Standardize for numerical stability
    mx = mean(x2); sx = std(x2); if sx==0, sx=1; end
    my = mean(y2); sy = std(y2); if sy==0, sy=1; end
    xz = (x2 - mx) / sx;
    yz = (y2 - my) / sy;
    m_ts     = theil_sen_hig(xz, yz);
    b_ts     = median(yz - m_ts * xz);
    slope    = m_ts * sy / sx;
    % SE from OLS residuals on standardised scale
    Xz   = [ones(size(xz)) xz];
    ez   = yz - (b_ts + m_ts*xz);
    s2z  = (ez.'*ez) / max(1, numel(yz)-2);
    se_z = sqrt(s2z) .* sqrt(diag(pinv(Xz.'*Xz)));
    slope_se = se_z(2) * sy / sx;
    methodUsed = 'theilsen';
end

dimension    = slope;      % HFD = slope of log(L(k)) vs log(1/k); positive, in [1,2]
standard_dev = slope_se;

%% Pack output
out = struct();
out.Lk          = Lk;
out.logLk       = log(Lk);
out.logInvk     = log(1./(1:Kmax).');
out.fitMethod   = methodUsed;
out.nUsedScales = numel(validK);
out.flags = struct( ...
    'fewScales',  numel(validK) < opts.MinScales, ...
    'highSE',     ~isnan(standard_dev) && standard_dev > 0.05, ...
    'outOfRange', ~isnan(dimension)    && (dimension < 1 || dimension > 2));
end


%% ===================== Helpers ===================== %%
function out = make_empty_out(~, reason)
out = struct();
out.Lk          = [];
out.logLk       = [];
out.logInvk     = [];
out.fitMethod   = 'none';
out.nUsedScales = 0;
out.flags = struct( ...
    'fewScales',  any(strcmp(reason, {'fewScales','tooShort'})), ...
    'highSE',     false, ...
    'outOfRange', false);
end

function m = theil_sen_hig(x, y)
x = x(:); y = y(:);
N = numel(x);
if N < 2, m = 0; return; end
sl = zeros(N*(N-1)/2, 1);
t  = 0;
for i = 1:N-1
    dx = x(i+1:end) - x(i);
    dy = y(i+1:end) - y(i);
    ok = dx ~= 0 & isfinite(dx) & isfinite(dy);
    nn = sum(ok);
    if nn > 0
        sl(t+1:t+nn) = dy(ok) ./ dx(ok);
        t = t + nn;
    end
end
if t == 0
    m = (y(end)-y(1)) / max(eps, x(end)-x(1));
else
    m = median(sl(1:t));
end
end