function [D, SD, info] = compute_fractal_volatility(data, varargin)
% compute_fractal_volatility  Fractal dimension via box counting with outlier-aware options.
%
%   [D, SD, info] = compute_fractal_volatility(data, 'Parallel', true, 'RejectBursts', true, ...)
%
%   Inputs
%     - data          : matrix [n_channels x n_samples]. If a vector, treated as 1 channel.
%     Name-Value pairs:
%       'RejectBursts' : true|false, sliding MAD z-score burst rejection (default true)
%       'WinFrac'      : window length for burst detection as fraction of [0,1] x-axis (default 0.02)
%       'ZThresh'      : MAD z threshold to flag bursts (default 6)
%       'RobustFit'    : 'huber' | 'theilsen' | 'ols' (default 'huber')
%       'ScaleTrimIQR' : true|false, trim scale-wise slope outliers by IQR rule (default true)
%       'MinBoxesPerCol': minimum vertical boxes to accept per column (default 1)
%       'MinJ'         : minimum dyadic level to consider, [] for auto (default [])
%       'MaxJ'         : maximum dyadic level to consider, [] for auto (default [])
%       'Parallel'     : true|false, enable parfor over channels (default false)
%       'Verbose'      : true|false, print per-channel info and progress (default true)
%
%   Outputs
%     - D   : fractal dimension per channel [n_channels x 1]
%     - SD  : standard error of slope per channel [n_channels x 1]
%     - info: struct with diagnostics
%         .r, .n, .validScales, .usedX, .usedY per channel in cell arrays
%         .cleanMask, .Jvec, .fitMethod, .rejectBursts, .params
%
%   Rationale and references
%     Short, high-power bursts can distort scale-based estimates even if they occupy a
%     tiny fraction of the recording. Pre-removal of such outliers and robust fitting of
%     scaling relationships stabilizes the slope estimate. Methodological guidance taken
%     from wavelet-leader multifractal outlier handling and robust estimation ideas:
%       Dumeur, M., Palva, J. M., & Ciuciu, P. (2025). Outlier detection and removal in
%       multifractal analysis of electrophysiological brain signals. EURASIP Journal on
%       Advances in Signal Processing, 2025(1), 35.
%
%   Copyright (C) Cedric Cannard 2025
%   Escape EEGLAB Plugin  https://github.com/amisepa/Escape
%
%   This code adapts general recommendations about detecting and downweighting outlier
%   segments before scale-exponent fitting. It does not implement the paper's full
%   p-leader segmentation, but mirrors its intent for box counting.

% 1. Parse inputs
p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ismatrix(x));
p.addParameter('RejectBursts', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('WinFrac',      0.02, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('ZThresh',      6,    @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('RobustFit',    'huber', @(s) ischar(s) || isstring(s));
p.addParameter('ScaleTrimIQR', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('MinBoxesPerCol', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('MinJ', [], @(x) isempty(x) || (isscalar(x) && x >= 1));
p.addParameter('MaxJ', [], @(x) isempty(x) || (isscalar(x) && x >= 1));
p.addParameter('Parallel', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Verbose',  true,  @(x) islogical(x) || isnumeric(x));
p.parse(data, varargin{:});
opts = p.Results;

% 2. Ensure [n_channels x n_samples]
if size(data,1) > size(data,2)
    data = data.';  % transpose to [n_channels x n_samples]
end
[nchan, ~] = size(data);
D  = nan(nchan,1);
SD = nan(nchan,1);

% Pre-alloc info containers
info = struct();
info.r           = cell(nchan,1);
info.n           = cell(nchan,1);
info.validScales = cell(nchan,1);
info.usedX       = cell(nchan,1);
info.usedY       = cell(nchan,1);
info.cleanMask   = cell(nchan,1);
info.Jvec        = cell(nchan,1);
info.fitMethod   = repmat({lower(char(opts.RobustFit))}, nchan, 1);
info.rejectBursts= repmat(logical(opts.RejectBursts), nchan, 1);
info.params      = opts;

if opts.Verbose
    fprintf('Fractal volatility via box counting on %d channels...\n', nchan);
end

% 3. Iterate channels, parallel or serial
if opts.Parallel
    parfor ch = 1:nchan %#ok<PFOUS>
        [D(ch), SD(ch), chInfo] = compute_fractal_volatility_single(data(ch,:), opts);
        info.r{ch}           = chInfo.r;
        info.n{ch}           = chInfo.n;
        info.validScales{ch} = chInfo.validScales;
        info.usedX{ch}       = chInfo.usedX;
        info.usedY{ch}       = chInfo.usedY;
        info.cleanMask{ch}   = chInfo.cleanMask;
        info.Jvec{ch}        = chInfo.Jvec;
        if opts.Verbose
            fprintf('Channel %3d/%3d  D=%6.4f  SE=%6.4f\n', ch, nchan, D(ch), SD(ch));
        end
    end
else
    progressbar('compute_fractal_volatility (serial)') % requires progressbar in path
    for ch = 1:nchan
        [D(ch), SD(ch), chInfo] = compute_fractal_volatility_single(data(ch,:), opts);
        info.r{ch}           = chInfo.r;
        info.n{ch}           = chInfo.n;
        info.validScales{ch} = chInfo.validScales;
        info.usedX{ch}       = chInfo.usedX;
        info.usedY{ch}       = chInfo.usedY;
        info.cleanMask{ch}   = chInfo.cleanMask;
        info.Jvec{ch}        = chInfo.Jvec;
        if opts.Verbose
            fprintf('Channel %3d/%3d  D=%6.4f  SE=%6.4f\n', ch, nchan, D(ch), SD(ch));
        end
        progressbar(ch / nchan);
    end
end
end

%% Single-channel worker
function [dimension, standard_dev, out] = compute_fractal_volatility_single(sig, opts)
sig = sig(:).';           % row
N   = numel(sig);

% 1. Build [x,y] embedded in unit square
x = (1:N).';
y = sig(:);
x = (x - min(x)) / max(eps, (max(x)-min(x)));
y = (y - min(y)) / max(eps, (max(y)-min(y)));

% 2. Optional burst rejection using sliding MAD z-score on y
maskGood = true(N,1);
if opts.RejectBursts
    win = max(5, round(opts.WinFrac * N));
    y_med = movmedian(y, win, 'Endpoints','shrink');
    y_mad = movmad(y, win, 1, 'Endpoints','shrink');
    mad_eps = median(y_mad(y_mad>0));
    if isempty(mad_eps) || mad_eps==0, mad_eps = 1e-6; end
    z = abs(y - y_med) ./ max(y_mad, mad_eps);
    maskGood = z <= opts.ZThresh;
end

xg = x(maskGood); yg = y(maskGood);
[xg, idx] = sort(xg, 'ascend'); yg = yg(idx);
if numel(xg) < 10
    % fall back to raw if cleaning removed too much
    xg = x; yg = y; maskGood = true(N,1);
end

% 3. Determine dyadic scales based on min spacing of x
dx = diff(xg); dx(dx<=0) = [];
if isempty(dx), dx = 1/N; end
minwidth = floor(abs(log2(max(eps, min(dx))))) - 1;
minwidth = max(minwidth, 2);
Jmax_auto = minwidth; Jmin_auto = 1;
Jmin = Jmin_auto; if ~isempty(opts.MinJ), Jmin = max(1, opts.MinJ); end
Jmax = Jmax_auto; if ~isempty(opts.MaxJ), Jmax = min(opts.MaxJ, Jmax_auto); end
Jvec = Jmin:Jmax;

% 4. Box counting over columns for each scale
n = zeros(numel(Jvec),1);
r = zeros(numel(Jvec),1);
for ii = 1:numel(Jvec)
    j = Jvec(ii);
    width = 2^(-j);
    n(ii) = count_boxes_in_columns(xg, yg, width, opts.MinBoxesPerCol);
    r(ii) = width;
end

% 5. Assemble log arrays and drop degenerate scales
valid = n > 0 & r > 0;
x2 = log(r(valid));
y2 = log(n(valid));

% Optional trimming of scale-wise slopes by IQR
if opts.ScaleTrimIQR && numel(x2) >= 5
    s = -gradient(log(n(valid))) ./ gradient(log(r(valid)));
    IQRs = iqr(s);
    keep = abs(s - median(s)) <= IQRs/2;
    if numel(keep) == numel(x2)
        x2 = x2(keep); y2 = y2(keep);
    end
end

% 6. Fit slope robustly
fitMethod = lower(char(opts.RobustFit));
switch fitMethod
    case 'theilsen'
        m = theil_sen(x2, y2);
        b = median(y2 - m*x2);
        yhat = m*x2 + b;
        X = [ones(size(x2)) x2];
        e = y2 - yhat;
        sigma2 = (e'*e) / max(1, numel(y2)-2);
        C = pinv(X.'*X);
        se = sqrt(sigma2) * sqrt(diag(C));
        slope = m; slope_se = se(2);
    case 'huber'
        try
            [bvec, stats] = robustfit(x2, y2, 'huber');
            slope = bvec(2); slope_se = stats.se(2);
        catch
            % fallback to OLS if robustfit not available
            X = [ones(size(x2)) x2];
            beta = X \ y2;
            yhat = X*beta;
            e = y2 - yhat;
            C = pinv(X.'*X);
            sigma2 = (e.'*e) / max(1, numel(y2)-2);
            se = sqrt(sigma2) .* sqrt(diag(C));
            slope = beta(2); slope_se = se(2);
        end
    otherwise
        X = [ones(size(x2)) x2];
        beta = X \ y2;
        yhat = X*beta;
        e = y2 - yhat;
        C = pinv(X.'*X);
        sigma2 = (e.'*e) / max(1, numel(y2)-2);
        se = sqrt(sigma2) .* sqrt(diag(C));
        slope = beta(2); slope_se = se(2);
end

dimension    = -slope;
standard_dev = slope_se;

% 7. Pack diagnostics
out = struct();
out.r = r;
out.n = n;
out.validScales = valid;
out.usedX = x2;
out.usedY = y2;
out.cleanMask = maskGood;
out.Jvec = Jvec;

end

%% Helpers
function boxcount = count_boxes_in_columns(xg, yg, width, minBoxes)
% Count vertical boxes needed inside each x column of width.
nCols = 2^round(-log2(width));
col = min(floor(xg / width) + 1, nCols);
boxcount = 0;
for c = 1:nCols
    idx = (col == c);
    if ~any(idx), continue; end
    yc = yg(idx);
    if numel(yc) == 1
        boxcount = boxcount + 1;
    else
        rawcount = (max(yc) - min(yc)) / width;
        rawcount = rawcount + rem(min(yc), width);
        cnt = ceil(rawcount);
        if cnt >= minBoxes
            boxcount = boxcount + cnt;
        end
    end
end
end

function m = theil_sen(x, y)
x = x(:); y = y(:);
N = numel(x);
if N < 2
    m = 0; return
end
K = N*(N-1)/2;
sl = zeros(K,1);
t = 0;
for i = 1:N-1
    dx = x(i+1:end) - x(i);
    dy = y(i+1:end) - y(i);
    k = dx ~= 0 & isfinite(dx) & isfinite(dy);
    nn = sum(k);
    if nn > 0
        sl(t+1:t+nn) = dy(k) ./ dx(k);
        t = t + nn;
    end
end
if t == 0
    m = (y(end)-y(1)) / max(eps, (x(end)-x(1)));
else
    m = median(sl(1:t));
end
end
