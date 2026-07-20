function SampEn = compute_SampEn(data, varargin)
% compute_SampEn  Computes Sample Entropy (SampEn) across multichannel data.
%
%   SampEn = compute_SampEn(data, 'm', 2, 'tau', 1, 'r', .15, ...
%       'ZScore', true, 'BlockSize', 2000, 'pdistMaxGB', 2.0, ...
%       'Parallel', true, 'Progress', true)
%
% Inputs:
%   data        : EEG data matrix [n_channels x n_samples]
%   'm'         : embedding dimension (default = 2)
%   'tau'       : time lag (default = 1)
%   'r'         : similarity bound (default = 0.15)
%   'ZScore'    : z-score each channel before computing SampEn (default = true).
%                 Set false when the caller has already normalised the data and
%                 'r' is an absolute Chebyshev threshold (used by compute_MSE).
%   'BlockSize' : block size for exact blocked fallback (default = 2000)
%   'pdistMaxGB': max GB allowed for pdist temporary vector (default = 2.0)
%   'Parallel'  : logical true/false to enable parfor over channels (default = true)
%   'Progress'  : logical true/false to show progress
%
% Output:
%   SampEn      : [n_channels x 1] Sample Entropy per channel
%
% Notes:
%   • When ZScore is true, data is z-scored per channel across time and r is a
%     fraction of the (unit) SD. When false, data is used as-is and r is the
%     absolute Chebyshev tolerance.
%   • Pairwise Chebyshev match fractions are computed exactly using pdist
%     when memory allows, otherwise an exact blocked fallback is used.

% ---------------- Parse inputs ----------------
p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ndims(x) == 2);
p.addParameter('m', 2,           @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('tau', 1,         @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('r', 0.15,        @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('ZScore', true,   @(x) islogical(x) && isscalar(x));
p.addParameter('BlockSize', 2000,@(x) isnumeric(x) && isscalar(x) && x >= 100);
p.addParameter('pdistMaxGB', 2.0,@(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Parallel', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true, @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
tau          = p.Results.tau;
r            = p.Results.r;
doZScore     = p.Results.ZScore;
blockSize    = p.Results.BlockSize;
pdistMaxGB   = p.Results.pdistMaxGB;
parallelMode = p.Results.Parallel;
showProgress = p.Results.Progress;

if size(data,1) > size(data,2)
    data = data.';
end
[nchan, ~] = size(data);

% ---------------- Z-score per channel (optional) ----------------
data_z = data;
if doZScore
    for c = 1:nchan
        x  = data(c,:);
        mu = mean(x, 'omitnan');
        sd = std(x, 0, 'omitnan');
        if ~isfinite(sd) || sd == 0
            data_z(c,:) = 0;
        else
            data_z(c,:) = (x - mu) ./ sd;
        end
    end
end

SampEn = nan(nchan, 1);

% ---------------- Progress header ----------------
if showProgress
    if parallelMode && ~isempty(ver('parallel'))
        fprintf('SampEn: %d channel(s) | m=%g, tau=%g, r=%g | parallel=on\n', ...
            nchan, m, tau, r);
    else
        fprintf('SampEn: %d channel(s) | m=%g, tau=%g, r=%g | parallel=off\n', ...
            nchan, m, tau, r);
    end
end

useWB = ~parallelMode && usejava('desktop') && showProgress;
hWB = [];
if useWB
    try
        hWB = waitbar(0, 'Computing Sample Entropy...', 'Name', 'compute_SampEn');
    catch
        hWB = [];
    end
end

useDQ = parallelMode && ~isempty(ver('parallel')) && showProgress;
if useDQ
    dq = parallel.pool.DataQueue;
    nDone = 0;
    afterEach(dq, @notifyProgress);
end

% ---------------- Compute per channel ----------------
if parallelMode && ~isempty(ver('parallel'))
    parfor iChan = 1:nchan
        SampEn(iChan) = compute_SampEn_single(data_z(iChan,:), m, r, tau, blockSize, pdistMaxGB);
        if useDQ
            send(dq, 1);
        end
    end
else
    for iChan = 1:nchan
        SampEn(iChan) = compute_SampEn_single(data_z(iChan,:), m, r, tau, blockSize, pdistMaxGB);

        if ~isempty(hWB) && isvalid(hWB)
            try
                waitbar(iChan/nchan, hWB, sprintf('Computing Sample Entropy... (%d/%d)', iChan, nchan));
            catch
            end
        end
    end

    if ~isempty(hWB) && isvalid(hWB)
        try
            close(hWB);
        catch
        end
    end
end

    function notifyProgress(~)
        nDone = nDone + 1;
        step = max(1, round(0.05 * nchan));
        if nDone == 1 || nDone == nchan || mod(nDone, step) == 0
            fprintf('  progress: ch %d/%d\n', nDone, nchan);
        end
    end
end

% =========================================================================
function SampEn = compute_SampEn_single(signal, m, r, tau, blockSize, pdistMaxGB)

signal = signal(isfinite(signal));
N = numel(signal);

if N <= m + 1
    SampEn = NaN;
    return
end

Xm  = embed_tau(signal, m,   tau);
Xm1 = embed_tau(signal, m+1, tau);

if size(Xm,1) < 2 || size(Xm1,1) < 2
    SampEn = NaN;
    return
end

Bm  = match_fraction_cheby_exact(Xm,  r, blockSize, pdistMaxGB);
Am1 = match_fraction_cheby_exact(Xm1, r, blockSize, pdistMaxGB);

if Bm > 0 && Am1 > 0 && isfinite(Bm) && isfinite(Am1)
    SampEn = -log(Am1 / Bm);
else
    SampEn = NaN;
end
end

% =========================================================================
function X = embed_tau(signal, m, tau)
signal = signal(:).';
N = numel(signal);
L = N - (m - 1) * tau;
if L <= 0
    X = zeros(0, m);
    return
end
idx  = (0:(m-1)) * tau;
rows = (1:L).';
X = signal(rows + idx);
end

% =========================================================================
function p = match_fraction_cheby_exact(V, r, blockSize, pdistMaxGB)
M = size(V,1);
if M < 2
    p = NaN;
    return
end

pairs = double(M) * double(M - 1) / 2;
gbNeeded = (pairs * 8) / (1024^3);

if gbNeeded <= pdistMaxGB
    d = pdist(V, 'chebychev');
    p = mean(d <= r);
else
    hits = pairwise_match_count_cheby_blocked(V, r, blockSize);
    p = hits / pairs;
end
end

% =========================================================================
function hits = pairwise_match_count_cheby_blocked(V, r, blockSize)

m = size(V,1);
pdim = size(V,2);
hits = 0;

for i0 = 1:blockSize:m-1
    i1 = min(i0 + blockSize - 1, m - 1);
    Vi = V(i0:i1,:);
    nb = size(Vi,1);

    % within-block
    for ii = 1:nb-1
        D = max(abs(Vi(ii+1:end,:) - Vi(ii,:)), [], 2);
        hits = hits + nnz(D <= r);
    end

    % across-blocks
    for j0 = i1+1:blockSize:m
        j1 = min(j0 + blockSize - 1, m);
        Vj = V(j0:j1,:);
        nj = size(Vj,1);

        Dmax = zeros(nb, nj);
        for k = 1:pdim
            Dmax = max(Dmax, abs(Vi(:,k) - Vj(:,k)'));
        end

        hits = hits + nnz(Dmax <= r);
    end
end
end