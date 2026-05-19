function [RCmvMFE, scales] = compute_RCmvMFE(data, varargin)
% compute_RCmvMFE  Refined Composite Multivariate Multiscale Fuzzy Entropy.
%
%   [RCmvMFE, scales] = compute_RCmvMFE(data, 'm', 2, 'r', 0.15, 'tau', 1, ...
%       'n', 2, 'coarsing', 'mean', 'num_scales', 15, ...
%       'Parallel', false, 'Progress', true, ...
%       'blockSize', 2000, 'pdistMaxGB', 2.0)
%
% Inputs:
%   data        : multichannel data [n_channels x n_samples] or EEGLAB struct
%   'm'         : scalar embedding dimension (default = 2)
%   'r'         : scalar threshold scaling factor sr (default = 0.15)
%   'tau'       : scalar time lag (default = 1)
%   'n'         : fuzzy exponent (default = 2)
%   'coarsing'  : 'mean' 'median' 'trimmed mean' 'std' or 'var' (default = 'std')
%   'num_scales': number of scales (default = 15)
%   'Parallel'  : parallelize offsets within a scale when useful (default = false)
%   'Progress'  : show progress (default = true)
%   'blockSize' : block size for exact blocked fallback (default = 2000)
%   'pdistMaxGB': max GB allowed for pdist temporary vector (default = 2.0)
%
% Outputs:
%   RCmvMFE     : [1 x num_scales] refined composite multivariate multiscale fuzzy entropy
%   scales      : 1:num_scales
%
% Notes:
%   • Channels are z-score normalized across time.
%   • The similarity bound is scaled as r = sr * sum(std(X_k)).
%   • Coarse-graining uses the external coarsegrain() helper.
%   • Pairwise similarities are computed exactly via fuzzy_pairmean_cheby_exact(),
%     which uses pdist when memory allows and an exact blocked fallback otherwise.
%   • This implementation follows the original RCmvMFE structure where, at each scale,
%     offset-specific phi_m and phi_m+1 terms are summed across offsets and the entropy
%     is computed as log(sum(phi_m) / sum(phi_m+1)).

% ---------- Parse inputs ----------
p = inputParser;
p.addRequired('data', @(x) (isstruct(x) && isfield(x,'data')) || isnumeric(x));
p.addParameter('m', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('r', 0.15, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 2);
p.addParameter('tau', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('n', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('coarsing', 'mean', @(s) any(strcmpi(s, {'mean','median','trimmed mean','trimmed','trimmean','std','sd','standard deviation','var','variance'})));
p.addParameter('num_scales', 15, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('Parallel', false, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true, @(x) islogical(x) && isscalar(x));
p.addParameter('blockSize', 2000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
p.addParameter('pdistMaxGB', 2.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.parse(data, varargin{:});

sm           = p.Results.m;
sr           = p.Results.r;
stau         = p.Results.tau;
n_exp        = p.Results.n;
coarseType   = lower(p.Results.coarsing);
num_scales   = max(1, floor(p.Results.num_scales));
parallelMode = p.Results.Parallel;
showProg     = p.Results.Progress;
blockSize    = p.Results.blockSize;
pdistMaxGB   = p.Results.pdistMaxGB;

if num_scales < 5
    error('compute_RCmvMFE:TooFewScales', ...
        'Multiscale entropy cannot be computed with less than 5 scales.');
end

% ---------- Coerce data ----------
if isstruct(data)
    X = double(data.data);
else
    X = double(data);
end
if size(X,1) > size(X,2)
    X = X.';
end
[nChan, ~] = size(X);

% ---------- Normalize per channel ----------
X = zscore(X, 0, 2);

% ---------- Original r scaling ----------
r = sr * sum(std(X, 0, 2));

% ---------- Embedding vectors ----------
M   = sm   * ones(1, nChan);
Tau = stau * ones(1, nChan);

% ---------- Outputs ----------
% Exclude scale 1 for std/var operators (std of single-sample bin = 0)
ct = lower(strtrim(coarseType));
stdvarOp = any(strcmp(ct, {'std','sd','standard deviation','var','variance'}));
if stdvarOp
    scales   = 2:num_scales;
    scaleVec = 2:num_scales;
else
    scales   = 1:num_scales;
    scaleVec = 1:num_scales;
end
RCmvMFE = nan(1, numel(scales));

% ---------- Header ----------
doPar = parallelMode && ~isempty(ver('parallel'));
if showProg
    parStr = 'off';
    if doPar
        parStr = 'on';
    end
    fprintf('RCmvMFE: channels=%d | m=%g, tau=%g, r=%g, n=%g | coarse=%s | scales=%s | parallel=%s\n', ...
        nChan, sm, stau, r, n_exp, coarseType, mat2str(scales), parStr);
end

parMinOffsets = 6;

if showProg
    progressbar('Computing RCmvMFE for each scale factor')
    disp('NOTE: the lower scale factors take much longer than the higher ones.')
end

% ---------- Main loop ----------
for si = 1:numel(scaleVec)
    iScale = scaleVec(si);
    if showProg
        fprintf('  computing scale %d/%d...\n', iScale, num_scales);
        tScale = tic;
    end

    fi_m1_vec = zeros(1, iScale);
    fi_m2_vec = zeros(1, iScale);

    % Precompute coarse-grained series on client side
    XsCell = cell(iScale, 1);
    for offset = 1:iScale
        xoff = X(:, offset:end);

        % Match original coarse-graining assumptions:
        % truncate to full blocks before calling column-wise coarsegrain helper
        [nvar, L] = size(xoff);
        J = floor(L / iScale);

        if J == 0
            XsCell{offset} = zeros(nvar, 0);
            continue
        end

        xoff = xoff(:, 1:J*iScale);
        xoff = reshape(xoff, nvar, iScale, J);   % nvar x scale x J

        Xs = zeros(nvar, J);
        for ch = 1:nvar
            Xs(ch, :) = coarsegrain(squeeze(xoff(ch, :, :)), coarseType);
        end
        XsCell{offset} = Xs;
    end

    if doPar && iScale >= parMinOffsets
        parfor offset = 1:iScale
            Xs = XsCell{offset};
            [~, fi_m1, fi_m2] = compute_mvFE_local(Xs, M, r, n_exp, Tau, blockSize, pdistMaxGB);
            fi_m1_vec(offset) = fi_m1;
            fi_m2_vec(offset) = fi_m2;
        end
    else
        for offset = 1:iScale
            Xs = XsCell{offset};
            [~, fi_m1, fi_m2] = compute_mvFE_local(Xs, M, r, n_exp, Tau, blockSize, pdistMaxGB);
            fi_m1_vec(offset) = fi_m1;
            fi_m2_vec(offset) = fi_m2;
        end
    end

    AA = sum(fi_m1_vec);
    BB = sum(fi_m2_vec);

    if BB > 0 && isfinite(AA) && isfinite(BB)
        RCmvMFE(si) = log(AA / BB);
    else
        RCmvMFE(si) = NaN;
    end

    if showProg
        fprintf('     --> finished scale %d/%d (%.1fs)\n', iScale, num_scales, toc(tScale));
        progressbar(si / numel(scaleVec))
    end
end
end

% =========================================================================
function [mvFE_Out, fi_m1, fi_m2] = compute_mvFE_local(X, M, r, n, tau, blockSize, pdistMaxGB)
% Local multivariate fuzzy entropy core used inside RCmvMFE.
% Matches original assumption that N is the same in m and m+1 normalizations:
%   N = nsamp - max(M)*max(tau)

if iscolumn(M),   M   = M.';   end
if iscolumn(tau), tau = tau.'; end

[nChan, nSamp] = size(X);

max_M   = max(M);
max_tau = max(tau);
nn      = max_M * max_tau;

N = nSamp - nn;
if N < 2
    mvFE_Out = NaN;
    fi_m1 = NaN;
    fi_m2 = NaN;
    return
end

% m embedding
A = embdN_fast(M, tau, X, N);
fi_m1 = fuzzy_pairmean_cheby_exact(A, r, n, blockSize, pdistMaxGB);

% m+1 pooled construction across subspaces
M2 = repmat(M, nChan, 1) + eye(nChan);
dim2 = sum(M) + 1;

B = zeros(nChan * N, dim2, 'double');
row0 = 0;
for h = 1:nChan
    Bh = embdN_fast(M2(h,:), tau, X, N);
    B(row0 + (1:N), :) = Bh;
    row0 = row0 + N;
end

fi_m2 = fuzzy_pairmean_cheby_exact(B, r, n, blockSize, pdistMaxGB);

mvFE_Out = log(fi_m1 / fi_m2);
end

% =========================================================================
function A = embdN_fast(M, tau, ts, N)
% Fixed-length multivariate delay embedding with exactly N rows.

dim = sum(M);
A = zeros(N, dim, 'double');

col0 = 0;
for j = 1:numel(M)
    mj = M(j);
    tj = tau(j);

    cols = col0 + (1:mj);
    col0 = col0 + mj;

    for k = 1:mj
        idx = (1:N) + (k - 1) * tj;
        if idx(end) > size(ts, 2)
            error('embdN_fast:IndexOverflow', ...
                'Not enough samples for requested embedding and N.');
        end
        A(:, cols(k)) = ts(j, idx).';
    end
end
end