function [RCmvMFE, scales] = compute_RCmvMFE(data, varargin)
% compute_RCmvMFE  Refined Composite Multivariate Multiscale Fuzzy Entropy (Azami & Escudero 2016)
% Exact, matches original RCmvMFE_mu and mvFE assumptions:
% 1) zscore per channel across time
% 2) r = sr * sum(std per channel)
% 3) coarse graining mean uses floor and plain mean (no omitnan)
% 4) mvFE uses the same N for m and m+1 normalizations (N = nsamp - max(M)*max(tau))
% 5) pairwise distances are exact (pdist if safe, else exact blocked summation), all in double
%
% Name value options:
%   'm'          : embedding dimension scalar sm (default 2)
%   'r'          : similarity bound scalar sr (default 0.15)
%   'tau'        : time lag scalar stau (default 1)
%   'n'          : fuzzy exponent (default 2)
%   'coarsing'   : 'mean' or 'var' (default 'mean')
%   'num_scales' : number of scales (default 15)
%   'Parallel'   : true or false, parallelize offsets within a scale (default false)
%   'Progress'   : true or false (default true)
%   'blockSize'  : block size for exact blocked fallback (default 2000)
%   'pdistMaxGB' : max GB allowed for pdist temporary vector (default 2.0)

% 1. Parse inputs
p = inputParser;
p.addRequired('data', @(x) (isstruct(x) && isfield(x,'data')) || isnumeric(x));
p.addParameter('m', 2, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('r', 0.15, @(x) isnumeric(x) && isscalar(x) && x>0 && x<2);
p.addParameter('tau', 1, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('n', 2, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('coarsing','mean', @(s) any(strcmpi(s,{'mean','var'})));
p.addParameter('num_scales', 15, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('parallel', false, @(x) islogical(x) && isscalar(x));
p.addParameter('progress', true, @(x) islogical(x) && isscalar(x));
p.addParameter('blockSize', 2000, @(x) isnumeric(x) && isscalar(x) && x>=100);
p.addParameter('pdistMaxGB', 2.0, @(x) isnumeric(x) && isscalar(x) && x>0);
p.parse(data, varargin{:});

sm          = p.Results.m;
sr          = p.Results.r;
stau        = p.Results.tau;
n_exp       = p.Results.n;
coarseType  = lower(p.Results.coarsing);
num_scales  = max(1, floor(p.Results.num_scales));
parallelMode= p.Results.parallel;
showProg    = p.Results.progress;
blockSize   = p.Results.blockSize;
pdistMaxGB  = p.Results.pdistMaxGB;

% 2. Coerce data
if isstruct(data)
    X = double(data.data);
else
    X = double(data);
end
if size(X,1) > size(X,2)
    X = X.';
end
[nChan, ~] = size(X);

% 3. Match original normalization (zscore per channel across time)
X = zscore(X, 0, 2);

% 4. Match original r scaling
r = sr * sum(std(X, 0, 2));

% 5. Embedding and tau vectors
M   = sm   * ones(1, nChan);
Tau = stau * ones(1, nChan);

% 6. Outputs
scales = 1:num_scales;
RCmvMFE = nan(1, num_scales);

% 7. Header
doPar = parallelMode && ~isempty(ver('parallel'));
if showProg
    parStr = 'off';
    if doPar, parStr = 'on'; end
    fprintf('RCmvMFE: channels=%d | m=%g, tau=%g, r=%g, n=%g | coarse=%s | scales=%d | parallel=%s\n', ...
        nChan, sm, stau, r, n_exp, coarseType, num_scales, parStr);
end

% Parallel only helps when there are enough offsets; otherwise overhead dominates
parMinOffsets = 6;

% 8. Compute: scales serial, offsets optionally parallel (no prints in workers)
for iScale = 1:num_scales

    if iScale == 1
        disp("NOTE: the 1st scale factor takes much longer than the rest of the scales.")
    end

    if showProg
        fprintf('  computing scale %d/%d...\n', iScale, num_scales);
        tScale = tic;
    end

    fi_m1_vec = zeros(1, iScale);
    fi_m2_vec = zeros(1, iScale);

    if doPar && iScale >= parMinOffsets
        parfor offset = 1:iScale
            Xs = coarsegrain_original(X(:, offset:end), iScale, coarseType);
            [~, fi_m1, fi_m2] = compute_mvFE_exact_originalN(Xs, M, r, n_exp, Tau, blockSize, pdistMaxGB);
            fi_m1_vec(offset) = fi_m1;
            fi_m2_vec(offset) = fi_m2;
        end
    else
        for offset = 1:iScale
            Xs = coarsegrain_original(X(:, offset:end), iScale, coarseType);
            [~, fi_m1, fi_m2] = compute_mvFE_exact_originalN(Xs, M, r, n_exp, Tau, blockSize, pdistMaxGB);
            fi_m1_vec(offset) = fi_m1;
            fi_m2_vec(offset) = fi_m2;
        end
    end

    AA = sum(fi_m1_vec);
    BB = sum(fi_m2_vec);

    if BB > 0 && isfinite(AA) && isfinite(BB)
        RCmvMFE(iScale) = log(AA / BB);
    else
        RCmvMFE(iScale) = NaN;
    end

    if showProg
        fprintf('     --> finished scale %d/%d (%.1fs)\n', iScale, num_scales, toc(tScale));
    end
end

end

% =====================================================================
function Y = coarsegrain_original(Data, S, method)
% Matches original Multi (mean) and Multi_var (var)

[nvar, L] = size(Data);
J = floor(L / S);
if J == 0
    Y = zeros(nvar, 0);
    return
end

Data = Data(:, 1:J*S);
Data = reshape(Data, nvar, S, J);   % nvar x S x J

switch lower(method)
    case 'mean'
        % original Multi: mean over each S block
        Y = squeeze(mean(Data, 2));         % nvar x J
    case 'var'
        % original Multi_var: var over each S block, MATLAB default normalization (N-1)
        Y = squeeze(var(Data, 0, 2));       % nvar x J
    otherwise
        error('method must be ''mean'' or ''var''.');
end
end

% =====================================================================
function [mvFE_Out, fi_m1, fi_m2] = compute_mvFE_exact_originalN(X, M, r, n, tau, blockSize, pdistMaxGB)
% Exact mvFE, respecting original assumption that N is the same in m and m+1 normalizations:
% N = nsamp - max(M)*max(tau)
% Uses exact pairwise Chebyshev distances, computed via pdist when safe, else exact blocked sum.

if iscolumn(M),   M   = M.';   end
if iscolumn(tau), tau = tau.'; end

[nvar, nsamp] = size(X);

max_M   = max(M);
max_tau = max(tau);
nn      = max_M * max_tau;

N = nsamp - nn;
if N < 2
    mvFE_Out = NaN; fi_m1 = NaN; fi_m2 = NaN;
    return
end

% m embedding with exactly N rows
A = embdN_fast(M, tau, X, N);

% fi_m1
fi_m1 = fuzzy_pairmean_cheby_exact(A, r, n, blockSize, pdistMaxGB);

% m+1 embeddings across subspaces, each with exactly N rows, then concatenate
M2 = repmat(M, nvar, 1) + eye(nvar);
dim2 = sum(M) + 1;

B = zeros(nvar*N, dim2, 'double');
row0 = 0;
for h = 1:nvar
    Bh = embdN_fast(M2(h,:), tau, X, N);
    B(row0 + (1:N), :) = Bh;
    row0 = row0 + N;
end

% fi_m2 (normalization uses nvar*N as in original)
fi_m2 = fuzzy_pairmean_cheby_exact(B, r, n, blockSize, pdistMaxGB);

mvFE_Out = log(fi_m1 / fi_m2);
end

% =====================================================================
function fi = fuzzy_pairmean_cheby_exact(V, r, n, blockSize, pdistMaxGB)
% fi = 2/(M(M-1)) * sum_{i<j} exp(-(d_ij^n)/r) exactly

M = size(V,1);
if M < 2
    fi = NaN;
    return
end

pairs = double(M) * double(M-1) / 2;
gbNeeded = (pairs * 8) / (1024^3); % pdist returns double

if gbNeeded <= pdistMaxGB
    d = pdist(V, 'chebychev');
    if n == 2
        s = sum(exp(-(d.*d)/r));
    else
        s = sum(exp(-(d.^n)/r));
    end
else
    s = pairwise_fuzzy_sum_cheby_blocked_double(V, r, n, blockSize);
end

fi = s * 2 / (M*(M-1));
end

% =====================================================================
function S = pairwise_fuzzy_sum_cheby_blocked_double(V, r, n, blockSize)
% Exact sum_{i<j} exp(-(d_ij^n)/r), Chebyshev distance, chunked, in double.

m = size(V,1);
p = size(V,2);

S = 0;

for i0 = 1:blockSize:m-1
    i1 = min(i0 + blockSize - 1, m-1);
    Vi = V(i0:i1,:);
    nb = size(Vi,1);

    for ii = 1:nb-1
        D = max(abs(Vi(ii+1:end,:) - Vi(ii,:)), [], 2);
        if n == 2
            S = S + sum(exp(-(D.*D)/r));
        else
            S = S + sum(exp(-(D.^n)/r));
        end
    end

    for j0 = i1+1:blockSize:m
        j1 = min(j0 + blockSize - 1, m);
        Vj = V(j0:j1,:);
        nj = size(Vj,1);

        Dmax = zeros(nb, nj);
        for k = 1:p
            Dmax = max(Dmax, abs(Vi(:,k) - Vj(:,k)'));
        end

        if n == 2
            S = S + sum(exp(-(Dmax.*Dmax)/r), 'all');
        else
            S = S + sum(exp(-(Dmax.^n)/r), 'all');
        end
    end
end
end

% =====================================================================
function A = embdN_fast(M, tau, ts, N)
% Fast, fixed-length multivariate delay embedding with exactly N rows.
% Matches original assumption that N is fixed and based on nn = max(M)*max(tau).

[nvar, nsamp] = size(ts);
dim = sum(M);

A = zeros(N, dim, 'double');

col0 = 0;
for j = 1:nvar
    mj = M(j);
    tj = tau(j);

    cols = col0 + (1:mj);
    col0 = col0 + mj;

    for k = 1:mj
        idx = (1:N) + (k-1)*tj;
        if idx(end) > nsamp
            error('embdN_fast:IndexOverflow', 'Not enough samples for requested embedding and N.');
        end
        A(:, cols(k)) = ts(j, idx).';
    end
end
end
