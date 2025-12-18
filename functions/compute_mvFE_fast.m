function [mvFE_Out, fi_m1, fi_m2] = compute_mvFE_fast(X, M, r, n, tau, varargin)
%
% This function calculates multivariate fuzzy entropy (mvFE) of a multivariate signal
% Accelerated, memory safe version (no pdist).
%
% Inputs:
%   X   : multivariate signal [nvar x nsamp]
%   M   : embedding vector [1 x nvar] or [nvar x 1]
%   r   : scalar threshold
%   n   : fuzzy power (typically 2)
%   tau : time lag vector [1 x nvar] or [nvar x 1]
%
% Name value options:
%   'blockSize' : block size for pairwise accumulation (default 2000)
%   'useSingle' : cast embedded matrices to single (default true)
%   'progress'  : print lightweight progress (default false)
%
% Outputs:
%   mvFE_Out : scalar mvFE of X
%   fi_m1    : scalar global quantity in dimension m
%   fi_m2    : scalar global quantity in dimension m+1
%
% Ref:
%   H. Azami and J. Escudero, Physica A, 2016.
%

% 1. Parse inputs
p = inputParser;
p.addParameter('blockSize', 2000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
p.addParameter('useSingle', true, @(x) islogical(x) && isscalar(x));
p.addParameter('progress',  false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});

blockSize = p.Results.blockSize;
useSingle = p.Results.useSingle;
showProg  = p.Results.progress;

% 2. Basic checks and shapes
if iscolumn(M),   M   = M.';   end
if iscolumn(tau), tau = tau.'; end

[nvar, nsamp] = size(X);

max_M   = max(M);
max_tau = max(tau);

% 3. Dimension m embedding
nn1 = max_M * max_tau;
N1_est = nsamp - nn1;

if N1_est < 2
    mvFE_Out = NaN; fi_m1 = NaN; fi_m2 = NaN;
    return
end

A = embd(M, tau, X);
N1 = size(A,1);
if N1 < 2
    mvFE_Out = NaN; fi_m1 = NaN; fi_m2 = NaN;
    return
end

if useSingle, A = single(A); end

if showProg
    fprintf('compute_mvFE_fast: m stage, A = %d x %d\n', size(A,1), size(A,2));
end

S1 = pairwise_fuzzy_sum_cheby(A, r, n, blockSize, showProg);
fi_m1 = S1 * 2 / (N1*(N1-1));

% 4. Dimension m+1 embedding across nvar subspaces
M2 = repmat(M, nvar, 1) + eye(nvar);

Bcell = cell(nvar,1);
N2 = inf;

for h = 1:nvar
    Bh = embd(M2(h,:), tau, X);
    if useSingle, Bh = single(Bh); end
    Bcell{h} = Bh;
    N2 = min(N2, size(Bh,1));
end

if N2 < 2
    mvFE_Out = NaN; fi_m2 = NaN;
    return
end

% 5. Truncate to common length to avoid off by 1 mismatches
for h = 1:nvar
    Bcell{h} = Bcell{h}(1:N2, :);
end
B = vertcat(Bcell{:});   % (nvar*N2) x dim2

NN = nvar * N2;

if showProg
    fprintf('compute_mvFE_fast: m+1 stage, B = %d x %d\n', size(B,1), size(B,2));
end

S2 = pairwise_fuzzy_sum_cheby(B, r, n, blockSize, showProg);
fi_m2 = S2 * 2 / (NN*(NN-1));

% 6. Final mvFE
mvFE_Out = log(fi_m1 / fi_m2);

end


% ========================================================================
function S = pairwise_fuzzy_sum_cheby(V, r, n, blockSize, showProg)
% Sum over i<j of exp(-(d_ij^n)/r), Chebyshev distance
%
% V: [m x p]

m = size(V,1);
p = size(V,2);

S = 0;

% 1. Blocked accumulation
for i0 = 1:blockSize:m-1
    i1 = min(i0 + blockSize - 1, m-1);
    Vi = V(i0:i1,:);
    nb = size(Vi,1);

    if showProg && i0 == 1
        fprintf('  pairwise blocks: %d rows, blockSize=%d\n', m, blockSize);
    end

    % 2. Within block upper triangle
    for ii = 1:nb-1
        D = max(abs(Vi(ii+1:end,:) - Vi(ii,:)), [], 2);
        S = S + sum(exp(-(D.^n)/r), 'omitnan');
    end

    % 3. Cross blocks
    for j0 = i1+1:blockSize:m
        j1 = min(j0 + blockSize - 1, m);
        Vj = V(j0:j1,:);
        nj = size(Vj,1);

        Dmax = zeros(nb, nj, 'like', V);
        for k = 1:p
            Dmax = max(Dmax, abs(Vi(:,k) - Vj(:,k)'));
        end

        S = S + sum(exp(-(Dmax.^n)/r), 'all', 'omitnan');
    end
end

end
