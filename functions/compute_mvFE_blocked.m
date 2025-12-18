function [mvFE_Out, fi_m1, fi_m2] = compute_mvFE_blocked(X, M, r, n, tau, blockSize)
% Memory-safe mvFE with blocked pairwise sums (Chebyshev)
% Handles m and m+1 potentially having different N.

if nargin < 7 || isempty(blockSize)
    blockSize = 2000;
end

[nvar, nsamp] = size(X);

% 1) Dimension m
max_M   = max(M);
max_tau = max(tau);
nn1     = max_M * max_tau;
N1      = nsamp - nn1;

if N1 < 2
    mvFE_Out = NaN; fi_m1 = NaN; fi_m2 = NaN;
    return
end

A = embd(M, tau, X);                 % should be N1 x dim
A = single(A);
N1 = size(A,1);                      % trust embd output

S1 = pairwise_fuzzy_sum_cheby(A, r, n, blockSize);
fi_m1 = S1 * 2 / (N1*(N1-1));

% 2) Dimension m+1 across nvar subspaces
M2 = repmat(M, nvar, 1) + eye(nvar);

Bcell = cell(nvar,1);
N2 = inf;

for h = 1:nvar
    Bh = embd(M2(h,:), tau, X);      % may be (N1 or N1-1) x dim2 depending on horizon
    Bcell{h} = single(Bh);
    N2 = min(N2, size(Bh,1));
end

if N2 < 2
    mvFE_Out = NaN; fi_m2 = NaN;
    return
end

% Truncate all to common N2, then concatenate
for h = 1:nvar
    Bcell{h} = Bcell{h}(1:N2, :);
end
B = vertcat(Bcell{:});               % (nvar*N2) x dim2

NN = nvar * N2;

S2 = pairwise_fuzzy_sum_cheby(B, r, n, blockSize);
fi_m2 = S2 * 2 / (NN*(NN-1));

mvFE_Out = log(fi_m1 / fi_m2);
end


function S = pairwise_fuzzy_sum_cheby(V, r, n, blockSize)
% S = sum_{i<j} exp(-(d_ij^n)/r), d_ij Chebyshev distance

m = size(V,1);
p = size(V,2);
S = 0;

for i0 = 1:blockSize:m-1
    i1 = min(i0 + blockSize - 1, m-1);
    Vi = V(i0:i1,:);
    nb = size(Vi,1);

    % within-block upper triangle
    for ii = 1:nb-1
        D = max(abs(Vi(ii+1:end,:) - Vi(ii,:)), [], 2);
        S = S + sum(exp(-(D.^n)/r), 'omitnan');
    end

    % cross blocks
    for j0 = i1+1:blockSize:m
        j1 = min(j0 + blockSize - 1, m);
        Vj = V(j0:j1,:);
        nj = size(Vj,1);

        Dmax = zeros(nb, nj, 'single');
        for k = 1:p
            Dmax = max(Dmax, abs(Vi(:,k) - Vj(:,k)'));
        end

        S = S + sum(exp(-(Dmax.^n)/r), 'all', 'omitnan');
    end
end
end
