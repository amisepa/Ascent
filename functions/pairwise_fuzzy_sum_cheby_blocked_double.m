function S = pairwise_fuzzy_sum_cheby_blocked_double(V, r, n, blockSize, kernelType)
% Exact sum_{i<j} kernel(d_ij), with Chebyshev distance, using blocks.

if nargin < 5 || isempty(kernelType)
    kernelType = 'exponential';
end
kernelType = lower(kernelType);

m = size(V, 1);
p = size(V, 2);

S = 0;

for i0 = 1:blockSize:m-1
    i1 = min(i0 + blockSize - 1, m - 1);
    Vi = V(i0:i1, :);
    nb = size(Vi, 1);

    % Within-block upper triangle
    for ii = 1:nb-1
        D = max(abs(Vi(ii+1:end, :) - Vi(ii, :)), [], 2);
        S = S + apply_kernel_sum(D, r, n, kernelType);
    end

    % Across-blocks
    for j0 = i1+1:blockSize:m
        j1 = min(j0 + blockSize - 1, m);
        Vj = V(j0:j1, :);
        nj = size(Vj, 1);

        Dmax = zeros(nb, nj);
        for k = 1:p
            Dmax = max(Dmax, abs(Vi(:, k) - Vj(:, k)'));
        end

        S = S + apply_kernel_sum(Dmax, r, n, kernelType);
    end
end
end

% -------------------------------------------------------------------------
function s = apply_kernel_sum(D, r, n, kernelType)
switch kernelType
    case 'exponential'
        if n == 2
            s = sum(exp(-(D .* D) / r), 'all');
        else
            s = sum(exp(-(D .^ n) / r), 'all');
        end
    case 'gaussian'
        s = sum(exp(-(D .* D) / (2 * r * r)), 'all');
    otherwise
        error('pairwise_fuzzy_sum_cheby_blocked_double:BadKernel', ...
            'Unknown kernelType "%s".', kernelType);
end
end