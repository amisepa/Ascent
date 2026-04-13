function fi = fuzzy_pairmean_cheby_exact(V, r, n, blockSize, pdistMaxGB, kernelType)
% fuzzy_pairmean_cheby_exact
% Exact mean fuzzy similarity over all unique pairs using Chebyshev distance.
%
%   fi = 2/(M*(M-1)) * sum_{i<j} kernel(d_ij)
%
% Inputs:
%   V          : [nVec x nDim] matrix of embedded vectors
%   r          : similarity bound
%   n          : fuzzy exponent (used for exponential kernel)
%   blockSize  : block size for exact blocked fallback
%   pdistMaxGB : max temporary memory allowed for pdist vector, in GB
%   kernelType : 'exponential' or 'gaussian' (default = 'exponential')

if nargin < 6 || isempty(kernelType)
    kernelType = 'exponential';
end
kernelType = lower(kernelType);

M = size(V, 1);
if M < 2
    fi = NaN;
    return
end

pairs = double(M) * double(M - 1) / 2;
gbNeeded = (pairs * 8) / (1024^3);   % pdist returns doubles

if gbNeeded <= pdistMaxGB
    d = pdist(V, 'chebychev');
    switch kernelType
        case 'exponential'
            if n == 2
                s = sum(exp(-(d .* d) / r));
            else
                s = sum(exp(-(d .^ n) / r));
            end
        case 'gaussian'
            s = sum(exp(-(d .* d) / (2 * r * r)));
        otherwise
            error('fuzzy_pairmean_cheby_exact:BadKernel', ...
                'Unknown kernelType "%s".', kernelType);
    end
else
    s = pairwise_fuzzy_sum_cheby_blocked_double(V, r, n, blockSize, kernelType);
end

fi = s * 2 / (M * (M - 1));
end