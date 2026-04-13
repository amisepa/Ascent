function [fe, pm, pm1] = fuzz_engine_raw(signal, m, r, n_exp, tau, kernelType, doNormalize, blockSize, pdistMaxGB)
% fuzz_engine_raw  Low-level Fuzzy Entropy backend for a single series.
%
%   [fe, pm, pm1] = fuzz_engine_raw(signal, m, r, n_exp, tau, kernelType, doNormalize)
%   [fe, pm, pm1] = fuzz_engine_raw(signal, m, r, n_exp, tau, kernelType, doNormalize, blockSize, pdistMaxGB)
%
% Inputs:
%   signal      : numeric vector
%   m           : embedding dimension
%   r           : similarity bound
%   n_exp       : fuzzy exponent
%   tau         : time lag
%   kernelType  : 'exponential' or 'gaussian'
%   doNormalize : true/false, z-score this signal before computing FuzzEn
%   blockSize   : block size for exact blocked fallback (default = 2000)
%   pdistMaxGB  : max GB allowed for pdist temporary vector (default = 2.0)
%
% Outputs:
%   fe   : Fuzzy Entropy = log(pm / pm1)
%   pm   : mean fuzzy similarity for embedding dimension m
%   pm1  : mean fuzzy similarity for embedding dimension m+1

if nargin < 8 || isempty(blockSize)
    blockSize = 2000;
end
if nargin < 9 || isempty(pdistMaxGB)
    pdistMaxGB = 2.0;
end

if ~isrow(signal)
    signal = signal(:).';
end

x = double(signal);
x = x(isfinite(x));

if doNormalize
    sd = std(x, 0, 'omitnan');
    if ~isfinite(sd) || sd == 0
        x = zeros(size(x));
    else
        x = (x - mean(x, 'omitnan')) ./ sd;
    end
end

N = numel(x);
nVec_m  = N - (m   - 1) * tau;
nVec_m1 = N - (m+1 - 1) * tau;

if nVec_m < 2 || nVec_m1 < 2
    fe = NaN;
    pm = NaN;
    pm1 = NaN;
    return
end

Xm  = embed_uni_local(x, m,   tau, nVec_m);
Xm1 = embed_uni_local(x, m+1, tau, nVec_m1);

pm  = fuzzy_pairmean_cheby_exact(Xm,  r, n_exp, blockSize, pdistMaxGB, kernelType);
pm1 = fuzzy_pairmean_cheby_exact(Xm1, r, n_exp, blockSize, pdistMaxGB, kernelType);

if pm > 0 && pm1 > 0 && isfinite(pm) && isfinite(pm1)
    fe = log(pm / pm1);
else
    fe = NaN;
end
end

% -------------------------------------------------------------------------
function X = embed_uni_local(signal, m, tau, nVec)
X = zeros(nVec, m, 'double');
for k = 1:m
    idx = (1:nVec) + (k-1) * tau;
    X(:, k) = signal(idx).';
end
end