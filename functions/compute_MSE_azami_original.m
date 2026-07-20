% Computes multiscale sample entropy (MSE), reproducing the coarse-graining
% and similarity-threshold logic of the original compute_mse.m /
% compute_se.m pair (Cedric Cannard, 2022), which does NOT rescale r at
% each scale (r is held fixed at whatever value is passed in, applied to
% the z-normalized signal).
%
% This is a direct, line-for-line port of compute_se.m as a local
% subfunction (sample_entropy_original below), NOT a reimplementation
% from a general SampEn definition. This matters because compute_se.m
% has specific indexing and normalization conventions (a double
% division by (n-m), a strict "<" rather than "<=" match criterion,
% index bounds tied to n-m rather than n-k) that a fresh reimplementation
% can easily get subtly wrong. Porting it directly removes that risk.
%
% Purpose: isolate whether the divergence between ASCENT's compute_MSE.m
% and this pre-existing implementation is driven specifically by the
% per-scale rescaling of r (ASCENT rescales; this does not), independent
% of any other implementation differences.
%
% INPUTS:
%   signal      - univariate signal, a vector 1 x N
%   m           - embedding dimension (default = 2)
%   r           - similarity threshold, applied to the z-normalized
%                 signal and held FIXED across all scales (default = 0.15)
%   tau         - time lag (default = 1). Passed to compute_se's own
%                 downsampling step, matching the original convention.
%   coarseType  - 'Mean', 'SD', or 'Variance' (default = 'Mean')
%   nScales     - number of scale factors (default = 20)
%   verbose     - if true, prints per-scale diagnostics (default = false)
%
% OUTPUTS:
%   mse     - 1 x nScales vector of entropy values. NaN at scale 1 for
%             SD/Variance coarse-graining (mathematically degenerate,
%             consistent with ASCENT's own documented exclusion of
%             scale 1 for these operators), and NaN wherever no matches
%             are found at either embedding dimension.
%   scales  - the scale factors corresponding to each mse value
%
% Cedric Cannard, 2026 (diagnostic/validation function; ported from
% compute_mse.m and compute_se.m, Cedric Cannard, 2022)

function [mse, scales] = compute_MSE_azami_original(signal, m, r, tau, coarseType, nScales, verbose)

if ~exist('m','var') || isempty(m)
    m = 2;
end
if ~exist('r','var') || isempty(r)
    r = 0.15;
end
if ~exist('tau','var') || isempty(tau)
    tau = 1;
end
if ~exist('coarseType','var') || isempty(coarseType)
    coarseType = 'Mean';
end
if ~exist('nScales','var') || isempty(nScales)
    nScales = 20;
end
if ~exist('verbose','var') || isempty(verbose)
    verbose = false;
end

% Resolve coarse-graining label (matches compute_mse.m convention)
cl = lower(strtrim(coarseType));
if contains(cl, 'standard') || strcmpi(cl,'sd')
    coarseType = 'SD';
elseif strcmpi(cl,'variance') || strcmpi(cl,'var')
    coarseType = 'Variance';
else
    coarseType = 'Mean';
end

% Center and normalize signal to SD = 1 (matches compute_mse.m exactly:
% signal = signal - mean(signal); signal = signal./std(signal))
signal = signal - mean(signal);
signal = signal ./ std(signal);

mse    = nan(1, nScales);
scales = 1:nScales;

for iScale = 1:nScales

    sig = signal;

    % Coarse-grain via reshape, exactly as in compute_mse.m
    L = floor(length(sig)/iScale) * iScale;
    y = reshape(sig(1:L), iScale, []);

    switch coarseType
        case 'Mean'
            x = mean(y, 'omitnan');
        case 'SD'
            x = std(y, 'omitnan');
        case 'Variance'
            x = var(y, 'omitnan');
    end

    % r is NOT rescaled here, matching compute_mse.m -> compute_se.m,
    % where r is passed straight through unchanged at every scale.
    r_fixed = r;

    % Scale 1 is mathematically degenerate for SD/Variance coarse-graining
    % (SD or variance of a single-sample bin is identically zero), so we
    % skip it explicitly rather than let it silently return an
    % uninformative constant sequence. This mirrors the exclusion already
    % documented for ASCENT's own compute_MSE.m (Section 2.2 of the
    % manuscript).
    if iScale == 1 && ~strcmp(coarseType,'Mean')
        mse(iScale) = NaN;
        if verbose
            fprintf('  scale %2d | SKIPPED (SD/Variance coarse-graining undefined at scale 1)\n', iScale);
        end
        continue
    end

    [mse(iScale), Acount, Bcount] = sample_entropy_original(x, m, r_fixed, tau);

    if verbose
        fprintf('  scale %2d | x SD=%.4f | r=%.4f | matches(m)=%d | matches(m+1)=%d | mse=%.4f\n', ...
            iScale, std(x,'omitnan'), r_fixed, Bcount, Acount, mse(iScale));
    end

end

end

% =========================================================================
% Direct line-for-line port of compute_se.m (Cedric Cannard, 2022),
% based on Azami & Escudero (2016) and Richman & Moorman (2000).
% No re-normalization of the input is performed here; r is used exactly
% as passed in, matching the original.
%
% Returns Acount and Bcount as diagnostic proxies: since the original
% code accumulates a doubly-normalized quantity (count summed and
% divided by (n-m) twice) rather than a raw integer match count, Acount
% and Bcount here are the RAW total match counts (before normalization)
% for embedding dimensions m+1 and m respectively, provided purely for
% diagnostic/verbose printing. The returned entropy value itself follows
% the original function's exact arithmetic, using p(1)/p(2) as defined
% in compute_se.m, not a recomputation from Acount/Bcount.
% =========================================================================
function [entropy, Acount, Bcount] = sample_entropy_original(signal, m, r, tau)

if tau > 1
    signal = signal(1:tau:end);  % equivalent to compute_se.m's downsamp with phase=0
end

n = length(signal);

if n <= m+1
    entropy = NaN; Acount = 0; Bcount = 0;
    return
end

p = zeros(1,2);
sMat = zeros(m+1, n-m);
for i = 1:m+1
    sMat(i,:) = signal(i:n-m+i-1);
end

rawCounts = zeros(1,2);  % raw match totals for k=m and k=m+1, for diagnostics only

for k = m:m+1
    count = zeros(1, n-m);
    tempMat = sMat(1:k,:);

    for i = 1:n-k
        % Exact port of compute_se.m's distance and match logic
        dist = max(abs(tempMat(:,i+1:n-m) - repmat(tempMat(:,i),1,n-m-i)));
        D = (dist < r);
        count(i) = sum(D) / (n-m);
        rawCounts(k-m+1) = rawCounts(k-m+1) + sum(D);
    end

    p(k-m+1) = sum(count) / (n-m);
end

Bcount = rawCounts(1);  % raw matches at embedding dimension m
Acount = rawCounts(2);  % raw matches at embedding dimension m+1

if p(1) == 0 || p(2) == 0
    entropy = NaN;
else
    entropy = log(p(1)/p(2));
end

end