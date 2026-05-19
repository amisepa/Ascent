function [mvFuzzEn, phi_m, phi_m1] = compute_mvFuzzEn(data, varargin)
% compute_mvFuzzEn  Computes uniscale Multivariate Fuzzy Entropy (mvFuzzEn).
%
%   [mvFuzzEn, phi_m, phi_m1] = compute_mvFuzzEn(data, ...
%       'm', 2, 'tau', 1, 'r', .15, 'n', 2, ...
%       'Kernel', 'exponential', 'BlockSize', 2000, 'pdistMaxGB', 2.0, ...
%       'Parallel', true, 'Progress', true)
%
% Inputs:
%   data        : multichannel data [n_channels x n_samples]
%   'm'         : embedding dimension per channel, or scalar (default = 2)
%   'tau'       : time lag per channel, or scalar (default = 1)
%   'r'         : similarity bound (default = 0.15)
%   'n'         : fuzzy exponent (default = 2)
%   'Kernel'    : 'exponential' [default] | 'gaussian'
%   'BlockSize' : block size for bounded-memory pairwise accumulation (default = 2000)
%   'pdistMaxGB': max GB allowed for pdist temporary vector (default = 2.0)
%   'Parallel'  : use parallel computation for m+1 subspaces when available (default = true)
%   'Progress'  : display progress (default = true)
%
% Outputs:
%   mvFuzzEn    : scalar multivariate fuzzy entropy
%   phi_m       : mean fuzzy similarity in embedding dimension m
%   phi_m1      : mean fuzzy similarity in embedding dimension m+1
%
% Notes:
%   • Data are z-score normalized per channel prior to embedding.
%   • Pairwise similarities are computed exactly using pdist when memory allows,
%     otherwise an exact blocked fallback is used.
%   • This implementation follows the multivariate fuzzy entropy formulation
%     used in the RCmvMFE reference, where phi_m1 is averaged across the K
%     augmented (m+1) subspaces.
%   • Entropy values are not constrained to be positive; negative values are valid.
%
% References: 
%	 Fuzzy entropy formulation, kernel, and phi_m+1 averaging:
%	Azami, H & Escudero, J 2017, 'Refined Composite Multivariate Generalized 
% Multiscale Fuzzy Entropy: A Tool for Complexity Analysis of Multichannel Signals',
% Physica a-Statistical mechanics and its applications, vol. 465, pp. 261-276
%
%   Multivariate embedding structure:
%	Ahmed, M. U., & Mandic, D. P. (2011). Multivariate multiscale entropy: 
%A tool for complexity analysis of multichannel data. Physical Review 
%E—Statistical, Nonlinear, and Soft Matter Physics, 84(6), 061918.
%
% Adapted and improved for the ASCENT EEGLAB plugin by Cedric Cannard, 2025.


% ---------- Parse inputs ----------
p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ndims(x) == 2);
p.addParameter('m', 2,                  @(x) isnumeric(x) && all(x(:) > 0));
p.addParameter('tau', 1,                @(x) isnumeric(x) && all(x(:) > 0));
p.addParameter('r', 0.15,               @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('n', 2,                  @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Kernel', 'exponential', @(s) any(strcmpi(s, {'exponential','gaussian'})));
p.addParameter('BlockSize', 2000,       @(x) isnumeric(x) && isscalar(x) && x >= 100);
p.addParameter('pdistMaxGB', 2.0,       @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Parallel', true,        @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,        @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m_in       = p.Results.m;
tau_in     = p.Results.tau;
r          = p.Results.r;
n_exp      = p.Results.n;
kernelType = lower(p.Results.Kernel);
blockSize  = p.Results.BlockSize;
pdistMaxGB = p.Results.pdistMaxGB;
doPar      = p.Results.Parallel && ~isempty(ver('parallel'));
showProg   = p.Results.Progress;

% ---------- Enforce [n_channels x n_samples] ----------
if size(data,1) > size(data,2)
    data = data.';
end
[nvar, nsamp] = size(data);

% ---------- Expand m and tau ----------
M   = expand_param_vec(m_in,   nvar, 'm');
tau = expand_param_vec(tau_in, nvar, 'tau');

% ---------- Z-score per channel ----------
X = data;
for ch = 1:nvar
    x  = X(ch,:);
    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');
    if ~isfinite(sd) || sd == 0
        X(ch,:) = 0;
    else
        X(ch,:) = (x - mu) ./ sd;
    end
end

% ---------- Embedding horizon ----------
N = nsamp - max(M) * max(tau);
if N < 3
    mvFuzzEn = NaN;
    phi_m    = NaN;
    phi_m1   = NaN;
    return
end

if showProg
    fprintf('mvFuzzEn: %d ch | N=%d | sum(m)=%d | r=%g | n=%g | kernel=%s\n', ...
        nvar, N, sum(M), r, n_exp, kernelType);
end

if ~strcmp(kernelType, 'exponential')
    error('compute_mvFuzzEn:KernelNotSupported', ...
        'Shared exact helper currently supports only the exponential kernel.');
end

% ---------- phi_m ----------
A = embed_mv(X, M, tau, N);
A = A - mean(A, 2);
phi_m = fuzzy_pairmean_cheby_exact(A, r, n_exp, blockSize, pdistMaxGB);

% ---------- phi_m1 ----------
Mplus = repmat(M, nvar, 1) + eye(nvar);
phi_m1_vec = nan(nvar, 1);

if doPar
    Bcells = cell(nvar, 1);
    parfor h = 1:nvar
        Bh = embed_mv(X, Mplus(h,:), tau, N);
        Bcells{h} = Bh - mean(Bh, 2);
    end

    for h = 1:nvar
        phi_m1_vec(h) = fuzzy_pairmean_cheby_exact(Bcells{h}, r, n_exp, blockSize, pdistMaxGB);
        if showProg
            fprintf('  subspace %3d/%3d: phi_m1=%.6f\n', h, nvar, phi_m1_vec(h));
        end
    end
else
    for h = 1:nvar
        Bh = embed_mv(X, Mplus(h,:), tau, N);
        Bh = Bh - mean(Bh, 2);
        phi_m1_vec(h) = fuzzy_pairmean_cheby_exact(Bh, r, n_exp, blockSize, pdistMaxGB);
        if showProg
            fprintf('  subspace %3d/%3d: phi_m1=%.6f\n', h, nvar, phi_m1_vec(h));
        end
    end
end

phi_m1 = mean(phi_m1_vec(isfinite(phi_m1_vec)));

% ---------- Entropy ----------
if phi_m > 0 && phi_m1 > 0 && isfinite(phi_m) && isfinite(phi_m1)
    mvFuzzEn = log(phi_m / phi_m1);
else
    mvFuzzEn = NaN;
end
end

% =========================================================================
function v = expand_param_vec(x, n, name)
x = x(:).';
if isscalar(x)
    v = repmat(x, 1, n);
elseif numel(x) == n
    v = x;
else
    error('compute_mvFuzzEn:Bad%s', name, ...
        '%s must be scalar or have length n_channels.', name);
end
end

% =========================================================================
function A = embed_mv(X, M, tau, N)
A = zeros(N, sum(M));

col = 0;
for ch = 1:numel(M)
    for k = 0:(M(ch)-1)
        col = col + 1;
        idx = (1:N) + k * tau(ch);
        A(:, col) = X(ch, idx).';
    end
end
end