function [mvFuzzEn, phi_m, phi_m1] = compute_mvFuzzEn(data, varargin)
% compute_mvFuzzEn  Computes uniscale Multivariate Fuzzy Entropy (mvFuzzEn).
%
%   [mvFuzzEn, phi_m, phi_m1] = compute_mvFuzzEn(data, ...
%       'M', 2, 'tau', 1, 'r', .15, 'n', 2, ...
%       'Kernel','exponential', 'BlockSize', 256, ...
%       'Parallel', false, 'Progress', false)
%
% Inputs:
%   data        : multichannel data [n_channels x n_samples] (numeric)
%   'm'         : embedding vector per channel OR scalar (default = 2)
%   'tau'       : time lag vector per channel OR scalar (default = 1)
%   'r'         : similarity bound (default = .15)
%   'n'         : fuzziness exponent for exponential kernel (default = 2)
%   'Kernel'    : 'exponential' [default] | 'gaussian' (note: Azami-style uses exponential)
%   'BlockSize' : pair-block size to bound memory (default = 256; >=128)
%   'Parallel'  : logical true/false to enable parfor for building m+1 embeddings (default = false)
%   'Progress'  : logical true/false to show text progress (default = false)
%
% Outputs:
%   mvFuzzEn    : scalar mvFuzzEn of the multivariate signal
%   phi_m       : mean fuzzy similarity in dimension m
%   phi_m1      : mean fuzzy similarity in dimension m+1 (averaged across nvar subspaces)
%
% Notes:
%   • Data is z-score normalized per channel so std = 1 (stable r across recordings).
%   • Implements the uniscale mvFE structure from Azami & Escudero (Physica A, 2016),
%     but computed with bounded-memory blockwise pair accumulation (no pdist).
%   • The entropy estimate is mvFuzzEn = log(phi_m / phi_m1).
%
% References:
%   Azami, H. & Escudero, J. (2016). Refined Composite Multivariate Generalized
%       Multiscale Fuzzy Entropy. Physica A.
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Ascent plugin — Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

% ---------------- Parse inputs ----------------
p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ndims(x) == 2);
p.addParameter('m', 2,                  @(x) isnumeric(x) && all(x(:) > 0));
p.addParameter('tau', 1,                @(x) isnumeric(x) && all(x(:) > 0));
p.addParameter('r', .15,                @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('n', 2,                  @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Kernel','exponential',  @(s) any(strcmpi(s,{'exponential','gaussian'})));
p.addParameter('BlockSize', 256,        @(x) isnumeric(x) && isscalar(x) && x >= 128 && x <= 4096);
p.addParameter('Parallel', true,        @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,        @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m_in         = p.Results.m;
tau_in       = p.Results.tau;
r            = p.Results.r;
n_exp        = p.Results.n;
kernelType   = lower(p.Results.Kernel);
blk          = p.Results.BlockSize;
doPar        = p.Results.Parallel && ~isempty(ver('parallel'));
showProg     = p.Results.Progress;

% Enforce [n_channels x n_samples]
if size(data,1) > size(data,2)
    data = data.';  % [n_channels x n_samples]
end
[nvar, nsamp] = size(data);

% Expand M and tau to per-channel vectors
M   = expand_param_vec(m_in,   nvar, 'M');
tau = expand_param_vec(tau_in, nvar, 'tau');

% Z-score normalization per channel
data_z = data;
for c = 1:nvar
    x  = data(c,:);
    sd = std(x,0,'omitnan');
    if ~isfinite(sd) || sd == 0
        data_z(c,:) = 0;
    else
        mu = mean(x,'omitnan');
        data_z(c,:) = (x - mu)./sd;
    end
end

% Azami-style effective embedding horizon
max_M   = max(M);
max_tau = max(tau);
nn      = max_M * max_tau;
N       = nsamp - nn;

if N < 3
    mvFuzzEn = NaN; phi_m = NaN; phi_m1 = NaN;
    return;
end

if showProg
    fprintf('mvFuzzEn: %d ch | N=%d | sum(M)=%d | r=%g | n=%g | kernel=%s\n', ...
        nvar, N, sum(M), r, n_exp, kernelType);
end

% ---------------- phi_m (dimension m) ----------------
A = embed_mv(data_z, M, tau, N);
phi_m = fuzzy_mean_similarity_block(A, r, n_exp, kernelType, blk);

% ---------------- phi_m1 (dimension m+1, averaged across nvar subspaces) ----------------
Mmat = repmat(M(:).', nvar, 1) + eye(nvar);  % each row increments one channel's embedding by 1

if doPar
    Bcells = cell(nvar,1);
    parfor h = 1:nvar
        Bcells{h} = embed_mv(data_z, Mmat(h,:), tau, N);
    end
    B = vertcat(Bcells{:});
else
    B = zeros(nvar*N, sum(M)+1);
    row0 = 0;
    for h = 1:nvar
        Bh = embed_mv(data_z, Mmat(h,:), tau, N);
        B((row0+1):(row0+N), :) = Bh;
        row0 = row0 + N;
        if showProg && (h==1 || h==nvar || mod(h, max(1,floor(nvar/10)))==0)
            fprintf('  building m+1 embeddings: %d/%d\n', h, nvar);
        end
    end
end

phi_m1 = fuzzy_mean_similarity_block(B, r, n_exp, kernelType, blk);

% ---------------- mvFuzzEn ----------------
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
else
    if numel(x) ~= n
        error('compute_mvFuzzEn:Bad%s', name, ...
            '%s must be scalar or have length n_channels.', name);
    end
    v = x;
end
end

function A = embed_mv(X, M, tau, N)
% Create multivariate delay-embedding matrix A of size [N x sum(M)].
% X: [nvar x nsamp], M: [1 x nvar], tau: [1 x nvar], N: #rows
[nvar, ~] = size(X);
D = sum(M);
A = zeros(N, D);

col = 0;
for ch = 1:nvar
    m  = M(ch);
    t  = tau(ch);
    for k = 0:(m-1)
        col = col + 1;
        idx = (1:N) + k*t;
        A(:, col) = X(ch, idx).';
    end
end
end

function mu_mean = fuzzy_mean_similarity_block(X, r, n, kernelType, blk)
% Mean fuzzy similarity over unique pairs under Chebyshev (max) distance.
nVec = size(X,1);
if nVec < 2, mu_mean = NaN; return; end

sum_mu = 0;
num_p  = 0;
nDim   = size(X,2);

t0 = tic;
lastPrint = 0;
for i1 = 1:blk:nVec
    if toc(t0) - lastPrint > 5
        fprintf('    pairs progress: block start %d/%d (%.1f%%)\n', i1, nVec, 100*i1/nVec);
        lastPrint = toc(t0);
    end
    
    i2 = min(i1+blk-1, nVec);
    Xi = X(i1:i2,:);  bi = size(Xi,1);

    % (A) Within-block: upper triangle only
    if bi > 1
        for d1 = 1:blk:bi
            d2 = min(d1+blk-1, bi);
            Xd = Xi(d1:d2,:); bd = size(Xd,1);
            if bd > 1
                Dmax = zeros(bd, bd);
                for dim = 1:nDim
                    Dij = abs(Xd(:,dim) - Xd(:,dim).');
                    if dim == 1
                        Dmax = Dij;
                    else
                        Dmax = max(Dmax, Dij);
                    end
                end
                ut = triu(true(bd,bd), 1);
                dvec = Dmax(ut);
                sum_mu = sum_mu + sum(fuzzy_kernel(dvec, r, n, kernelType));
                num_p  = num_p  + nnz(ut);
            end
        end
    end

    % (B) Across-blocks: full rectangle
    for j1 = (i2+1):blk:nVec
        j2 = min(j1+blk-1, nVec);
        Xj = X(j1:j2,:);  bj = size(Xj,1);

        Dmax = zeros(bi, bj);
        for dim = 1:nDim
            Dij = abs(Xi(:,dim) - Xj(:,dim).');
            if dim == 1
                Dmax = Dij;
            else
                Dmax = max(Dmax, Dij);
            end
        end
        sum_mu = sum_mu + sum(fuzzy_kernel(Dmax(:), r, n, kernelType));
        num_p  = num_p  + bi * bj;
    end
end

if num_p == 0
    mu_mean = NaN;
else
    mu_mean = sum_mu / num_p;
end
end

function y = fuzzy_kernel(d, r, n, kernelType)
switch kernelType
    case 'exponential'
        % Azami-style kernel
        y = exp(-(d.^n) / r);
    case 'gaussian'
        y = exp(-(d.^2) / (2*r^2));
    otherwise
        y = exp(-(d.^n) / r);
end
end
