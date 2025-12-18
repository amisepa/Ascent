function [RCmvMFE, scales] = compute_RCmvMFE(data, varargin)
% compute_RCmvMFE  Refined Composite Multivariate Multiscale Fuzzy Entropy (Azami & Escudero 2016).
% paper: 
%
%   [RCmvMFE, scales] = compute_RCmvMFE(data, 'm', 2, 'r', 0.15, 'tau', 1, ...
%                                       'n', 2, 'coarsing','mean', 'num_scales', 15, ...
%                                       'Parallel', true, 'Progress', true)
%
% Inputs (name-value):
%   data       : multivariate signal [n_var x n_samp] OR EEGLAB EEG struct
%   'm'        : embedding dimension (default = 2)
%   'r'        : similarity bound (default = 0.15; signals z-scored per channel)
%   'tau'      : time lag (default = 1)
%   'n'        : fuzzy exponent (default = 2)
%   'coarsing' : 'mean' | 'var' [default = 'mean']
%   'num_scales': requested number of scales (default = 15)
%   'Parallel' : parfor over scales (default = true)
%   'Progress' : header / progress reporting (default = true)
%
% Output:
%   RCmvMFE : [1 x S] refined composite multivariate fuzzy entropy values
%   scales  : 1:S (simple numeric scale index)
%
% Reference:
%  Azami, H., & Escudero, J. (2017). Refined composite multivariate
% generalized multiscale fuzzy entropy: A tool for complexity analysis of 
% multichannel signals. Physica A: Statistical Mechanics and its Applications, 465, 261-276.
% https://www.sciencedirect.com/science/article/abs/pii/S0378437116305404
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Ascent plugin — Author: Cedric Cannard (adapted from Azami & Escudero)
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

% ---------- Parse inputs ----------
p = inputParser;
p.addRequired('data', @(x) (isstruct(x) && isfield(x,'data')) || isnumeric(x));
p.addParameter('m', 2, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('r', 0.15, @(x) isnumeric(x) && isscalar(x) && x>0 && x<2);
p.addParameter('tau', 1, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('n', 2, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('coarsing','mean', @(s) any(strcmpi(s,{'mean','var'})));
p.addParameter('num_scales', 15, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('Parallel', false, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true, @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
r_base       = p.Results.r;
tau          = p.Results.tau;
n_exp        = p.Results.n;
coarseType   = p.Results.coarsing;
num_scales  = p.Results.num_scales;
parallelMode = p.Results.Parallel;
showProg     = p.Results.Progress;

% ---------- Coerce data ----------
if isstruct(data), X = double(data.data); else, X = double(data); end
if size(X,1) > size(X,2), X = X.'; end
[nChan, nSamp] = size(X);

% Normalize channels to unit variance
X = zscore(X,0,2);

% Adjust r to multivariate scale
r = r_base * sum(std(X,0,2));

% Embedding & tau vectors
M   = m * ones(1,nChan);
Tau = tau * ones(1,nChan);

if showProg
    if parallelMode
        parStr = 'on';
    else
        parStr = 'off';
    end
    fprintf('RCmvMFE: channels=%d | m=%g, tau=%g, r=%g, n=%g | coarse=%s | scales=%d | parallel=%s\n', ...
        nChan, m, tau, r, n_exp, coarseType, num_scales, parStr);
end

% Initialize
scales = 1:num_scales; % output
RCmvMFE = nan(1,num_scales);

% % Compute
% isolates the per-scale computation so the same code can be executed 
% safely in either a for or parfor loop without duplicating logic or 
% violating MATLAB's parallel loop constraints.
if parallelMode && ~isempty(ver('parallel'))
    parfor iScale = 1:num_scales
        RCmvMFE(iScale) = one_scale_rcmvmfe(X, M, r, n_exp, Tau, iScale, coarseType, showProg, num_scales);
    end
else
    for iScale = 1:num_scales
        RCmvMFE(iScale) = one_scale_rcmvmfe(X, M, r, n_exp, Tau, iScale, coarseType, showProg, num_scales);
    end
end


% % Parallel computing mode
% if parallelMode && ~isempty(ver('parallel'))
%     parfor iScale = 1:num_scales
%         tScale = [];
%         if showProg
%             fprintf('  progress: computing for scale %d/%d...\n', iScale, num_scales);
%             tScale = tic;
%         end
% 
%         fi_m1_vec = zeros(1, iScale);
%         fi_m2_vec = zeros(1, iScale);
% 
%         for offset = 1:iScale
%             Xseg = X(:, offset:end);
%             Xs   = coarsegrain(Xseg, iScale, coarseType);
% 
%             [~, fi_m1, fi_m2] = compute_mvFE_fast(Xs, M, r, n_exp, Tau, ...
%                 'blockSize', 2000, 'useSingle', true, 'progress', false);
% 
%             fi_m1_vec(offset) = fi_m1;
%             fi_m2_vec(offset) = fi_m2;
%         end
% 
%         AA = sum(fi_m1_vec);
%         BB = sum(fi_m2_vec);
% 
%         if BB > 0 && isfinite(AA) && isfinite(BB)
%             RCmvMFE(iScale) = log(AA / BB);
%         else
%             RCmvMFE(iScale) = NaN;
%         end
% 
%         if showProg
%             fprintf('  progress: finished scale %d/%d (%.1fs)\n', iScale, num_scales, toc(tScale));
%         end
%     end
% 
% 
% else % No parallel computing
%     for iScale = 1:num_scales
%         if showProg
%             fprintf('  progress: computing for scale %d/%d...\n', iScale, num_scales);
%             tScale = tic;
%         end
% 
%         fi_m1_vec = zeros(1, iScale);
%         fi_m2_vec = zeros(1, iScale);
% 
%         for offset = 1:iScale
%             Xseg = X(:, offset:end);
%             Xs   = coarsegrain(Xseg, iScale, coarseType);
% 
%             [~, fi_m1, fi_m2] = compute_mvFE_fast(Xs, M, r, n_exp, Tau, ...
%                 'blockSize', 2000, 'useSingle', true, 'progress', false);
% 
%             fi_m1_vec(offset) = fi_m1;
%             fi_m2_vec(offset) = fi_m2;
%         end
% 
%         AA = sum(fi_m1_vec);
%         BB = sum(fi_m2_vec);
% 
%         if BB > 0 && isfinite(AA) && isfinite(BB)
%             RCmvMFE(iScale) = log(AA / BB);
%         else
%             RCmvMFE(iScale) = NaN;
%         end
% 
%         if showProg
%             fprintf('  progress: finished scale %d/%d (%.1fs)\n', iScale, num_scales, toc(tScale));
%         end
%     end
% end
    
end % main function



%% HELPERS

function val = one_scale_rcmvmfe(X, M, r, n_exp, Tau, iScale, coarseType, showProg, num_scales)
% isolates the per-scale computation so the same code can be executed 
% safely in either a for or parfor loop without duplicating logic or 
% violating MATLAB's parallel loop constraints.
tScale = [];
if showProg
    fprintf('  progress: computing for scale %d/%d...\n', iScale, num_scales);
    tScale = tic;
end
fi_m1_vec = zeros(1, iScale);
fi_m2_vec = zeros(1, iScale);
for offset = 1:iScale
    Xs = coarsegrain(X(:, offset:end), iScale, coarseType);

    % [~, fi_m1, fi_m2] = compute_mvFE_fast(Xs, M, r, n_exp, Tau, ...
    %     'blockSize', 2000, 'useSingle', true, 'progress', false);
    [~, fi_m1, fi_m2] = compute_mvFE_block(Xs, M, r, n_exp, Tau, 2000); % memory block method

    fi_m1_vec(offset) = fi_m1;
    fi_m2_vec(offset) = fi_m2;
end
AA = sum(fi_m1_vec);
BB = sum(fi_m2_vec);
if BB > 0 && isfinite(AA) && isfinite(BB)
    val = log(AA / BB);
else
    val = NaN;
end
if showProg && ~isempty(tScale)
    fprintf('  progress: finished scale %d/%d (%.1fs)\n', iScale, num_scales, toc(tScale));
end
end


function Y = coarsegrain(Data, S, method)
% Coarse graining of multichannel data by mean or variance
% Inputs:
%   Data   : nvar x nsamp
%   S      : scale factor (positive integer)
%   method : 'mean' or 'var'
% Output:
%   Y      : nvar x floor(nsamp/S)

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
        Y = squeeze(mean(Data, 2, 'omitnan'));
    case 'var'
        % If you want omitnan variance and your MATLAB supports it, use:
        % Y = squeeze(var(Data, 0, 2, 'omitnan'));
        Y = squeeze(var(Data, 0, 2));
    otherwise
        error('method must be ''mean'' or ''var''.');
end
end





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
p = inputParser;
p.addParameter('blockSize', 2000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
p.addParameter('useSingle', true, @(x) islogical(x) && isscalar(x));
p.addParameter('progress',  false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
blockSize = p.Results.blockSize;
useSingle = p.Results.useSingle;
showProg  = p.Results.progress;
if iscolumn(M),   M   = M.';   end
if iscolumn(tau), tau = tau.'; end
[nvar, nsamp] = size(X);
max_M   = max(M);
max_tau = max(tau);

% Dimension m embedding
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

% Dimension m+1 embedding across nvar subspaces
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

% Truncate to common length to avoid off by 1 mismatches
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

mvFE_Out = log(fi_m1 / fi_m2);
end



%--------------------------------------------------------------------------
function [mvFE_Out, fi_m1, fi_m2] = compute_mvFE_block(X, M, r, n, tau, blockSize)
% Memory-safe mvFE with blocked pairwise sums (Chebyshev)
% Handles m and m+1 potentially having different N.
if nargin < 7 || isempty(blockSize)
    blockSize = 2000;
end
[nvar, nsamp] = size(X);

% Dimension m
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

% Dimension m+1 across nvar subspaces
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
S2 = pairwise_fuzzy_sum_cheby_block(B, r, n, blockSize, false);
fi_m2 = S2 * 2 / (NN*(NN-1));
mvFE_Out = log(fi_m1 / fi_m2);
end




% ========================================================================
function S = pairwise_fuzzy_sum_cheby_block(V, r, n, blockSize, showProg)
% Sum over i<j of exp(-(d_ij^n)/r), Chebyshev distance
%
% V: [m x p]
m = size(V,1);
p = size(V,2);
S = 0;

%  Blocked accumulation
for i0 = 1:blockSize:m-1
    i1 = min(i0 + blockSize - 1, m-1);
    Vi = V(i0:i1,:);
    nb = size(Vi,1);
    if showProg && i0 == 1
        fprintf('  pairwise blocks: %d rows, blockSize=%d\n', m, blockSize);
    end

    %  Within block upper triangle
    for ii = 1:nb-1
        D = max(abs(Vi(ii+1:end,:) - Vi(ii,:)), [], 2);
        S = S + sum(exp(-(D.^n)/r), 'omitnan');
    end

    %  Cross blocks
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


function A=embd(M,tau,ts)
% This function creates multivariate delay embedded vectors with embedding
% vector parameter M and time lag vector parameter tau.
% M is a row vector [m1 m2 ...mnvar] and tau is also a row vector [tau1 tau2....taunvar] where nvar is the
% number of channels;
% ts is the multivariate time series-a matrix of size nvarxnsamp;
% Ref: M. U. Ahmed and D. P. Mandic, "Multivariate multiscale entropy
% analysis", IEEE Signal Processing Letters, vol. 19, no. 2, pp.91-94.2012
[nvar,nsamp]=size(ts);
A=[];
for j=1:nvar
    for i=1:nsamp-max(M)
        temp1(i,:)=ts(j,i:tau(j):i+M(j)-1);
    end
    A=horzcat(A,temp1);
    temp1=[];
end
end