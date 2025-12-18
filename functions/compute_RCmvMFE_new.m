function [RCmvMFE, scales] = compute_RCmvMFE_new(data, varargin)
% compute_RCmvMFE  Refined Composite Multivariate Multiscale Fuzzy Entropy (Azami & Escudero 2016).
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
%   Azami, H., & Escudero, J. (2016). Physica A, 465, 261–276.
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
p.addParameter('Parallel', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true, @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
r_base       = p.Results.r;
tau          = p.Results.tau;
n_exp        = p.Results.n;
coarseType   = p.Results.coarsing;
nScales_req  = p.Results.num_scales;
parallelMode = p.Results.Parallel;
showProg     = p.Results.Progress;

% ---------- Coerce data ----------
if isstruct(data), X = double(data.data); else, X = double(data); end
if size(X,1) > size(X,2), X = X.'; end
[nvar, nSamp] = size(X);

% Normalize channels to unit variance
X = zscore(X,0,2);

% Adjust r to multivariate scale
r = r_base * sum(std(X,0,2));

% Embedding & tau vectors
M   = m * ones(1,nvar);
Tau = tau * ones(1,nvar);

% ---------- Output ----------
S = max(1, floor(nScales_req));
scales = 1:S;
RCmvMFE = nan(1,S);

if showProg
    if parallelMode
        parStr = 'on';
    else
        parStr = 'off';
    end
    fprintf('RCmvMFE: %d vars | m=%g, tau=%g, r=%g, n=%g | coarse=%s | S=%d | parallel=%s\n', ...
        nvar, m, tau, r, n_exp, upper(coarseType), S, parStr);
end

% % ---------- Compute ----------
% if parallelMode && ~isempty(ver('parallel'))
%     parfor s = 2:S
%         RCmvMFE(s) = rc_mvfe_one_scale(X, M, r, n_exp, Tau, s, coarseType);
%     end
% else
for s = 2:S
    if showProg
        fprintf('  progress: starting scale %d/%d\n', s, S);
        tScale = tic;
    end

    RCmvMFE(s) = rc_mvfe_one_scale(X, M, r, n_exp, Tau, s, coarseType);

    if showProg
        fprintf('  progress: finished scale %d/%d (%.1fs)\n', s, S, toc(tScale));
    end
end
% end

end % main

% ========================================================================
function val = rc_mvfe_one_scale(X, M, r, n_exp, Tau, scale, coarseType)

fi_m1_vec = zeros(1, scale);
fi_m2_vec = zeros(1, scale);

for offset = 1:scale
    if strcmpi(coarseType,'mean')
        Xs = coarsegrain(X(:,offset:end), scale, 'mean');   % if you made the unified version
    else
        Xs = coarsegrain(X(:,offset:end), scale, 'var');
    end

    [~, fi_m1, fi_m2] = compute_mvFE_blocked(Xs, M, r, n_exp, Tau); % <- blocked call
    % [~, fi_m1, fi_m2] = compute_mvFE_fast(Xs, M, r, n_exp, Tau, ...
    %     'blockSize', 2000, 'useSingle', true, 'progress', false);

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
end

% ========================================================================
function Y = coarsegrain(Data, S, method)
% 1. Coarse graining of multichannel data by mean or variance
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