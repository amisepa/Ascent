function idx = surrogate_segment(X, maxTrim)
% surrogate_segment  Samples to keep so that a signal's two ends match.
%
%   idx = surrogate_segment(X)        % X: [channels x samples]
%   idx = surrogate_segment(X, 0.1)   % trim at most 10% (default)
%
% Fourier surrogates treat a signal as periodic, so a jump between its first
% and last samples adds broadband power to them, and a test against them
% rejects too often. This picks the segment, common to all channels, whose
% ends match best in value and slope (Theiler et al. 1992), cutting at most
% maxTrim/2 of the samples at each end. Channels are z-scored first so that
% each counts the same; channels without finite data are ignored.
%
%   Theiler, J., Eubank, S., Longtin, A., Galdrikian, B., & Farmer, J. D. (1992).
%     Testing for nonlinearity in time series: The method of surrogate data.
%     Physica D, 58, 77-94.
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).
if nargin < 2 || isempty(maxTrim), maxTrim = 0.1; end
if isvector(X), X = X(:).'; end
X = double(X);
N = size(X,2);
K = min(round(maxTrim/2 * N), 1000);   % candidate start and end points
idx = 1:N;
if K < 2, return, end
X = X(all(isfinite(X), 2), :);
if isempty(X), return, end
X = (X - mean(X,2)) ./ max(std(X,0,2), eps);
s = 1:K;                 % starts (the slope uses s and s+1)
e = N-K+1:N;             % ends (the slope uses e-1 and e)
A = [X(:,s); X(:,s+1) - X(:,s)];        % [2*nCh x K]
B = [X(:,e); X(:,e) - X(:,e-1)];        % [2*nCh x K]
% squared distance between every start and every end, as one product
err = sum(A.^2, 1).' + sum(B.^2, 1) - 2 * (A.' * B);   % [K x K]
[~, i] = min(err(:));
[is, ie] = ind2sub([K K], i);
idx = s(is):e(ie);
end
