function [Xz, bad] = prep_channels(X)
% prep_channels  Fill missing samples and z-score each channel (sample SD),
%   once, for the multiscale measures. X is [channels x samples].
%   bad marks channels with no usable data (all missing, or flat): their row
%   is 0 and the callers report NaN for them, so that one dead channel does
%   not look like a perfectly regular one.
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).
X(~isfinite(X)) = NaN;   % fillmissing leaves Inf in place
Xz  = zeros(size(X), 'like', X);
bad = false(size(X,1), 1);
for c = 1:size(X,1)
    x = X(c,:);
    if all(isnan(x)), bad(c) = true; continue, end
    if any(isnan(x)), x = fillmissing(x, 'linear', 'EndValues', 'nearest'); end
    sd = std(x, 0);
    if ~(sd > 0), bad(c) = true; continue, end
    Xz(c,:) = (x - mean(x)) ./ sd;
end
end
