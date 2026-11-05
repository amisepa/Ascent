function [SampEn, counts] = compute_SampEn(data, varargin)
% compute_SampEn  Computes Sample Entropy (SampEn) across multichannel data.
%
%   [SampEn, counts] = compute_SampEn(data, 'm', 2, 'tau', 1, 'r', .15, ...
%       'ZScore', true, 'BlockSize', 2000, 'pdistMaxGB', 2.0, ...
%       'Parallel', true, 'Progress', true)
%
% Inputs:
%   data        : EEG data matrix [n_channels x n_samples]
%   'm'         : embedding dimension (default = 2)
%   'tau'       : time lag (default = 1)
%   'r'         : similarity bound (default = 0.15)
%   'ZScore'    : z-score each channel before computing SampEn (default = true).
%                 Set false when the caller has already normalised the data and
%                 'r' is an absolute Chebyshev threshold (used by compute_MSE).
%   'BlockSize', 'pdistMaxGB': no longer used (accepted so that existing
%                 calls still work)
%   'Parallel'  : logical true/false to enable parfor over channels (default = true)
%   'Progress'  : logical true/false to show progress
%
% Output:
%   SampEn      : [n_channels x 1] Sample Entropy per channel
%   counts      : [n_channels x 4] the counts behind each value,
%                 [matches(m+1) pairs(m+1) matches(m) pairs(m)], so that
%                 SampEn = -log((c1/c2) / (c3/c4)). Summing counts over
%                 segments before the log gives a pooled estimate (used by
%                 compute_mMSE across starting points).
%
% Notes:
%   • When ZScore is true, data is z-scored per channel across time and r is a
%     fraction of the (unit) SD. When false, data is used as-is and r is the
%     absolute Chebyshev tolerance.
%   • Match counts are exact (sampen_exact: templates sorted by their first
%     coordinate, O(N) memory); they equal comparing every pair with pdist.
%
% References:
%   Richman, J. S., & Moorman, J. R. (2000). Physiological time-series analysis using
%     approximate entropy and sample entropy. Am J Physiol Heart Circ Physiol, 278,
%     H2039-H2049.
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).

% ---------------- Parse inputs ----------------
p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ndims(x) == 2);
p.addParameter('m', 2,           @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('tau', 1,         @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('r', 0.15,        @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('ZScore', true,   @(x) islogical(x) && isscalar(x));
p.addParameter('BlockSize', 2000,@(x) isnumeric(x) && isscalar(x) && x >= 100);
p.addParameter('pdistMaxGB', 2.0,@(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Parallel', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true, @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
tau          = p.Results.tau;
r            = p.Results.r;
doZScore     = p.Results.ZScore;
parallelMode = p.Results.Parallel;
showProgress = p.Results.Progress;

if size(data,1) > size(data,2)
    data = data.';
end
[nchan, ~] = size(data);

% ---------------- Z-score per channel (optional) ----------------
data_z = data;
if doZScore
    for c = 1:nchan
        x  = data(c,:);
        mu = mean(x, 'omitnan');
        sd = std(x, 0, 'omitnan');
        if ~isfinite(sd) || sd == 0
            data_z(c,:) = 0;
        else
            data_z(c,:) = (x - mu) ./ sd;
        end
    end
end

SampEn = nan(nchan, 1);
counts = zeros(nchan, 4);

% ---------------- Progress header ----------------
if showProgress
    if parallelMode && ~isempty(ver('parallel'))
        fprintf('SampEn: %d channel(s) | m=%g, tau=%g, r=%g | parallel=on\n', ...
            nchan, m, tau, r);
    else
        fprintf('SampEn: %d channel(s) | m=%g, tau=%g, r=%g | parallel=off\n', ...
            nchan, m, tau, r);
    end
end

useWB = ~parallelMode && usejava('desktop') && showProgress;
hWB = [];
if useWB
    try
        hWB = waitbar(0, 'Computing Sample Entropy...', 'Name', 'compute_SampEn');
    catch
        hWB = [];
    end
end

useDQ = parallelMode && ~isempty(ver('parallel')) && showProgress;
if useDQ
    dq = parallel.pool.DataQueue;
    nDone = 0;
    afterEach(dq, @notifyProgress);
else
    dq = [];   % placeholder: referenced in parfor but only sent to when useDQ
end

% ---------------- Compute per channel ----------------
if parallelMode && ~isempty(ver('parallel'))
    parfor iChan = 1:nchan
        [SampEn(iChan), counts(iChan,:)] = sampen_exact(data_z(iChan,:), m, r, tau);
        if useDQ
            send(dq, 1);
        end
    end
else
    for iChan = 1:nchan
        [SampEn(iChan), counts(iChan,:)] = sampen_exact(data_z(iChan,:), m, r, tau);

        if ~isempty(hWB) && isvalid(hWB)
            try
                waitbar(iChan/nchan, hWB, sprintf('Computing Sample Entropy... (%d/%d)', iChan, nchan));
            catch
            end
        end
    end

    if ~isempty(hWB) && isvalid(hWB)
        try
            close(hWB);
        catch
        end
    end
end

    function notifyProgress(~)
        nDone = nDone + 1;
        step = max(1, round(0.05 * nchan));
        if nDone == 1 || nDone == nchan || mod(nDone, step) == 0
            fprintf('  progress: ch %d/%d\n', nDone, nchan);
        end
    end
end
