function [MSE, scales] = compute_MSE(data, varargin)
% compute_MSE  Multiscale Entropy (Costa-style) across multichannel data.
%
%   [MSE, scales] = compute_MSE(data, 'm', 2, 'r', 0.15, 'tau', 1, ...
%                               'coarsing','std', 'num_scales', 20, ...
%                               'MinSamplesPerBin', 4, 'StableMinBins', 100, ...
%                               'Parallel', true, 'Progress', true)

% ---------------- Parse inputs ----------------
p = inputParser;
p.addRequired('data', @(x) (isstruct(x) && isfield(x,'data')) || (isnumeric(x) && ndims(x)==2));
p.addParameter('m', 2,                 @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('r', 0.15,              @(x) isnumeric(x) && isscalar(x) && x>0 && x<2);
p.addParameter('tau', 1,               @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('coarsing','std',       @(s) any(strcmpi(s,{'median','mean','trimmed mean','trimmed','trimmean','std','sd','standard deviation','var','variance'})));
p.addParameter('num_scales', 20,       @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('MinSamplesPerBin', 4,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('StableMinBins', 100,   @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('Parallel', true,       @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,       @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m             = p.Results.m;
r             = p.Results.r;
tau           = p.Results.tau;
coarseType    = p.Results.coarsing;
nScales_req   = p.Results.num_scales;
minBinsAll    = p.Results.MinSamplesPerBin;
minBinsStable = p.Results.StableMinBins;
parallelMode  = p.Results.Parallel;
showProgress  = p.Results.Progress;

% ---------------- Get numeric data ----------------
if isstruct(data)
    X = double(data.data);
else
    X = double(data);
end
if size(X,1) > size(X,2)
    X = X.';
end
[nch, nSamp] = size(X);

% ---------------- Cap S so every scale has >= minBinsAll coarse bins -----
S = max(1, floor(nScales_req));
maxS = floor(nSamp / max(1, minBinsAll));
if S > maxS
    warning('compute_MSE:ReducingScales', ...
        'Reducing num_scales from %d to %d to keep >=%d coarse bins at every scale.', ...
        S, maxS, minBinsAll);
    S = maxS;
end
if S < 1
    S = 1;
end
scales = 2:S;  % skip scale 1 (sampEn) that is not very comparable to the rest

% ---------------- Fill missing & z-score per channel ---------------------
for ch = 1:nch
    xi = X(ch,:);
    if any(~isfinite(xi))
        try
            xi = fillmissing(xi, 'linear', 'EndValues', 'nearest');
        catch
            isn = ~isfinite(xi);
            if any(isn)
                idx = find(~isn, 1, 'first');
                if ~isempty(idx), xi(1:idx-1) = xi(idx); end
                idx = find(~isn, 1, 'last');
                if ~isempty(idx), xi(idx+1:end) = xi(idx); end
                prev = [xi(1), xi(1:end-1)];
                next = [xi(2:end), xi(end)];
                bad  = isn & isfinite(prev) & isfinite(next);
                xi(bad) = 0.5 * (prev(bad) + next(bad));
            end
        end
        X(ch,:) = xi;
    end
end

Xz = X;
for c = 1:nch
    x  = X(c,:);
    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');
    if ~isfinite(sd) || sd == 0
        Xz(c,:) = 0;
    else
        Xz(c,:) = (x - mu) ./ sd;
    end
end

% ---------------- Outputs & progress header ------------------------------
% MSE = nan(nch, S);
MSE = nan(nch, S-1);


cl = lower(strtrim(coarseType));
if any(strcmp(cl, {'sd','std','standard deviation'}))
    coarseLabel = 'STD';
elseif any(strcmp(cl, {'var','variance'}))
    coarseLabel = 'VAR';
elseif strcmp(cl, 'mean')
    coarseLabel = 'MEAN';
elseif strcmp(cl, 'median')
    coarseLabel = 'MEDIAN';
elseif any(strcmp(cl, {'trimmed mean','trimmed','trimmean'}))
    coarseLabel = 'TRIM20';
else
    coarseLabel = upper(coarseType);
end

if showProgress
    if parallelMode && ~isempty(ver('parallel'))
        fprintf('MSE: %d ch | m=%g, tau=%g, r=%g | coarse=%s | S=%d | parallel=on\n', ...
            nch, m, tau, r, coarseLabel, S);
    else
        fprintf('MSE: %d ch | m=%g, tau=%g, r=%g | coarse=%s | S=%d | parallel=off\n', ...
            nch, m, tau, r, coarseLabel, S);
    end
end

% ---------------- Progress helpers ----------------
useWB = showProgress && ~parallelMode && usejava('desktop');
hWB = [];
if useWB
    try
        hWB = waitbar(0, 'Computing Multiscale Entropy...', 'Name', 'compute_MSE');
    catch
        hWB = [];
    end
end

useDQ = parallelMode && ~isempty(ver('parallel')) && showProgress;
if useDQ
    dq = parallel.pool.DataQueue;
    nDone = 0;
    afterEach(dq, @notifyProgress);
end

% ---------------- Compute per channel ------------------------------------
if parallelMode && ~isempty(ver('parallel'))
    parfor ch = 1:nch
        MSE(ch,:) = mse_one_channel(Xz(ch,:), m, r, tau, coarseType, ...
            S, minBinsAll, minBinsStable, false, ch, nch);
        if useDQ
            send(dq, 1);
        end
    end
else
    for ch = 1:nch
        MSE(ch,:) = mse_one_channel(Xz(ch,:), m, r, tau, coarseType, ...
            S, minBinsAll, minBinsStable, showProgress, ch, nch);

        if ~isempty(hWB) && isvalid(hWB)
            try
                waitbar(ch/nch, hWB, sprintf('Computing MSE... (%d/%d)', ch, nch));
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
        step = max(1, round(0.05 * nch));
        if nDone == 1 || nDone == nch || mod(nDone, step) == 0
            fprintf('  progress: ch %d/%d\n', nDone, nch);
        end
    end
end

% ========================================================================
function v = mse_one_channel(sig, m, r, tau, coarseType, S, ...
    minBinsAll, minBinsStable, showProgress, ch, nch)

% v = nan(1, S);
v = nan(1, S-1);   % scales 2:S => S-1 values

for s = 2:S
    % if s == 1
    %     v(1) = compute_SampEn(sig, ...
    %         'm', m, 'r', r, 'tau', tau, ...
    %         'Parallel', false, 'Progress', false);
    %     continue
    % end

    L = floor(numel(sig)/s) * s;
    nBins = L / s;
    minNeeded = max([minBinsAll, m+1, minBinsStable]);

    if nBins < minNeeded
        if showProgress
            fprintf('  [drop] ch %d: scale %d nBins=%d (<%d)\n', ch, s, nBins, minNeeded);
        end
        continue
    end

    Y  = reshape(sig(1:L), s, []);
    cg = coarsegrain(Y, coarseType);

    % v(s) = compute_SampEn(cg, 'm', m, 'r', r, 'tau', tau, ...
    %     'Parallel', false, 'Progress', false);
    v(s-1) = compute_SampEn(cg, 'm', m, 'r', r, 'tau', tau, ...
        'Parallel', false, 'Progress', false); % when ignoring scale 1
end

if showProgress
    fprintf('  ch %3d/%3d: done\n', ch, nch);
end
end