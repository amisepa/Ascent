function [CMFE, scales] = compute_CMFE(data, varargin)
% compute_CMFE  Composite Multiscale Fuzzy Entropy (CMFE).

p = inputParser;
p.addRequired('data', @(x) (isstruct(x) && isfield(x,'data')) || (isnumeric(x) && ndims(x)==2));
p.addParameter('m', 2,                 @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('r', 0.15,              @(x) isnumeric(x) && isscalar(x) && x>0 && x<2);
p.addParameter('tau', 1,               @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('n', 2,                 @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('coarsing','std',       @(s) any(strcmpi(s,{'median','mean','trimmed mean','trimmed','tmean','trim20','std','sd','standard deviation','var','variance'})));
p.addParameter('num_scales', 15,       @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('MinSamplesPerBin', 4,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('Parallel', true,       @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,       @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
r            = p.Results.r;
tau          = p.Results.tau;
n_exp        = p.Results.n;
coarseType   = p.Results.coarsing;
nScales_req  = p.Results.num_scales;
minBinsHard  = p.Results.MinSamplesPerBin;
parallelMode = p.Results.Parallel;
showProg     = p.Results.Progress;

if isstruct(data)
    X = double(data.data);
else
    X = double(data);
end
if size(X,1) > size(X,2)
    X = X.';
end
[nch, nSamp] = size(X);

S    = max(1, floor(nScales_req));
maxS = floor(nSamp / max(1, minBinsHard));
if S > maxS
    warning('compute_CMFE:ReducingScales', ...
        'Reducing num_scales from %d to %d to keep >=%d coarse bins at every scale.', ...
        S, maxS, minBinsHard);
    S = maxS;
end
if S < 1
    S = 1;
end
scales = 1:S;

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

cl = lower(strtrim(coarseType));
if any(strcmp(cl, {'sd','std','standard deviation'}))
    coarseLabel = 'STD';
elseif any(strcmp(cl, {'var','variance'}))
    coarseLabel = 'VAR';
elseif strcmp(cl, 'mean')
    coarseLabel = 'MEAN';
elseif strcmp(cl, 'median')
    coarseLabel = 'MEDIAN';
elseif any(strcmp(cl, {'trimmed mean','trimmed','tmean','trim20'}))
    coarseLabel = 'TRIM20';
else
    coarseLabel = upper(coarseType);
end

CMFE = nan(nch, S);
if showProg
    if parallelMode && ~isempty(ver('parallel'))
        fprintf('CMFE: %d ch | m=%g, tau=%g, r=%g, n=%g | coarse=%s | S=%d | parallel=on\n', ...
            nch, m, tau, r, n_exp, coarseLabel, S);
    else
        fprintf('CMFE: %d ch | m=%g, tau=%g, r=%g, n=%g | coarse=%s | S=%d | parallel=off\n', ...
            nch, m, tau, r, n_exp, coarseLabel, S);
    end
end

useWB = ~parallelMode && usejava('desktop') && showProg;
hWB = [];
if useWB
    try
        hWB = waitbar(0, 'Computing CMFE...', 'Name', 'compute_CMFE');
    catch
        hWB = [];
    end
end

useDQ = parallelMode && ~isempty(ver('parallel')) && showProg;
if useDQ
    dq = parallel.pool.DataQueue;
    nDone = 0;
    afterEach(dq, @notifyProgress);
end

if parallelMode && ~isempty(ver('parallel'))
    parfor ch = 1:nch
        sig = Xz(ch,:);
        v = nan(1, S);

        for s = 1:S
            if s == 1
                [fe, ~, ~] = fuzz_engine_raw(sig, m, r, n_exp, tau, 'exponential', false);
                v(1) = double(fe);
                continue
            end

            nBins = floor(numel(sig) / s);
            if nBins < max(minBinsHard, m+1)
                continue
            end

            fe_vals = nan(1, s);
            for off = 1:s
                xoff = sig(off:end);
                Loff = floor(numel(xoff) / s) * s;
                nBins_off = Loff / s;

                if nBins_off < max(minBinsHard, m+1)
                    continue
                end

                Y  = reshape(xoff(1:Loff), s, []);
                cg = coarsegrain(Y, coarseType);

                [fe, ~, ~] = fuzz_engine_raw(cg, m, r, n_exp, tau, 'exponential', false);
                fe_vals(off) = double(fe);
            end

            valid = isfinite(fe_vals);
            if any(valid)
                v(s) = mean(fe_vals(valid));
            end
        end

        CMFE(ch,:) = v;
        if useDQ
            send(dq, 1);
        end
    end
else
    for ch = 1:nch
        sig = Xz(ch,:);
        v = nan(1, S);

        for s = 1:S
            if s == 1
                [fe, ~, ~] = fuzz_engine_raw(sig, m, r, n_exp, tau, 'exponential', false);
                v(1) = double(fe);
                continue
            end

            nBins = floor(numel(sig) / s);
            if nBins < max(minBinsHard, m+1)
                continue
            end

            fe_vals = nan(1, s);
            for off = 1:s
                xoff = sig(off:end);
                Loff = floor(numel(xoff) / s) * s;
                nBins_off = Loff / s;

                if nBins_off < max(minBinsHard, m+1)
                    continue
                end

                Y  = reshape(xoff(1:Loff), s, []);
                cg = coarsegrain(Y, coarseType);

                [fe, ~, ~] = fuzz_engine_raw(cg, m, r, n_exp, tau, 'exponential', false);
                fe_vals(off) = double(fe);
            end

            valid = isfinite(fe_vals);
            if any(valid)
                v(s) = mean(fe_vals(valid));
            end
        end

        CMFE(ch,:) = v;

        if ~isempty(hWB) && isvalid(hWB)
            try
                waitbar(ch/nch, hWB, sprintf('Computing CMFE... (%d/%d)', ch, nch));
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
        step = max(1, round(0.05*nch));
        if nDone == 1 || nDone == nch || mod(nDone, step) == 0
            fprintf('  progress: ch %d/%d\n', nDone, nch);
        end
    end
end