function [CMFE, scales] = compute_CMFE(data, varargin)
% compute_CMFE  Composite Multiscale Fuzzy Entropy (CMFE).
%
%   [CMFE, scales] = compute_CMFE(data, 'm', 2, 'r', 0.15, 'tau', 1, ...
%                                 'n', 2, 'coarsing', 'mean', ...
%                                 'Mode', 'local', 'num_scales', 30, 'zNorm', 0)
%
% Inputs
%   data : EEGLAB EEG struct with .data OR numeric [n_ch x n_samp]
%
% Name-value parameters
%   'm'               : embedding dimension (default = 2)
%   'r'               : similarity bound, fraction of SD (default = 0.15)
%   'tau'             : embedding delay for FuzEn (default = 1)
%   'n'               : fuzzy exponent (default = 2)
%   'coarsing'        : 'mean' | 'median' | 'std' | 'var' (default = 'mean')
%   'Mode'            : 'local' | 'global' (default = 'local')
%   'num_scales'      : requested number of scales (default = 30)
%   'zNorm'           : varying-tolerance rescaling, integer 0-4 (default 0 = OFF)
%                       0  fixed tolerance across scales (classic)
%                       1  r*std / 2  r*var / 3  r*mad(mean) / 4  r*mad(median)
%                       A SINGLE tolerance is computed per scale (from the
%                       unshifted coarse-grained series) and shared across all
%                       offsets, so the composite average stays coherent.
%   'IncludeScale1'   : include scale 1 (= plain FuzzEn) (default = false).
%                       Ignored for 'std'/'var' (single-point spread is undefined).
%   'MinSamplesPerBin': minimum coarse-grained bins for a valid FE estimate (default = 4)
%   'Parallel'        : true/false (default = true)
%   'Progress'        : true/false (default = true)
%
% Outputs
%   CMFE   : [n_channels x numel(scales)] composite multiscale fuzzy entropy
%   scales : scales retained
%
% Notes
%   • CMFE averages Fuzzy Entropy values across shifted coarse-grained series.
%   • This differs from RCMFE, which averages phi_m and phi_m+1 first, then
%     computes the log-ratio.
%   • The signal is z-scored once per channel; the coarse-grained series is
%     NOT re-normalized. With zNorm=0 the tolerance r is fixed across scales;
%     with zNorm>0 it is rescaled per scale (shared across offsets).
%   • Scale 1 is skipped by default so the scale range matches the other
%     multiscale measures (MSE, MFE, mMSE). Set IncludeScale1=true to add it
%     for mean/median coarse-graining.
%   • For best agreement with the original FuzEn-style literature, use
%     'Mode','local'. Use 'global' for a SampEn-like variant.

p = inputParser;
p.addRequired('data', @(x) (isstruct(x) && isfield(x,'data')) || (isnumeric(x) && ndims(x)==2));
p.addParameter('m', 2,                 @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('r', 0.15,              @(x) isnumeric(x) && isscalar(x) && x>0 && x<2);
p.addParameter('tau', 1,               @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('n', 2,                 @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('coarsing','mean',      @(s) any(strcmpi(s,{'median','mean','std','sd','standard deviation','var','variance'})));
p.addParameter('Mode','local',         @(s) any(strcmpi(s,{'local','global'})));
p.addParameter('num_scales', 30,       @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('zNorm', 0,             @(x) isnumeric(x) && isscalar(x) && ismember(x,0:4));
p.addParameter('IncludeScale1', false, @(x) islogical(x) && isscalar(x));
p.addParameter('MinSamplesPerBin', 4,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('Parallel', true,       @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,       @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
r            = p.Results.r;
tau          = p.Results.tau;
n_exp        = p.Results.n;
coarseType   = p.Results.coarsing;
modeType     = lower(p.Results.Mode);
nScales_req  = p.Results.num_scales;
zNorm        = p.Results.zNorm;
incScale1    = p.Results.IncludeScale1;
minBinsHard  = p.Results.MinSamplesPerBin;
parallelMode = p.Results.Parallel;
showProg     = p.Results.Progress;

if isstruct(data), X = double(data.data); else, X = double(data); end
if size(X,1) > size(X,2), X = X.'; end
[nch, nSamp] = size(X);

S    = max(1, floor(nScales_req));
maxS = floor(nSamp / max(1, minBinsHard));
if S > maxS
    warning('compute_CMFE:ReducingScales', ...
        'Reducing num_scales from %d to %d to keep >=%d coarse bins at every scale.', ...
        S, maxS, minBinsHard);
    S = maxS;
end
if S < 1, S = 1; end

% Scale range: skip scale 1 by default (matches MSE/MFE/mMSE). Include it only
% when explicitly requested and only for mean/median coarse-graining.
isSpread = ismember(lower(strtrim(coarseType)), {'sd','std','standard deviation','var','variance'});
if incScale1 && ~isSpread
    scales = 1:S;
else
    if incScale1 && isSpread
        warning('compute_CMFE: IncludeScale1 ignored for spread coarse-graining (single-point spread is undefined).');
    end
    scales = 2:S;
end
if isempty(scales)
    error('compute_CMFE: not enough samples/scales to compute even scale 2. Provide longer data or increase num_scales.');
end

Xz = zscore_channels_local(X);
CMFE = nan(nch, numel(scales));

if showProg
    state = ternary(parallelMode && ~isempty(ver('parallel')), 'on', 'off');
    fprintf('CMFE: %d ch | m=%g, tau=%g, r=%g, n=%g | coarse=%s | mode=%s | zNorm=%d | scales=%d:%d | parallel=%s\n', ...
        nch, m, tau, r, n_exp, upperLabel(coarseType), upper(modeType), zNorm, scales(1), scales(end), state);
end

useWB = ~parallelMode && usejava('desktop') && showProg;
hWB = [];
if useWB
    try, hWB = waitbar(0, 'Computing CMFE...', 'Name', 'compute_CMFE'); catch, hWB = []; end
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
        v = nan(1, numel(scales));

        for ii = 1:numel(scales)
            s   = scales(ii);
            r_s = scale_tolerance(sig, s, coarseType, r, zNorm, minBinsHard, m);

            fe_vals = nan(1, s);
            for off = 1:s
                xoff = sig(off:end);
                Loff = floor(numel(xoff) / s) * s;
                nBins_off = Loff / s;
                if nBins_off < max(minBinsHard, m+1), continue; end

                Y  = reshape(xoff(1:Loff), s, []);
                cg = coarsegrain(Y, coarseType);

                [fe, ~, ~] = fuzz_engine_raw(cg, m, r_s, n_exp, tau, ...
                    'exponential', false, 2000, 2.0, modeType);
                fe_vals(off) = double(fe);
            end

            valid = isfinite(fe_vals);
            if any(valid), v(ii) = mean(fe_vals(valid)); end
        end

        CMFE(ch,:) = v;
        if useDQ, send(dq,1); end
    end
else
    for ch = 1:nch
        sig = Xz(ch,:);
        v = nan(1, numel(scales));

        for ii = 1:numel(scales)
            s   = scales(ii);
            r_s = scale_tolerance(sig, s, coarseType, r, zNorm, minBinsHard, m);

            fe_vals = nan(1, s);
            for off = 1:s
                xoff = sig(off:end);
                Loff = floor(numel(xoff) / s) * s;
                nBins_off = Loff / s;
                if nBins_off < max(minBinsHard, m+1), continue; end

                Y  = reshape(xoff(1:Loff), s, []);
                cg = coarsegrain(Y, coarseType);

                [fe, ~, ~] = fuzz_engine_raw(cg, m, r_s, n_exp, tau, ...
                    'exponential', false, 2000, 2.0, modeType);
                fe_vals(off) = double(fe);
            end

            valid = isfinite(fe_vals);
            if any(valid), v(ii) = mean(fe_vals(valid)); end
        end

        CMFE(ch,:) = v;
        if ~isempty(hWB) && isvalid(hWB)
            try, waitbar(ch/nch, hWB, sprintf('Computing CMFE... (%d/%d)', ch, nch)); catch, end
        end
    end
    if ~isempty(hWB) && isvalid(hWB), try, close(hWB); catch, end, end
end

    function notifyProgress(~)
        nDone = nDone + 1;
        step = max(1, round(0.05*nch));
        if nDone == 1 || nDone == nch || mod(nDone, step) == 0
            fprintf('  progress: ch %d/%d\n', nDone, nch);
        end
    end
end

%% Local helpers

function Xz = zscore_channels_local(X)
Xz = X;
for c = 1:size(X,1)
    x = X(c,:);
    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');
    if ~isfinite(sd) || sd == 0
        Xz(c,:) = 0;
    else
        Xz(c,:) = (x - mu) ./ sd;
    end
end
end


function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end