function [mse, scales, info] = compute_mMSE(data, varargin)
% compute_mMSE  Multichannel Modified Multiscale Entropy (mMSE)
%               with narrowband filtskip (Kosciessa 2020) or statistical
%               coarse-graining (Costa 2002), time-resolved output,
%               and multiple coarse-graining statistics.
%
% Scale 1 is excluded: it is equivalent to uniscale SampEn and can be
% computed directly via compute_SampEn. Scales run from 2 to num_scales.
%
% Scale-wise r (Kosciessa et al. 2020, Modification 2): data are
% z-normalized once at the start, then each scale's signal (filtered or
% coarse-grained) is further normalized to unit variance before SampEn,
% so r is always relative to the scale-specific SD rather than the
% original broadband variance. This scale-wise normalization is the
% defining feature of the method and is applied once, explicitly, per
% scale. compute_SampEn is therefore called with 'ZScore', false so it
% does NOT re-normalize each segment on top of it. This matters most for
% the time-resolved output: without it, every time window would be
% renormalized to its own variance and windows would no longer be
% comparable to one another or to the whole-epoch value.
%
% Note that the different coarse-graining methods only apply to filter_mode = 'none'.
% the filt-skip method filters to isolate the frequency band and then decimates
% by selecting every s-th sample, so applying mean/std/variance operators to
% bins would change the measure's meaning and defeat the purpose of the spectral isolation.
%
% Author: Cedric Cannard, 2025
% EEGLAB ASCENT PLUGIN

% -------- Parse inputs
p = inputParser;
p.addParameter('Fs',               [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.addParameter('m',                2,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('tau',              1,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('r',               .15, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('num_scales',      30,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('coarsing',      'mean', @(s) any(strcmpi(s,{'median','mean','std','sd','standard deviation','var','variance'})));
p.addParameter('filter_mode', 'narrowband', @(s) any(strcmpi(s,{'none','narrowband'})));
p.addParameter('BlockSize',     2000, @(x) isnumeric(x) && isscalar(x) && x>=100);
p.addParameter('pdistMaxGB',     2.0, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('Parallel',       true, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress',       true, @(x) islogical(x) && isscalar(x));
p.addParameter('TimeWin',          10, @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.addParameter('TimeStep',         [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.addParameter('TimeOnly',      false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
o = p.Results;

% Validate TimeOnly
if o.TimeOnly && isempty(o.TimeWin)
    warning('compute_mMSE: TimeOnly=true but TimeWin is empty - whole-epoch MSE will be computed.');
    o.TimeOnly = false;
end

% Defaults
o.MinSamplesPerBin = max(4, o.m + 1);
o.MinWinSamples    = 100;
o.NBWiden          = 0.05;
o.PadLen           = [];
o.PlotTime         = false;

info = [];

% -------- Coerce input to [nChan x nSamples]
[dat, fs_from_struct] = coerce_data_matrix(data);
if isempty(o.Fs)
    if isempty(fs_from_struct)
        error('compute_mMSE: Fs must be provided or present in EEG.srate.');
    else
        o.Fs = fs_from_struct;
    end
end

doTime = ~isempty(o.TimeWin);
[nChan, nSamp] = size(dat);
ct = parse_coarse_type(o.coarsing);

% -------- Time grid
Tsec     = nSamp / o.Fs;
tCenters = [];
if doTime
    step = iff(isempty(o.TimeStep), o.TimeWin/2, o.TimeStep);
    t0   = o.TimeWin/2;
    t1   = Tsec - o.TimeWin/2;
    if t1 < t0
        warning('compute_mMSE: TimeWin > data duration; disabling time-resolved output.');
        doTime = false; o.TimeOnly = false;
    else
        tCenters = t0:step:t1;
    end
end
nTOI = numel(tCenters);

if doTime && nTOI > 5000
    warning(['compute_mMSE: %d time windows requested (TimeWin=%.3gs, TimeStep=%.3gs). ' ...
        'This may be very slow or crash. Consider increasing TimeWin/TimeStep ' ...
        'or setting TimeOnly=true.'], nTOI, o.TimeWin, o.TimeStep);
end

% -------- Cap number of scales
% S is the maximum scale value (inclusive). Actual scales run 2:S,
% giving S-1 values. All arrays are length S-1, indexed by (actualScale-1).
S    = max(2, floor(o.num_scales));   % enforce minimum of 2 so at least scale 2 exists
maxS = floor(nSamp / o.MinSamplesPerBin);
if S > maxS
    warning('compute_mMSE: Reducing num_scales from %d to %d (need >=%d coarse bins).', ...
        S, maxS, o.MinSamplesPerBin);
    S = maxS;
end
if S < 2
    error('compute_mMSE: not enough samples to compute even scale 2. Provide longer data.');
end

% -------- Build frequency bands (one entry per actual scale 2:S)
nyq    = o.Fs / 2;
bands  = cell(1, S-1);
scales = cell(1, S-1);
if strcmpi(o.filter_mode,'none')
    for s = 2:S
        bands{s-1}  = [NaN NaN];
        scales{s-1} = num2str(s);
    end
else
    for s = 2:S
        lo_raw = nyq / (s+1);
        hi_raw = nyq / s;
        lo_eff = max(0,            lo_raw * (1 - o.NBWiden));
        hi_eff = min(nyq*(1-1e-6), hi_raw * (1 + o.NBWiden));
        if hi_eff <= lo_eff
            bands{s-1}  = [NaN NaN];
            scales{s-1} = 'invalid';
        else
            bands{s-1}  = [lo_eff, hi_eff];
            scales{s-1} = sprintf('[%.3f %.3f] Hz', lo_eff, hi_eff);
        end
        if o.Progress && all(isfinite(bands{s-1}))
            fprintf('  Scale %2d band: %s\n', s, scales{s-1});
        end
    end
end

% -------- Pre-allocate (S-1 columns, one per actual scale 2:S)
mse      = nan(nChan, S-1);
cg_len   = zeros(1, S-1);
mse_time = [];
if doTime, mse_time = nan(nChan, S-1, nTOI); end

% -------- Progress setup
useParScales = o.Parallel && ~isempty(ver('parallel'));
hWB = []; dq = [];
if o.Progress
    toStr  = tern(o.TimeOnly, ' | TimeOnly=ON', '');
    twStr  = tern(doTime, sprintf('%.3gs', o.TimeWin), 'off');
    parStr = tern(useParScales, 'ON', 'OFF');
    fprintf('mMSE: %d scales (2:%d) | Filter=%s | Coarse=%s | Parallel=%s | TimeWin=%s | nTOI=%d%s\n', ...
        S-1, S, o.filter_mode, ct, parStr, twStr, nTOI, toStr);
    if useParScales
        dq = parallel.pool.DataQueue;
        afterEach(dq, @(si) fprintf('  scale %2d done\n', si+1));
    else
        try
            hWB = waitbar(0, 'Computing mMSE...', 'Name', 'compute_mMSE');
        catch
            hWB = [];
        end
    end
end

% -------- Fill missing values
for ch = 1:nChan
    x = dat(ch,:);
    if any(~isfinite(x))
        try
            dat(ch,:) = fillmissing(x, 'linear', 'EndValues', 'nearest');
        catch
            isn  = ~isfinite(x);
            idx  = find(~isn,1,'first'); if ~isempty(idx), x(1:idx-1)   = x(idx); end
            idx  = find(~isn,1,'last');  if ~isempty(idx), x(idx+1:end) = x(idx); end
            prev = [x(1), x(1:end-1)];
            next = [x(2:end), x(end)];
            bad  = isn & isfinite(prev) & isfinite(next);
            x(bad) = 0.5*(prev(bad)+next(bad));
            dat(ch,:) = x;
        end
    end
end

% -------- Z-normalize each channel (Kosciessa modification 2 requires this
%          as baseline so that r is interpretable as a fraction of SD)
for ch = 1:nChan
    mu = mean(dat(ch,:), 'omitnan');
    sg = std(dat(ch,:), 0, 2, 'omitnan');
    if sg > 0
        dat(ch,:) = (dat(ch,:) - mu) / sg;
    end
end

% =========================================================================
%  SCALE LOOP — PARALLEL
%  Loop variable s is the array index (1:S-1).
%  Actual scale value is actualScale = s+1.
% =========================================================================
if useParScales

    parfor s = 1:S-1
        actualScale = s + 1;
        feff        = o.Fs / actualScale;
        band        = bands{s};
        Y           = dat;
        cg_len_s    = 0;
        mse_s       = nan(nChan,1);
        mseslice    = [];

        if strcmpi(o.filter_mode,'none') || any(isnan(band))

            [CG, nBins] = coarsegrain_stat(Y, actualScale, ct);
            cg_len_s    = nBins;
            if nBins >= max([o.MinSamplesPerBin, o.m+1, o.MinWinSamples])
                % Scale-wise r: normalize CG to unit variance per channel
                CG_std = std(CG, 0, 2);
                CG_std(CG_std < 1e-10) = 1;
                CG = CG ./ CG_std;
                if ~o.TimeOnly
                    for ch = 1:nChan
                        mse_s(ch) = compute_SampEn(CG(ch,:), 'm', o.m, 'tau', o.tau, ...
                            'r', o.r, 'ZScore', false, 'BlockSize', o.BlockSize, 'pdistMaxGB', o.pdistMaxGB, ...
                            'Parallel', false, 'Progress', false);
                    end
                end
                if doTime
                    mseslice = time_windows_one_scale(CG, feff, o, nTOI, tCenters);
                end
            end

        else
            isHPonly = isinf(band(2));
            Y = fir_zero_phase(Y, band, o.Fs, o.PadLen);
            % Scale-wise r: normalize filtered signal to unit variance per
            % channel so r is relative to scale-specific SD (Kosciessa 2020)
            Y_std = std(Y, 0, 2);
            Y_std(Y_std < 1e-10) = 1;
            Y = Y ./ Y_std;

            if isHPonly
                nStarts   = 1;
                nBins_eff = size(Y,2);
            else
                nStarts   = actualScale;
                nBins_eff = floor(size(Y,2)/actualScale);
            end
            cg_len_s = nBins_eff;

            if ~isHPonly
                pbw_guard = band(2) - band(1);
                tw_guess  = max(0.15*max(pbw_guard,eps), o.Fs/size(Y,2));
                if pbw_guard < 2*tw_guess
                    if ~isempty(dq) && o.Progress, send(dq,s); end
                    cg_len(s) = cg_len_s; %#ok<PFOUS>
                    continue
                end
            end
            if nBins_eff < max(4, o.m+1)
                if ~isempty(dq) && o.Progress, send(dq,s); end
                cg_len(s) = cg_len_s; %#ok<PFOUS>
                continue
            end

            if ~o.TimeOnly
                vals = nan(nChan, nStarts);
                for is = 1:nStarts
                    CG = iff(nStarts==1, Y, Y(:, is:actualScale:is+(nBins_eff-1)*actualScale));
                    tmp = nan(nChan,1);
                    for ch = 1:nChan
                        tmp(ch) = compute_SampEn(CG(ch,:), 'm', o.m, 'tau', o.tau, ...
                            'r', o.r, 'ZScore', false, 'BlockSize', o.BlockSize, 'pdistMaxGB', o.pdistMaxGB, ...
                            'Parallel', false, 'Progress', false);
                    end
                    vals(:,is) = tmp;
                end
                mse_s = mean(vals, 2, 'omitnan');
            end

            if doTime
                feff_eff = tern(nStarts==1, o.Fs, feff);
                mseslice = time_windows_filtskip_one_scale(Y, actualScale, nStarts, nBins_eff, feff_eff, o, nTOI, tCenters);
            end
        end

        mse(:,s)  = mse_s;    %#ok<PFBNS>
        cg_len(s) = cg_len_s; %#ok<PFBNS>
        if doTime && ~isempty(mseslice)
            mse_time(:,s,:) = mseslice; %#ok<PFBNS>
        end
        if ~isempty(dq) && o.Progress, send(dq,s); end
    end

else
    % =========================================================================
    %  SCALE LOOP — SERIAL
    %  Same convention: s is array index (1:S-1), actualScale = s+1.
    % =========================================================================
    for s = 1:S-1
        actualScale = s + 1;
        feff        = o.Fs / actualScale;
        Y           = dat;
        band        = bands{s};

        if o.Progress
            fprintf('  scale %2d/%2d\n', actualScale, S);
            if ~isempty(hWB) && isvalid(hWB)
                try waitbar(s/(S-1), hWB, sprintf('Computing mMSE... (%d/%d)', actualScale, S)); catch, end
            end
        end

        if strcmpi(o.filter_mode,'none') || any(isnan(band))

            [CG, nBins] = coarsegrain_stat(Y, actualScale, ct);
            cg_len(s)   = nBins;
            if nBins >= max([o.MinSamplesPerBin, o.m+1, o.MinWinSamples])
                % Scale-wise r: normalize CG to unit variance per channel
                CG_std = std(CG, 0, 2);
                CG_std(CG_std < 1e-10) = 1;
                CG = CG ./ CG_std;
                if ~o.TimeOnly
                    for ch = 1:nChan
                        mse(ch,s) = compute_SampEn(CG(ch,:), 'm', o.m, 'tau', o.tau, ...
                            'r', o.r, 'ZScore', false, 'BlockSize', o.BlockSize, 'pdistMaxGB', o.pdistMaxGB, ...
                            'Parallel', false, 'Progress', false);
                    end
                end
                if doTime
                    mse_time(:,s,:) = time_windows_one_scale(CG, feff, o, nTOI, tCenters);
                end
            end
            continue
        end

        isHPonly = isinf(band(2));
        Y = fir_zero_phase(Y, band, o.Fs, o.PadLen);
        % Scale-wise r: normalize filtered signal to unit variance per
        % channel so r is relative to scale-specific SD (Kosciessa 2020)
        Y_std = std(Y, 0, 2);
        Y_std(Y_std < 1e-10) = 1;
        Y = Y ./ Y_std;

        if isHPonly
            nStarts   = 1;
            nBins_eff = size(Y,2);
        else
            nStarts   = actualScale;
            nBins_eff = floor(size(Y,2)/actualScale);
        end
        cg_len(s) = nBins_eff;

        if ~isHPonly
            pbw_guard = band(2) - band(1);
            tw_guess  = max(0.15*max(pbw_guard,eps), o.Fs/size(Y,2));
            if pbw_guard < 2*tw_guess
                if o.Progress
                    fprintf('  [drop] scale %d: FIR band too narrow (pbw=%.4g Hz < 2*tw=%.4g Hz)\n', ...
                        actualScale, pbw_guard, tw_guess);
                end
                continue
            end
        end
        if nBins_eff < max(4, o.m+1)
            if o.Progress
                fprintf('  [drop] scale %d: nBins=%d (<%d)\n', actualScale, nBins_eff, max(4,o.m+1));
            end
            continue
        end

        if ~o.TimeOnly
            vals = nan(nChan, nStarts);
            for is = 1:nStarts
                CG = iff(nStarts==1, Y, Y(:, is:actualScale:is+(nBins_eff-1)*actualScale));
                tmp = nan(nChan,1);
                for ch = 1:nChan
                    tmp(ch) = compute_SampEn(CG(ch,:), 'm', o.m, 'tau', o.tau, ...
                        'r', o.r, 'ZScore', false, 'BlockSize', o.BlockSize, 'pdistMaxGB', o.pdistMaxGB, ...
                        'Parallel', false, 'Progress', false);
                end
                vals(:,is) = tmp;
            end
            mse(:,s) = mean(vals, 2, 'omitnan');
        end

        if doTime
            feff_eff = tern(nStarts==1, o.Fs, feff);
            mse_time(:,s,:) = time_windows_filtskip_one_scale(Y, actualScale, nStarts, nBins_eff, feff_eff, o, nTOI, tCenters);
        end
    end
end

if ~isempty(hWB) && isvalid(hWB)
    try close(hWB); catch, end
end

if doTime
    info.mse_time = mse_time;
    info.time_sec = tCenters;
    if o.PlotTime
        scIdx = unique(round(linspace(1, S-1, min(6, S-1))));
        plot_mMSE_timecourse(mse_time, tCenters, scales, 'ScaleIdx', scIdx);
    end
end
info.TimeOnly = o.TimeOnly;
info.cg_len   = cg_len;

end

%% =========================================================================
% LOCAL HELPERS
% =========================================================================

function y = iff(cond, a, b)
if cond, y = a; else, y = b; end
end

function y = tern(cond, a, b)
if cond, y = a; else, y = b; end
end

function [X, fs] = coerce_data_matrix(data)
fs = [];
if isnumeric(data)
    X = double(data);
elseif isstruct(data) && isfield(data,'data')
    X = double(data.data);
    if isfield(data,'srate') && ~isempty(data.srate)
        fs = double(data.srate);
    end
else
    error('compute_mMSE: data must be [nChan x nSamples] numeric or an EEGLAB EEG struct.');
end
if ndims(X) ~= 2
    error('compute_mMSE: epoched/3-D data is no longer supported. Provide continuous 2-D data [nChan x nSamples].');
end
if ~isreal(X)
    warning('compute_mMSE: complex data; using real part.');
    X = real(X);
end
end

function ct = parse_coarse_type(token)
ct = 'mean';
if isempty(token), return; end
if isnumeric(token)
    switch round(token)
        case 0,    ct = 'mean';
        case 1,    ct = 'std';
        case 2,    ct = 'variance';
        otherwise, ct = 'mean';
    end
else
    t = lower(regexprep(strtrim(char(token)),'[^a-z]',''));
    if     any(strcmp(t,{'std','sd','standard deviation'})), 	ct = 'std';
    elseif any(strcmp(t,{'variance','var'})),         			ct = 'variance';
    elseif any(strcmp(t,{'mean','average'})),       			ct = 'mean';
    elseif any(strcmp(t,{'median','med'})),                    	ct = 'median';
    end
end
end

function [CG, nBins] = coarsegrain_stat(Y, s, ct)
nChan = size(Y,1);
nFull = floor(size(Y,2)/s)*s;
if nFull == 0
    CG = nan(nChan,0);
    nBins = 0;
    return
end
nBins = nFull / s;
Z = reshape(Y(:,1:nFull), nChan, s, nBins);
CG = nan(nChan, nBins);
for ch = 1:nChan
    CG(ch,:) = coarsegrain(squeeze(Z(ch,:,:)), ct);
end
end

function cg = coarsegrain(Y, ct)
% Coarse-grain down dim 1 (over the s samples of each window) -> 1 x nBins.
switch lower(strtrim(ct))
    case {'mean'}
        cg = mean(Y, 1);
    case {'median'}
        cg = median(Y, 1);
    case {'std','sd','standard deviation'}
        cg = std(Y, 0, 1);
    case {'var','variance'}
        cg = var(Y, 0, 1);
    otherwise
        cg = mean(Y, 1);
end
end

function M = time_windows_one_scale(CG, feff, o, nTOI, tCenters)
nChan = size(CG,1);
M     = nan(nChan, nTOI);
W     = max(1, round(o.TimeWin * feff));
if W < max(o.m+1, o.MinWinSamples), return; end
centers_bins = round(tCenters * feff) + 1;
nBins  = size(CG,2);
starts = centers_bins - floor(W/2);
stops  = starts + W - 1;
ki     = find(starts>=1 & stops<=nBins);
for ii = 1:numel(ki)
    seg = CG(:, starts(ki(ii)):stops(ki(ii)));
    for ch = 1:nChan
        M(ch,ki(ii)) = compute_SampEn(seg(ch,:), 'm', o.m, 'tau', o.tau, ...
            'r', o.r, 'ZScore', false, 'BlockSize', o.BlockSize, 'pdistMaxGB', o.pdistMaxGB, ...
            'Parallel', false, 'Progress', false);
    end
end
end

function M = time_windows_filtskip_one_scale(Y, s, nStarts, nBins_eff, feff_eff, o, nTOI, tCenters)
nChan  = size(Y,1);
M      = nan(nChan, nTOI);
W      = max(1, round(o.TimeWin * feff_eff));
if W < max(4, o.m+1), return; end
centers_bins = round(tCenters * feff_eff) + 1;
nAvail = tern(nStarts==1, size(Y,2), nBins_eff);
starts = centers_bins - floor(W/2);
stops  = starts + W - 1;
ki     = find(starts>=1 & stops<=nAvail);
for ii = 1:numel(ki)
    tiAbs = ki(ii);
    acc   = zeros(nChan,1);
    good  = true;
    for is = 1:nStarts
        if nStarts == 1
            idxv = starts(tiAbs):stops(tiAbs);
        else
            idx0 = is + (starts(tiAbs)-1)*s;
            idxv = idx0:s:idx0+(W-1)*s;
            idxv(idxv>size(Y,2)) = [];
        end
        if numel(idxv) < max(4,o.m+1), good = false; break; end
        seg = Y(:,idxv);
        for ch = 1:nChan
            acc(ch) = acc(ch) + compute_SampEn(seg(ch,:), 'm', o.m, 'tau', o.tau, ...
                'r', o.r, 'ZScore', false, 'BlockSize', o.BlockSize, 'pdistMaxGB', o.pdistMaxGB, ...
                'Parallel', false, 'Progress', false);
        end
    end
    if good, M(:,tiAbs) = acc ./ nStarts; end
end
end

%% =========================================================================
% FILTERS
% =========================================================================

function Y = fir_zero_phase(Y, band, fs, padLen)
maxOrder = 2000;
nyq = fs/2;
lo  = band(1); hi = band(2);
if isinf(hi)
    Wn = max(eps, lo/nyq);                         ftype = 'high';
elseif lo <= 0
    Wn = min(0.999999, hi/nyq);                    ftype = 'low';
else
    Wn = [max(eps,lo/nyq), min(0.999999,hi/nyq)]; ftype = 'bandpass';
end
T       = size(Y,2);
pbw     = iff(isinf(hi), nyq-lo, hi-lo);
transHz = max(0.5, 0.15*max(pbw, eps));
tw      = max(transHz, fs/T);
ord     = min(maxOrder, max(10, ceil(3.3*fs/tw)));
if mod(ord,2)==1, ord = ord+1; end
b   = fir1(ord, Wn, ftype, 'scale');
L   = length(b)-1;
pad = iff(isempty(padLen), min(max(100,9*L),floor((T-1)/2)), min(padLen,floor((T-1)/2)));
if pad > 0
    Yp = [fliplr(Y(:,1:pad)), Y, fliplr(Y(:,end-pad+1:end))];
    Yp = filtfilt(b, 1, Yp.').';
    Y  = Yp(:, pad+1:end-pad);
else
    Y = filtfilt(b, 1, Y.').';
end
end