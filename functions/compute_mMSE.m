function [mse, scales, info] = compute_mMSE(data, varargin)
% compute_mMSE  Multichannel Modified Multiscale Entropy (mMSE)
%               with low-pass or band-pass filt-skip (Kosciessa 2020) or
%               statistical coarse-graining (Costa 2002), time-resolved
%               output, and multiple coarse-graining statistics.
%
% Scale 1 is excluded: it is equivalent to uniscale SampEn and can be
% computed directly via compute_SampEn. Scales run from 2 to num_scales.
%
% Scale-wise r (Kosciessa et al. 2020): data are
% z-normalized once at the start, then each scale's signal (filtered or
% coarse-grained) is further normalized to unit variance before SampEn,
% so r is always relative to the scale-specific SD rather than the
% original broadband variance. It is applied once, explicitly, per
% scale. compute_SampEn is therefore called with 'ZScore', false so it
% does NOT re-normalize each segment on top of it. This matters most for
% the time-resolved output: without it, every time window would be
% renormalized to its own variance and windows would no longer be
% comparable to one another or to the whole-epoch value.
%
% filter_mode (scale s, Nyquist frequency f_Nyq = Fs/2):
%   'lowpass'    (default) Kosciessa et al.'s low-pass version: low-pass at f_Nyq/s,
%                then keep every s-th sample from each of the s starting
%                points (filt-skip). Scale s keeps all activity below
%                f_Nyq/s, so slow, high-power activity (1/f slope, alpha)
%                is present at every scale and dominates the fine scales.
%   'narrowband' band-pass over [f_Nyq/(s+1), f_Nyq/s] (edges widened by 5%),
%                then filt-skip. Kosciessa et al. treat band-pass entropy as a
%                parallel tool, not as a replacement: it removes slow, high-power
%                activity from the fine scales, but loses sensitivity to the
%                broadband slope. 'bandpass' is accepted as an alias. The band
%                reaches above the Nyquist frequency left after skipping
%                (f_Nyq/s), so the top of each band folds back into it: about
%                12% of the band at scale 2 and 38% at scale 30 (256 Hz).
%   'none'       Costa coarse-graining (mean, median, SD or variance of
%                non-overlapping windows), one starting point. The
%                coarse-graining options apply only to this mode: filt-skip
%                filters to isolate the frequency range and then keeps every
%                s-th sample, so applying mean/std/variance operators to bins
%                would change the measure's meaning.
%
% FilterType ('lowpass' only):
%   'fir'    (default) zero-phase windowed-sinc FIR (fir1, Hamming), -6 dB
%            per pass at f_Nyq/s, transition max(0.5 Hz, 15% of the cutoff).
%   'butter' zero-phase 6th-order Butterworth, as Kosciessa et al. (applied
%            in second-order sections, stable at low cutoffs).
%   Band-pass always uses the FIR (Kosciessa et al. use Chebyshev Type I
%   filters and a passband of the cutoff +-5% for their band-pass check).
%
% StartPoints (filt-skip modes):
%   'pooled'  (default) match counts are summed over the s starting points
%             before the logarithm, as in Kosciessa et al. (equivalent to
%             refined composite MSE). Counts are normalized by the number of
%             pairs of each dimension, as in compute_SampEn.
%   'average' SampEn of each starting point, then averaged (composite MSE).
%
% The output info struct holds filter_mode, FilterType, StartPoints, Fs and
% bands ([S-1 x 2] Hz; [0 f_Nyq/s] for low-pass, NaN for 'none') so that
% plots can label each scale with the frequencies it keeps.
%
% References:
%   Kosciessa, J. Q., Kloosterman, N. A., & Garrett, D. D. (2020). Standard multiscale
%     entropy reflects neural dynamics at mismatched temporal scales: What's signal
%     irregularity got to do with it? PLoS Comput Biol, 16, e1007885.
%   Kloosterman, N. A., Kosciessa, J. Q., Lindenberger, U., Fahrenfort, J. J., &
%     Garrett, D. D. (2020). Boosts in brain signal variability track liberal shifts in
%     decision bias. eLife, 9, e54201.
%   Costa, M., Goldberger, A. L., & Peng, C.-K. (2002). Multiscale entropy analysis of
%     complex physiologic time series. Phys Rev Lett, 89, 068102.
%   Wu, S.-D., Wu, C.-W., Lin, S.-G., Lee, K.-Y., & Peng, C.-K. (2014). Analysis of
%     complex time series using refined composite multiscale entropy. Phys Lett A, 378,
%     1369-1374.
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).
% EEGLAB ASCENT PLUGIN

% -------- Parse inputs
p = inputParser;
p.addParameter('Fs',               [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.addParameter('m',                2,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('tau',              1,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('r',               .15, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('num_scales',      30,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('coarsing',      'mean', @(s) any(strcmpi(s,{'median','mean','std','sd','standard deviation','var','variance'})));
p.addParameter('filter_mode', 'lowpass', @(s) any(strcmpi(s,{'lowpass','narrowband','bandpass','none'})));
p.addParameter('FilterType',     'fir', @(s) any(strcmpi(s,{'fir','butter'})));
p.addParameter('StartPoints', 'pooled', @(s) any(strcmpi(s,{'pooled','average'})));
p.addParameter('BlockSize',     2000, @(x) isnumeric(x) && isscalar(x) && x>=100);
p.addParameter('pdistMaxGB',     2.0, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('Parallel',       true, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress',       true, @(x) islogical(x) && isscalar(x));
p.addParameter('TimeWin',          10, @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.addParameter('TimeStep',         [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.addParameter('TimeOnly',      false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
o = p.Results;
o.filter_mode = lower(o.filter_mode);
if strcmp(o.filter_mode, 'bandpass'), o.filter_mode = 'narrowband'; end
o.FilterType  = lower(o.FilterType);
o.StartPoints = lower(o.StartPoints);
if strcmp(o.FilterType, 'butter') && ~strcmp(o.filter_mode, 'lowpass')
    warning('compute_mMSE: FilterType ''butter'' applies to the low-pass mode only; using the FIR.');
    o.FilterType = 'fir';
end

% Validate TimeOnly
if o.TimeOnly && isempty(o.TimeWin)
    warning('compute_mMSE: TimeOnly=true but TimeWin is empty - whole-epoch MSE will be computed.');
    o.TimeOnly = false;
end

% Defaults
o.MinSamplesPerBin = max(4, o.m + 1);
o.MinWinSamples    = 100;
o.NBWiden          = 0.05;
o.ButterOrder      = 6;      % Kosciessa et al. 2020
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

% -------- Frequency range kept at each scale (one entry per actual scale 2:S)
nyq    = o.Fs / 2;
bands  = cell(1, S-1);
scales = cell(1, S-1);
for s = 2:S
    switch o.filter_mode
        case 'none'
            bands{s-1}  = [NaN NaN];
            scales{s-1} = num2str(s);
        case 'lowpass'
            bands{s-1}  = [0, nyq / s];
            scales{s-1} = sprintf('[0 %.3f] Hz', nyq / s);
        case 'narrowband'
            lo_eff = max(0,            nyq / (s+1) * (1 - o.NBWiden));
            hi_eff = min(nyq*(1-1e-6), nyq / s     * (1 + o.NBWiden));
            if hi_eff <= lo_eff
                bands{s-1}  = [NaN NaN];
                scales{s-1} = 'invalid';
            else
                bands{s-1}  = [lo_eff, hi_eff];
                scales{s-1} = sprintf('[%.3f %.3f] Hz', lo_eff, hi_eff);
            end
    end
    if o.Progress && ~strcmp(o.filter_mode,'none') && all(isfinite(bands{s-1}))
        fprintf('  Scale %2d band: %s\n', s, scales{s-1});
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
    if strcmp(o.filter_mode, 'none')
        fltStr = sprintf('none | Coarse=%s', ct);
    else
        fltStr = sprintf('%s (%s) | StartPoints=%s', o.filter_mode, ...
            tern(strcmp(o.filter_mode,'lowpass'), o.FilterType, 'fir'), o.StartPoints);
    end
    fprintf('mMSE: %d scales (2:%d) | Filter=%s | Parallel=%s | TimeWin=%s | nTOI=%d%s\n', ...
        S-1, S, fltStr, parStr, twStr, nTOI, toStr);
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
        x(~isfinite(x)) = NaN;   % fillmissing leaves Inf in place
        try
            dat(ch,:) = fillmissing(x, 'linear', 'EndValues', 'nearest');
        catch
            isn  = isnan(x);
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

% -------- Z-normalize each channel, so that r is a fraction of the SD
for ch = 1:nChan
    mu = mean(dat(ch,:), 'omitnan');
    sg = std(dat(ch,:), 0, 2, 'omitnan');
    if sg > 0
        dat(ch,:) = (dat(ch,:) - mu) / sg;
    end
end

% A channel with no usable samples gives NaN and does not stop the others
% (the filters need finite input)
badChan = any(~isfinite(dat), 2);
dat(badChan,:) = 0;

% =========================================================================
%  SCALE LOOP. Loop variable s is the array index (1:S-1); the actual
%  scale value is s+1.
% =========================================================================
if useParScales
    parfor s = 1:S-1
        [mse_s, cg_len_s, mseslice, note] = one_scale(dat, s+1, bands{s}, o, ct, doTime, nTOI, tCenters);
        mse(:,s)  = mse_s;
        cg_len(s) = cg_len_s;
        if doTime && ~isempty(mseslice)
            mse_time(:,s,:) = mseslice;
        end
        if o.Progress && ~isempty(note), fprintf('  [drop] scale %d: %s\n', s+1, note); end
        if ~isempty(dq) && o.Progress, send(dq,s); end
    end
else
    for s = 1:S-1
        if o.Progress
            fprintf('  scale %2d/%2d\n', s+1, S);
            if ~isempty(hWB) && isvalid(hWB)
                try waitbar(s/(S-1), hWB, sprintf('Computing mMSE... (%d/%d)', s+1, S)); catch, end
            end
        end
        [mse(:,s), cg_len(s), mseslice, note] = one_scale(dat, s+1, bands{s}, o, ct, doTime, nTOI, tCenters);
        if doTime && ~isempty(mseslice)
            mse_time(:,s,:) = mseslice;
        end
        if o.Progress && ~isempty(note), fprintf('  [drop] scale %d: %s\n', s+1, note); end
    end
end

if ~isempty(hWB) && isvalid(hWB)
    try close(hWB); catch, end
end

mse(badChan,:) = NaN;
if doTime
    mse_time(badChan,:,:) = NaN;
    info.mse_time = mse_time;
    info.time_sec = tCenters;
    if o.PlotTime
        scIdx = unique(round(linspace(1, S-1, min(6, S-1))));
        plot_mMSE_timecourse(mse_time, tCenters, scales, 'ScaleIdx', scIdx);
    end
end
info.TimeOnly    = o.TimeOnly;
info.cg_len      = cg_len;
info.filter_mode = o.filter_mode;
info.FilterType  = tern(strcmp(o.filter_mode,'lowpass'), o.FilterType, tern(strcmp(o.filter_mode,'none'), '', 'fir'));
info.StartPoints = tern(strcmp(o.filter_mode,'none'), '', o.StartPoints);
info.coarsing    = tern(strcmp(o.filter_mode,'none'), ct, '');
info.Fs          = o.Fs;
info.bands       = vertcat(bands{:});

end

%% =========================================================================
% LOCAL HELPERS
% =========================================================================

function [mse_s, cg_len_s, mseslice, note] = one_scale(dat, s, band, o, ct, doTime, nTOI, tCenters)
% Entropy of every channel at one scale s (whole epoch and, if requested,
% time-resolved). note is non-empty when the scale is skipped.
nChan    = size(dat,1);
mse_s    = nan(nChan,1);
cg_len_s = 0;
mseslice = [];
note     = '';
feff     = o.Fs / s;

% ---- Costa coarse-graining, one starting point
if strcmp(o.filter_mode,'none') || any(isnan(band))
    [CG, nBins] = coarsegrain_stat(dat, s, ct);
    cg_len_s    = nBins;
    if nBins >= max([o.MinSamplesPerBin, o.m+1, o.MinWinSamples])
        % Scale-wise r: normalize CG to unit variance per channel
        CG_std = std(CG, 0, 2);
        CG_std(CG_std < 1e-10) = 1;
        CG = CG ./ CG_std;
        if ~o.TimeOnly
            for ch = 1:nChan
                mse_s(ch) = sampen(CG(ch,:), o);
            end
        end
        if doTime
            mseslice = time_windows_one_scale(CG, feff, o, nTOI, tCenters);
        end
    else
        note = sprintf('nBins=%d (<%d)', nBins, max([o.MinSamplesPerBin, o.m+1, o.MinWinSamples]));
    end
    return
end

% ---- Filt-skip: filter, normalize, keep every s-th sample from s starts
if size(dat,2) < 2*(o.m + 2)
    note = sprintf('%d samples are too few to filter', size(dat,2));
    return
end
if strcmp(o.FilterType,'butter')
    Y = butter_zero_phase(dat, band(2), o.Fs, o.PadLen, o.ButterOrder);
else
    Y = fir_zero_phase(dat, band, o.Fs, o.PadLen);
end
% Scale-wise r: normalize filtered signal to unit variance per
% channel so r is relative to scale-specific SD (Kosciessa 2020)
Y_std = std(Y, 0, 2);
Y_std(Y_std < 1e-10) = 1;
Y = Y ./ Y_std;

nBins_eff = floor(size(Y,2)/s);
cg_len_s  = nBins_eff;

pbw_guard = band(2) - band(1);
tw_guess  = max(0.15*max(pbw_guard,eps), o.Fs/size(Y,2));
if pbw_guard < 2*tw_guess
    note = sprintf('FIR band too narrow (pbw=%.4g Hz < 2*tw=%.4g Hz)', pbw_guard, tw_guess);
    return
end
if nBins_eff < max(4, o.m+1)
    note = sprintf('nBins=%d (<%d)', nBins_eff, max(4,o.m+1));
    return
end

if ~o.TimeOnly
    vals = nan(nChan, s);
    cnt  = zeros(nChan, 4);
    for is = 1:s
        CG = Y(:, is:s:is+(nBins_eff-1)*s);
        for ch = 1:nChan
            [vals(ch,is), c] = sampen(CG(ch,:), o);
            cnt(ch,:) = cnt(ch,:) + c;
        end
    end
    if strcmp(o.StartPoints,'pooled')
        mse_s = pooled_sampen(cnt);
    else
        mse_s = mean(vals, 2, 'omitnan');
    end
end

if doTime
    mseslice = time_windows_filtskip_one_scale(Y, s, nBins_eff, feff, o, nTOI, tCenters);
end
end

function [se, cnt] = sampen(x, o)
% SampEn of an already normalized series (r is an absolute bound here)
[se, cnt] = compute_SampEn(x, 'm', o.m, 'tau', o.tau, 'r', o.r, 'ZScore', false, ...
    'BlockSize', o.BlockSize, 'pdistMaxGB', o.pdistMaxGB, 'Parallel', false, 'Progress', false);
end

function se = pooled_sampen(cnt)
% cnt: [nChan x 4] summed [matches(m+1) pairs(m+1) matches(m) pairs(m)]
A  = cnt(:,1) ./ cnt(:,2);
B  = cnt(:,3) ./ cnt(:,4);
se = -log(A ./ B);
se(~(A > 0 & B > 0)) = NaN;
end

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
    t = regexprep(lower(strtrim(char(token))),'[^a-z]','');
    if     any(strcmp(t,{'std','sd','standarddeviation'})), 	ct = 'std';
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
    CG(ch,:) = coarsegrain(reshape(Z(ch,:,:), s, nBins), ct);
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
starts = min(centers_bins - floor(W/2), nBins - W + 1);   % rounding can pass the end by one bin
stops  = starts + W - 1;
ki     = find(starts>=1 & stops<=nBins);
for ii = 1:numel(ki)
    seg = CG(:, starts(ki(ii)):stops(ki(ii)));
    for ch = 1:nChan
        M(ch,ki(ii)) = sampen(seg(ch,:), o);
    end
end
end

function M = time_windows_filtskip_one_scale(Y, s, nBins_eff, feff, o, nTOI, tCenters)
% Time-resolved filt-skip entropy: each window keeps the samples of all s
% starting points, pooled or averaged as for the whole epoch.
nChan  = size(Y,1);
M      = nan(nChan, nTOI);
W      = max(1, round(o.TimeWin * feff));
if W < max(4, o.m+1), return; end
centers_bins = round(tCenters * feff) + 1;
starts = min(centers_bins - floor(W/2), nBins_eff - W + 1);   % rounding can pass the end by one bin
stops  = starts + W - 1;
ki     = find(starts>=1 & stops<=nBins_eff);
pooled = strcmp(o.StartPoints,'pooled');
for ii = 1:numel(ki)
    tiAbs = ki(ii);
    acc   = zeros(nChan,1);
    cnt   = zeros(nChan,4);
    good  = true;
    for is = 1:s
        idx0 = is + (starts(tiAbs)-1)*s;
        idxv = idx0:s:idx0+(W-1)*s;
        idxv(idxv>size(Y,2)) = [];
        if numel(idxv) < max(4,o.m+1), good = false; break; end
        seg = Y(:,idxv);
        for ch = 1:nChan
            [se, c] = sampen(seg(ch,:), o);
            acc(ch)   = acc(ch) + se;
            cnt(ch,:) = cnt(ch,:) + c;
        end
    end
    if good
        if pooled, M(:,tiAbs) = pooled_sampen(cnt);
        else,      M(:,tiAbs) = acc ./ s;
        end
    end
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
% On short data the transition widens: filtfilt needs more than 3*ord samples,
% and mirror padding gives at most 2*T
maxOrd = floor((2*T - 4)/3);
if ord > maxOrd, ord = maxOrd - mod(maxOrd,2); end
b   = fir1(ord, Wn, ftype, 'scale');
L   = length(b)-1;
pad = iff(isempty(padLen), min(max(100,9*L),floor((T-1)/2)), min(padLen,floor((T-1)/2)));
Y   = filtfilt_mirror(Y, pad, @(X) filtfilt(b, 1, X));
end

function Y = butter_zero_phase(Y, hi, fs, padLen, order)
% Zero-phase Butterworth low-pass in second-order sections (the transfer
% function form loses precision at low normalized cutoffs).
[z, p, k] = butter(order, min(0.999999, hi/(fs/2)), 'low');
[sos, g]  = zp2sos(z, p, k);
T   = size(Y,2);
pad = iff(isempty(padLen), min(max(100, 3*round(fs/hi)*order), floor((T-1)/2)), min(padLen, floor((T-1)/2)));
Y   = filtfilt_mirror(Y, pad, @(X) filtfilt(sos, g, X));
end

function Y = filtfilt_mirror(Y, pad, ff)
% Y: [nChan x T]; ff filters down the columns of a [T x nChan] matrix
if pad > 0
    Yp = [fliplr(Y(:,1:pad)), Y, fliplr(Y(:,end-pad+1:end))];
    Yp = ff(Yp.').';
    Y  = Yp(:, pad+1:end-pad);
else
    Y = ff(Y.').';
end
end
