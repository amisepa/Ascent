function [mse, scales, info] = compute_mMSE(data, varargin)
% compute_mMSE  Multichannel Modified Multiscale Entropy (mMSE) with
%               narrowband filtskip (Kosciessa 2020) or statistical
%               coarse-graining (Costa 2002), time-resolved output,
%               and multiple coarse-graining statistics.
%
% SUMMARY
% -------
% • Two processing modes
%   (A) Statistical coarse-graining (filter_mode='none'): Costa-like block
%       reduction using the chosen statistic ('mean' | 'std' | 'variance' | 'median').
%   (B) Spectral coarse-graining (filter_mode='narrowband'): Kosciessa-style
%       filtskip — filter to a scale-specific annular band, decimate by
%       scale s, and average SampEn across all s start phases.
%
% Narrowband behavior (Kosciessa 2020):
%   • Scale 1: high-pass only at ~Nyq/2 (widened downward), NO decimation.
%   • Scales >1: bandpass annuli [Nyq/(s+1), Nyq/s] widened by ±NBWiden,
%     then skip every s samples and average across all s start phases.
%   • Auto-switch to IIR for very narrow bands; else zero-phase FIR.
%
% Scales are auto-dropped when infeasible (too few coarse bins, or FIR
% band unrealizable). Scale 1 with filter_mode='none' and coarsing='std'
% is also skipped (std of a single-element window is trivially zero).
%
%
% DIFFERENCES relative to CMSE/RCMSE (Azami)
% ------------------------------------------
% • Filtskip mode here is composite across phases by averaging entropy values
%   from all s start phases (CMSE-like). RCMSE pools the match counts
%   (phi_m, phi_m+1) across phases before taking -ln, further reducing
%   variance. This function does NOT implement that refined pooling.
% • Tolerance r is forwarded to SampEn as a fraction of the SD of each
%   coarse-grained or filtered series.
%
%
% INPUTS (name-value)
% -------------------
%   'Fs'              : sampling rate in Hz (required unless EEG.srate present)
%   'm'               : embedding dimension (default 2)
%   'tau'             : lag (default 1)
%   'r'               : tolerance fraction of SD for SampEn (default 0.15)
%   'num_scales'      : requested number of scales (default 20)
%   'coarsing'        : 'mean' | 'median' | 'std' | 'variance' (default 'std')
%   'filter_mode'     : 'none' | 'narrowband' (default 'narrowband')
%                       'none'       — Costa (2002) statistical coarse-graining
%                       'narrowband' — Kosciessa (2020) annular bandpass filtskip
%   'PadLen'          : reflection padding samples for filtfilt (default [], auto)
%   'MinSamplesPerBin': minimum coarse bins per scale (default max(4, m+1))
%   'NBWiden'         : fractional band widening for filtskip annuli
%                       (default 0.05; range [0, 0.25))
%   'Parallel'        : parallelize across scales (default true)
%   'Progress'        : console/waitbar feedback (default true)
%   'TimeWin'         : sliding window length in seconds for time-resolved
%                       output (default 10; [] to disable)
%   'TOI'             : vector of window center times in seconds; if empty,
%                       auto grid with step = TimeWin/2 (default [])
%   'TimeStep'        : step between window centers in seconds
%                       (default TimeWin/2 when TimeWin is set)
%   'TimeOnly'        : skip whole-epoch SampEn, only compute time-resolved
%                       output — faster when time course is all you need
%                       (default false; requires TimeWin to be set)
%   'PlotTime'        : auto-plot 6 evenly spaced scales when TimeWin is
%                       set (default true)
%
% Accepted data:
%   data : [nChan x nSamples] numeric OR EEGLAB EEG struct (.data, .srate)
%          Only continuous / 2-D data are supported.

% OUTPUTS
% -------
%   mse    : [nChan x nScales] whole-epoch SampEn (NaN when TimeOnly=true)
%   scales : 1 x nScales cellstr — integer labels ('none') or Hz band labels
%   info   : struct with fields:
%            .mse_time  [nChan x nScales x nTime]  (if TimeWin set)
%            .time_sec  [1 x nTime]                window centers in seconds
%            .TimeOnly  logical
%            .cg_len    [1 x nScales]              coarse-grained length per scale
%
% USAGE EXAMPLES
% --------------
% 1) Narrowband mMSE (default)
%    [mse, scales] = compute_mMSE(EEG, 'num_scales', 20);
%
% 2) Statistical path, mean coarse-graining (Costa-style)
%    [mse, scales] = compute_mMSE(EEG, 'filter_mode','none', 'coarsing','mean');
%
% 3) Time-resolved, 10-s windows every 5 s (auto-plots 6 scales)
%    [mse, scales, info] = compute_mMSE(EEG, 'TimeWin', 10, 'TimeStep', 5);
%
% 4) Time-course only (skip whole-epoch, faster)
%    [~, scales, info] = compute_mMSE(EEG, 'TimeWin', 10, 'TimeOnly', true);
%
% 5) Manual plot with custom scale selection
%    [~, scales, info] = compute_mMSE(EEG, 'TimeWin', 10, 'PlotTime', false);
%    plot_mMSE_timecourse(info.mse_time, info.time_sec, scales, 'ScaleIdx', [1 3 6 9 12 15]);
%
%
% Notes
% -----
% • Scale 1 with filter_mode='none' and coarsing='std' is skipped because
%   the std of a single-element window is always zero.
% • In narrowband mode, scale 1 is HP-only (no decimation). Scales >1 use
%   annular bandpass with filtskip. Very narrow bands auto-switch to IIR.
% • TimeOnly=true has no effect unless TimeWin is also set.
%
%
% References
% ----------
%   Costa M, Goldberger AL, Peng CK (2002). Multiscale entropy analysis of
%       complex physiologic time series. PRL, 89(6), 068102.
%
%   Kosciessa JQ et al. (2020). Standard multiscale entropy reflects neural
%       dynamics at mismatched temporal scales. PLoS Comput Biol, 16, e1007885.
%
%   Grandy TH et al. (2016). On the estimation of brain signal entropy from
%       sparse neuroimaging data. Sci Rep, 6:23073.
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Ascent plugin — Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------


% -------- Parse inputs
p = inputParser;
p.addParameter('Fs',               [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.addParameter('m',                2,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('tau',              1,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('r',               .15, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('num_scales',      20,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('coarsing',       'std');
p.addParameter('filter_mode', 'narrowband', @(s) any(strcmpi(s,{'none','narrowband'})));
p.addParameter('Parallel',          true, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress',          true, @(x) islogical(x) && isscalar(x));
p.addParameter('TimeWin',           10, @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.addParameter('TimeStep',          [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.addParameter('TimeOnly',          false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
o = p.Results;

% Validate TimeOnly
if o.TimeOnly && isempty(o.TimeWin)
    warning('compute_mMSE: TimeOnly=true but TimeWin is empty — whole-epoch MSE will be computed.');
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
    % if isempty(o.TOI)
        step = iff(isempty(o.TimeStep), o.TimeWin/2, o.TimeStep);
        t0   = o.TimeWin/2;
        t1   = Tsec - o.TimeWin/2;
        if t1 < t0
            warning('compute_mMSE: TimeWin > data duration; disabling time-resolved output.');
            doTime = false; o.TimeOnly = false;
        else
            tCenters = t0:step:t1;
        end
    % else
    %     tCenters = o.TOI(:).';
    %     half     = o.TimeWin/2;
    %     tCenters = tCenters(tCenters>=half & tCenters<=Tsec-half);
    %     if isempty(tCenters), doTime = false; o.TimeOnly = false; end
    % end
end
nTOI = numel(tCenters);

% Sanity check on nTOI
if doTime && nTOI > 5000
    warning(['compute_mMSE: %d time windows requested (TimeWin=%.3gs, TimeStep=%.3gs). ' ...
        'This may be very slow or crash. Consider increasing TimeWin/TimeStep ' ...
        'or setting TimeOnly=true.'], nTOI, o.TimeWin, o.TimeStep);
end

% -------- Cap number of scales
S    = max(1, floor(o.num_scales));
maxS = floor(nSamp / o.MinSamplesPerBin);
if S > maxS
    warning('compute_mMSE: Reducing num_scales from %d to %d (need >=%d coarse bins).', ...
        S, maxS, o.MinSamplesPerBin);
    S = maxS;
end
if S < 1, S = 1; o.filter_mode = 'none'; end

% -------- Build frequency bands
%   filter_mode='none'       -> NaN (statistical path for all scales)
%   filter_mode='narrowband' -> annular bands [Nyq/(s+1), Nyq/s] +/- NBWiden
%                               Scale 1: HP-only [lo_eff, Inf]
nyq    = o.Fs / 2;
bands  = cell(1, S);
scales = cell(1, S);
if strcmpi(o.filter_mode,'none')
    for s = 1:S
        bands{s}  = [NaN NaN];
        scales{s} = num2str(s);
    end
else
    for s = 1:S
        lo_raw = nyq / (s+1);
        hi_raw = nyq / s;
        if s == 1
            lo_eff   = max(0, lo_raw * (1 - o.NBWiden));
            bands{s} = [lo_eff, Inf];
            scales{s}= sprintf('HP>%.3f Hz', lo_eff);
        else
            lo_eff = max(0,            lo_raw * (1 - o.NBWiden));
            hi_eff = min(nyq*(1-1e-6), hi_raw * (1 + o.NBWiden));
            if hi_eff <= lo_eff
                bands{s}  = [NaN NaN];
                scales{s} = 'invalid';
            else
                bands{s}  = [lo_eff, hi_eff];
                scales{s} = sprintf('[%.3f %.3f] Hz', lo_eff, hi_eff);
            end
        end
        if o.Progress && all(isfinite(bands{s}))
            fprintf('  Scale %2d band: %s\n', s, scales{s});
        end
    end
end

% -------- Pre-allocate
mse      = nan(nChan, S);
cg_len   = zeros(1, S);
mse_time = [];
if doTime, mse_time = nan(nChan, S, nTOI); end

% -------- Progress setup
useParScales = o.Parallel && ~isempty(ver('parallel'));
hWB = []; dq = [];
if o.Progress
    toStr  = tern(o.TimeOnly, ' | TimeOnly=ON', '');
    twStr  = tern(doTime, sprintf('%.3gs', o.TimeWin), 'off');
    parStr = tern(useParScales, 'ON', 'OFF');
    fprintf('mMSE: %d scales | Filter=%s | Coarse=%s | Parallel=%s | TimeWin=%s | nTOI=%d%s\n', ...
        S, o.filter_mode, ct, parStr, twStr, nTOI, toStr);
    if useParScales
        dq = parallel.pool.DataQueue;
        afterEach(dq, @(s) fprintf('  scale %2d done\n', s));
    else
        try
            hWB = waitbar(0, 'Computing mMSE...', 'Name', 'compute_mMSE');
        catch, hWB = []; end
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

% =========================================================================
%  SCALE LOOP — PARALLEL
% =========================================================================
if useParScales

    parfor s = 1:S
        feff     = o.Fs / s;
        Y        = dat;
        band     = bands{s};
        cg_len_s = 0;
        mse_s    = nan(nChan,1);
        mseslice = [];

        % --- Statistical path
        if strcmpi(o.filter_mode,'none') || any(isnan(band))

            % Scale 1 + std coarsing is trivially zero — skip
            if s == 1 && strcmpi(ct,'std')
                if o.Progress
                    fprintf('  [skip] scale 1: std of single-element window is zero.\n');
                end
            else
                [CG, nBins] = coarsegrain_stat(Y, s, ct);
                cg_len_s    = nBins;
                if nBins >= max(o.MinSamplesPerBin, o.m+1)
                    if ~o.TimeOnly
                        for ch = 1:nChan
                            mse_s(ch) = compute_SampEn(CG(ch,:), 'm', o.m, 'tau', o.tau, ...
                                'r', o.r, 'Parallel', false, 'Progress', false);
                        end
                    end
                    if doTime
                        mseslice = time_windows_one_scale(CG, feff, o, nTOI, tCenters);
                    end
                end
            end

        else
            % --- Narrowband filtskip path (Kosciessa 2020)
            isHPonly = isinf(band(2));

            % Filter
            Y = fir_zero_phase(Y, band, o.Fs, o.PadLen);

            % Decimation setup — scale 1 (HP-only) keeps full length
            if isHPonly
                nStarts   = 1;
                nBins_eff = size(Y,2);
            else
                nStarts   = s;
                nBins_eff = floor(size(Y,2)/s);
            end
            cg_len_s = nBins_eff;

            % FIR feasibility drop (bandpass only)
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

            % Whole-epoch SampEn — skipped when TimeOnly=true
            if ~o.TimeOnly
                vals = nan(nChan, nStarts);
                for is = 1:nStarts
                    CG = iff(nStarts==1, Y, Y(:, is:s:is+(nBins_eff-1)*s));
                    tmp = nan(nChan,1);
                    for ch = 1:nChan
                        tmp(ch) = compute_SampEn(CG(ch,:), 'm', o.m, 'tau', o.tau, ...
                            'r', o.r, 'Parallel', false, 'Progress', false);
                    end
                    vals(:,is) = tmp;
                end
                mse_s = mean(vals, 2, 'omitnan');
            end

            % Time-resolved — always runs when doTime
            if doTime
                feff_eff = tern(nStarts==1, o.Fs, feff);
                mseslice = time_windows_filtskip_one_scale(Y, s, nStarts, nBins_eff, feff_eff, o, nTOI, tCenters);
            end
        end

        % Write outputs
        mse(:,s)  = mse_s;    %#ok<PFBNS>
        cg_len(s) = cg_len_s; %#ok<PFBNS>
        if doTime && ~isempty(mseslice)
            mse_time(:,s,:) = mseslice; %#ok<PFBNS>
        end
        if ~isempty(dq) && o.Progress, send(dq,s); end
    end

    % =========================================================================
    %  SCALE LOOP — SERIAL
    % =========================================================================
else
    for s = 1:S
        if o.Progress
            fprintf('  scale %2d/%2d\n', s, S);
            if ~isempty(hWB) && isvalid(hWB)
                try waitbar(s/S, hWB, sprintf('Computing mMSE... (%d/%d)', s, S)); catch, end
            end
        end

        feff = o.Fs / s;
        Y    = dat;
        band = bands{s};

        % --- Statistical path
        if strcmpi(o.filter_mode,'none') || any(isnan(band))

            % Scale 1 + std coarsing is trivially zero — skip
            if s == 1 && strcmpi(ct,'std')
                if o.Progress
                    fprintf('  [skip] scale 1: std of single-element window is zero (use narrowband or coarsing=mean).\n');
                end
                continue
            end

            [CG, nBins] = coarsegrain_stat(Y, s, ct);
            cg_len(s)   = nBins;
            if nBins >= max([o.MinSamplesPerBin, o.m+1, o.MinWinSamples])
                if ~o.TimeOnly
                    for ch = 1:nChan
                        mse(ch,s) = compute_SampEn(CG(ch,:), 'm', o.m, 'tau', o.tau, ...
                            'r', o.r, 'Parallel', false, 'Progress', false);
                    end
                end
                if doTime
                    mse_time(:,s,:) = time_windows_one_scale(CG, feff, o, nTOI, tCenters);
                end
            end
            continue
        end

        % --- Narrowband filtskip path (Kosciessa 2020)
        isHPonly = isinf(band(2));

        % Filter
        Y = fir_zero_phase(Y, band, o.Fs, o.PadLen);

        % Decimation setup — scale 1 (HP-only) keeps full length
        if isHPonly
            nStarts   = 1;
            nBins_eff = size(Y,2);
        else
            nStarts   = s;
            nBins_eff = floor(size(Y,2)/s);
        end
        cg_len(s) = nBins_eff;

        % FIR feasibility drop (bandpass only)
        if ~isHPonly
            pbw_guard = band(2) - band(1);
            tw_guess  = max(0.15*max(pbw_guard,eps), o.Fs/size(Y,2));
            if pbw_guard < 2*tw_guess
                if o.Progress
                    fprintf('  [drop] scale %d: FIR band too narrow (pbw=%.4g Hz < 2*tw=%.4g Hz)\n', ...
                        s, pbw_guard, tw_guess);
                end
                continue
            end
        end
        if nBins_eff < max(4, o.m+1)
            if o.Progress
                fprintf('  [drop] scale %d: nBins=%d (<%d)\n', s, nBins_eff, max(4,o.m+1));
            end
            continue
        end

        % Whole-epoch SampEn — skipped when TimeOnly=true
        if ~o.TimeOnly
            vals = nan(nChan, nStarts);
            for is = 1:nStarts
                CG = iff(nStarts==1, Y, Y(:, is:s:is+(nBins_eff-1)*s));
                tmp = nan(nChan,1);
                for ch = 1:nChan
                    tmp(ch) = compute_SampEn(CG(ch,:), 'm', o.m, 'tau', o.tau, ...
                        'r', o.r, 'Parallel', false, 'Progress', false);
                end
                vals(:,is) = tmp;
            end
            mse(:,s) = mean(vals, 2, 'omitnan');
        end

        % Time-resolved — always runs when doTime
        if doTime
            feff_eff = tern(nStarts==1, o.Fs, feff);
            mse_time(:,s,:) = time_windows_filtskip_one_scale(Y, s, nStarts, nBins_eff, feff_eff, o, nTOI, tCenters);
        end
    end
end

% -------- Close waitbar
if ~isempty(hWB) && isvalid(hWB), try close(hWB); catch, end, end

% -------- Pack info and auto-plot
if doTime
    info.mse_time = mse_time;
    info.time_sec = tCenters;
    if o.PlotTime
        scIdx = unique(round(linspace(1, S, min(6, S))));
        plot_mMSE_timecourse(mse_time, tCenters, scales, 'ScaleIdx', scIdx);
    end
end
info.TimeOnly = o.TimeOnly;
info.cg_len   = cg_len;

end % compute_mMSE


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
ct = 'std';
if isempty(token), return; end
if isnumeric(token)
    switch round(token)
        case 0,    ct = 'mean';
        case 1,    ct = 'std';
        case 2,    ct = 'variance';
        otherwise, ct = 'std';
    end
else
    t = lower(regexprep(strtrim(char(token)),'[^a-z]',''));
    if     any(strcmp(t,{'std','sigma','standarddeviation'})), ct = 'std';
    elseif any(strcmp(t,{'variance','var','sigma2'})),         ct = 'variance';
    elseif any(strcmp(t,{'mean','avg','average','mu'})),       ct = 'mean';
    elseif any(strcmp(t,{'median','med'})),                    ct = 'median';
    end
end
end

function [CG, nBins] = coarsegrain_stat(Y, s, ct)
nChan = size(Y,1);
nFull = floor(size(Y,2)/s)*s;
if nFull == 0, CG = nan(nChan,0); nBins = 0; return; end
nBins = nFull / s;
Z = reshape(Y(:,1:nFull), nChan, s, nBins);
switch ct
    case 'std',      CG = squeeze(std(Z,  0, 2, 'omitnan'));
    case 'variance', CG = squeeze(var(Z,  0, 2, 'omitnan'));
    case 'mean',     CG = squeeze(mean(Z,    2, 'omitnan'));
    case 'median',   CG = squeeze(median(Z,  2, 'omitnan'));
    otherwise,       CG = squeeze(std(Z,  0, 2, 'omitnan'));
end
if isvector(CG), CG = CG(:)'; end
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
            'r', o.r, 'Parallel', false, 'Progress', false);
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
                'r', o.r, 'Parallel', false, 'Progress', false);
        end
    end
    if good, M(:,tiAbs) = acc ./ nStarts; end
end
end


%% =========================================================================
% FILTERS (FIR / IIR)
% =========================================================================

function Y = fir_zero_phase(Y, band, fs, padLen)
% hardcoded sensible defaults
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


function Y = ft_iir(Y, band, fs, ord, padLen)
nyq = fs/2; lo = band(1); hi = band(2);
if isinf(hi)
    [b,a] = butter(max(3,ceil(ord/2)), max(eps,lo/nyq), 'high');
    Y = iir_filtpad(Y, b, a, padLen);
else
    [bl,al] = butter(ord, min(0.999999,hi/nyq), 'low');
    [bh,ah] = butter(ord, max(eps,lo/nyq),      'high');
    Y = iir_filtpad(Y, bl, al, padLen);
    Y = iir_filtpad(Y, bh, ah, padLen);
end
end

function Y = iir_filtpad(Y, b, a, padLen)
T        = size(Y,2);
effOrder = max(length(a), length(b));
pad = iff(isempty(padLen), min(max(100,9*effOrder),floor((T-1)/2)), min(padLen,floor((T-1)/2)));
if pad > 0
    Yp = [fliplr(Y(:,1:pad)), Y, fliplr(Y(:,end-pad+1:end))];
    Yp = filtfilt(b, a, Yp.').';
    Y  = Yp(:, pad+1:end-pad);
else
    Y = filtfilt(b, a, Y.').';
end
end