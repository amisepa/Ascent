function [EEG, com] = ascent_compute(EEG, varargin)
% ascent_compute - EEGLAB plugin to compute entropy and complexity measures
%   on M/EEG or any biosignal time series (ECG, PPG, RR intervals, etc.).
%
% USAGE:
%   EEG = ascent_compute(EEG);                        % GUI mode
%   EEG = ascent_compute(EEG, 'measure', 'RCMFE');    % name-value mode
%
% REQUIRED INPUT:
%   EEG         - EEGLAB EEG structure
%
% OPTIONAL NAME-VALUE INPUTS (all have defaults):
%
%   General:
%     'measure'          - Complexity measure to compute. One of:
%                          'SampEn', 'FuzzEn', 'ExSEnt', 'FracDim',
%                          'HigFracDim', 'Aperiodic', 'MSE', 'mMSE',
%                          'MFE', 'CMFE', 'RCMFE', 'mvFuzzEn', 'RCmvMFE'
%     'domain'           - Data to run the measure on (default: 'channel'):
%                            'channel' - scalp channel signals (EEG.data)
%                            'ica'     - independent components (EEG.icaact).
%                          ICA mode requires an existing decomposition; IC
%                          results are back-projected to the scalp for topos
%                          using EEG.icawinv. 'chanlist' is ignored (all ICs).
%     'chanlist'         - Channel labels, e.g. {'Fz','Cz', 'POz'} or 'all' (default: all)
%                          Channel domain only.
%     'tau'              - Time lag for embedding (default: 1)
%     'm'                - Embedding dimension (default: 2)
%     'r'                - Similarity bound as fraction of SD (default: 0.15)
%     'vis'              - Plot outputs, logical (default: true)
%     'parallel'         - Use parallel computing, logical (default: true)
%     'progress'         - Track progress in console, logical (default: true)
%
%   Multiscale measures (MSE, mMSE, MFE, RCMFE, RCmvMFE):
%     'coarsing'         - Coarse-graining method: 'mean', 'median', 'std', 'var' (default: 'mean')
%     'num_scales'       - Number of scale factors (default: 30)
%     'zNorm'            - z-normalize per channel
%     'filter_mode'      - Scale filtering: 'none' or 'narrowband' (mMSE only)
%     'TimeOnly'         - compute time-resolved mMSE only
%     'TimeWin'          - Window length for time-resolved mMSE (default: [])
%     'TimeStep'         - Step between time windows for time-resolved mMSE
%                           (default: []). Setting this value to TimeWin/2
%                           gives 50% overlap (default).
%
%   Fuzzy measures (FuzzEn, MFE, RCMFE, RCmvMFE):
%     'n'                - Fuzzy power (default: 2)
%     'kernel'           - Fuzzy kernel type: 'exponential' (default)
%     'blocksize'        - Block size for fuzzy computation (default: 256)
%
%   Aperiodic (passed automatically from GUI, or set individually):
%     'aperiodicmode'    - Aperiodic fit mode: 'fixed' (default) or 'knee'
%     'maxpeaks'         - Max number of spectral peaks to fit (default: 6)
%     'minpeakheight'    - Min peak height above aperiodic (default: 0.05)
%     'peakthreshold'    - Peak detection threshold in SDs (default: 2.0)
%     'peakwidthlimits'  - [min max] peak width in Hz (default: [1 12])
%     'correctaperiodic' - Subtract aperiodic model from PSD, logical (default: true)
%                          Important: only valid with aperiodicmode = fixed.
%     'alphaband'        - [fMin fMax] Hz for the saved alpha power (default: [8 13])
%     'timeresolved'     - Compute sliding-window aperiodic fit, logical (default: false)
%                          slidWinSec (2 s) and slidOverlap (.5 for 50%).
%                          expose via GUI or name-value pairs in a future version.
%
% OUTPUTS:
%   Below, "signals" are scalp channels, or ICs when domain = 'ica'.
%
%   EEG.ascent.(measure).data               - computed measure [signals x scales]
%   EEG.ascent.(measure).domain             - 'channel' or 'ica'
%   EEG.ascent.(measure).electrode_labels   - signal labels ('IC1'.. in ICA mode)
%   EEG.ascent.(measure).electrode_locations- scalp channels (the ones ICs project onto)
%   EEG.ascent.(measure).icawinv            - mixing matrix        (ICA mode only)
%   EEG.ascent.(measure).icachansind        - channels ICA ran on  (ICA mode only)
%   EEG.ascent.(measure).scales             - scale vector (multiscale only)
%   EEG.ascent.aperiodic.data.exponent      - aperiodic exponent [signals x 1]
%   EEG.ascent.aperiodic.data.offset        - aperiodic offset [signals x 1]
%   EEG.ascent.aperiodic.data.psd           - raw PSD [signals x freqs]
%   EEG.ascent.aperiodic.data.freqs         - frequency vector
%   EEG.ascent.aperiodic.data.psd_corrected - aperiodic-corrected PSD (if requested).
%                                             RATIO to the 1/f fit (dimensionless).
%   EEG.ascent.aperiodic.data.alpha_raw     - total alpha power [signals x 1], µV²/Hz
%                                             (includes the aperiodic background)
%   EEG.ascent.aperiodic.data.alpha_osc     - alpha power above the 1/f fit
%                                             [signals x 1], µV²/Hz. NaN in knee mode.
%   EEG.ascent.aperiodic.data.exp_slid      - time-resolved exponent [channels x times]
%   EEG.ascent.aperiodic.data.off_slid      - time-resolved offset   [channels x times]
%   EEG.ascent.aperiodic.data.times_slid    - window centre times (s) [1 x times]
%   EEG.ascent.aperiodic.data.freqs_slid    - frequency vector for sliding PSD
%   EEG.ascent.aperiodic.data.psd_slid      - sliding PSD [channels x freqs x times]
%   EEG.ascent.aperiodic.data.psd_corr_slid - sliding corrected PSD (if requested)
%
% REFERENCES:
%   Sample Entropy:
%     Richman & Moorman (2000). Am J Physiol Heart Circ Physiol, 278(6).
%   Fuzzy Entropy:
%     Chen et al. (2009). IEEE Trans Biomed Eng, 56(5).
%   Multiscale Entropy (MSE):
%     Costa et al. (2002). Phys Rev Lett, 89(6).
%   Modified MSE (mMSE):
%     Kosciessa et al. (2020). NeuroImage, 211.
%   Multiscale Fuzzy Entropy (MFE / RCMFE):
%     Azami & Escudero (2016). Entropy, 18(10).
%   Aperiodic / 1/f parameterization (specparam):
%     Donoghue et al. (2020). Nature Neuroscience, 23.
%   Higuchi Fractal Dimension:
%     Higuchi (1988). Physica D, 31(2).
%   Fractal Dimension (box-counting):
%     Esteller et al. (2001). IEEE Trans Circuits Syst, 48(2).
%
% Copyright - Ascent EEGLAB plugin, Cedric Cannard, 2022

com    = '';
tstart = tic;

% Add path to subfolders
plugin_path = fileparts(which('eegplugin_ascent.m'));
addpath(fullfile(plugin_path, 'functions'));

% ----------------------------
% Input checks
% ----------------------------
if nargin < 1, help ascent_compute; return; end
if isempty(EEG.data),              error('Empty dataset.');          end
if isempty(EEG.chanlocs(1).labels),error('No channel labels.');      end
if ~isfield(EEG.chanlocs,'X') || isempty(EEG.chanlocs(1).X)
    error("Electrode locations are required. Go to 'Edit > Channel locations' and import coordinates.");
end
if isempty(EEG.ref)
    warning('EEG data not referenced. Referencing is strongly recommended (e.g., average-reference).');
end

% ----------------------------
% Get and validate all parameters
% ----------------------------
p = ascent_get_params(EEG, varargin{:});
if isempty(p)
    disp('User abort.');
    return
end

% Unpack struct into local variables for readability below
measure          = p.measure;
chanlist         = p.chanlist;
tau              = p.tau;
m                = p.m;
r                = p.r;
vis              = p.vis;
parallel         = p.parallel;
progress         = p.progress;
coarsing         = p.coarsing;
num_scales       = p.num_scales;
zNorm            = p.zNorm;
filter_mode      = p.filter_mode;
TimeWin          = p.TimeWin;
TimeStep         = p.TimeStep;
TimeOnly         = p.TimeOnly;
n                = p.n;
kernel_meth      = p.kernel_meth;
blocksize        = p.blocksize;
freqRange        = p.freqRange;
winSec           = p.winSec;
psdOverlap       = p.psdOverlap;
windowType       = p.windowType;
aperiodicMode    = p.aperiodicMode;
fitFreqRange     = p.fitFreqRange;
maxPeaks         = p.maxPeaks;
minPeakHeight    = p.minPeakHeight;
peakThreshold    = p.peakThreshold;
peakWidthLimits  = p.peakWidthLimits;
correctAperiodic = p.correctAperiodic;
alphaBand        = p.alphaBand;
timeResolved     = p.timeResolved;
domain           = p.domain;
slidWinSec       = p.slidWinSec;
slidOverlap      = p.slidOverlap;
slidNAvg         = p.slidNAvg;
slidPeakWidthLimits = p.slidPeakWidthLimits;

% ----------------------------
% Data domain validation (ICA)
% ----------------------------
% Done before any heavy lifting so a missing/inconsistent decomposition fails
% immediately with an actionable message rather than deep inside a measure.
isICA       = strcmpi(domain, 'ica');
icawinv     = [];
icachansind = [];
if isICA
    if isempty(EEG.icaweights) || isempty(EEG.icasphere)
        error(['ascent_compute: domain = ''ica'' requires an ICA decomposition, but ' ...
               'EEG.icaweights/EEG.icasphere are empty. Run ICA first ' ...
               '(e.g. EEG = pop_runica(EEG, ''icatype'', ''picard'')).']);
    end
    if isempty(EEG.icawinv)
        error('ascent_compute: EEG.icawinv is empty; cannot back-project ICs to the scalp.');
    end

    % Channels the decomposition was actually run on. EEGLAB allows ICA on a
    % subset, in which case icawinv rows map to these channels only -- never to
    % the full EEG.chanlocs.
    icachansind = EEG.icachansind;
    if isempty(icachansind), icachansind = 1:EEG.nbchan; end

    if size(EEG.icawinv,1) ~= numel(icachansind)
        error(['ascent_compute: EEG.icawinv has %d rows but ICA was run on %d channels. ' ...
               'The decomposition looks inconsistent with EEG.icachansind.'], ...
              size(EEG.icawinv,1), numel(icachansind));
    end

    % icaact is often not saved to disk to keep .set files small
    if isempty(EEG.icaact)
        disp('EEG.icaact is empty: reconstructing IC activations.');
        try
            EEG.icaact = eeg_getica(EEG);
        catch
            EEG.icaact = EEG.icaweights * EEG.icasphere * EEG.data(icachansind,:);
        end
    end
    icawinv = EEG.icawinv;
end

% ----------------------------
% Parpool management
% ----------------------------
p = gcp('nocreate');
if parallel
    if isempty(p)
        disp('Launching parallel pool...');
        parpool;
    end
else
    if ~isempty(p)
        disp('Shutting down parallel pool (parallel computing OFF).');
        delete(p);
    end
end

% ----------------------------
% Signal selection
% ----------------------------
% Extract data up front — avoids struct broadcast overhead in parfor.
% Note: from here on nChan counts *signals*, i.e. ICs when domain = 'ica'.
% chanlocs always describes scalp channels: in ICA mode it is the set the
% decomposition was run on, and is what icawinv back-projects onto.
fs = EEG.srate;
if isICA
    if numel(chanlist) < EEG.nbchan
        warning(['ascent_compute: channel selection is ignored when domain = ''ica''; ' ...
                 'using all %d ICs.'], size(EEG.icaact,1));
    end
    data       = EEG.icaact(:,:);          % flatten epochs, as in channel mode
    nChan      = size(data, 1);
    chanLabels = arrayfun(@(k) sprintf('IC%d',k), 1:nChan, 'UniformOutput', false);
    chanlocs   = EEG.chanlocs(icachansind);
else
    nChan = length(chanlist);
    if nChan > 1 && nChan < EEG.nbchan
        [~, chanIdx] = intersect({EEG.chanlocs.labels}, split(chanlist));
    else
        chanIdx = 1:EEG.nbchan;
    end
    data       = EEG.data(chanIdx, :);
    nChan      = size(data, 1);
    chanLabels = {EEG.chanlocs(chanIdx).labels};
    chanlocs   = EEG.chanlocs(chanIdx);
end

% Preallocate
if contains(lower(measure), {'mse','mmse','mfe','cmfe','rcmfe','rcmvmfe'})
    entropy = nan(nChan, num_scales);
else
    entropy = nan(nChan, 1);
end
scales = {};

% ----------------------------
% Compute
% ----------------------------
switch lower(measure)

    case 'sampen'
        % Sample Entropy — Richman & Moorman (2000)
        entropy = compute_SampEn(data, 'm', m, 'r', r, ...
            'Parallel', parallel, 'Progress', progress);

    case 'fuzzen'
        % Fuzzy Entropy — Chen et al. (2009)
        entropy = compute_FuzzEn(data, 'm', m, 'tau', tau, 'n', n, 'r', r, ...
            'Kernel', kernel_meth, ...
            'Parallel', parallel, 'Progress', progress);

    case 'exsent'
        % Extrema-Segmented Entropy
        [HD, HA, HDA, ~] = compute_ExSEnt(data, 'm', m, 'r', r, ...
            'lambda', 0.001, 'Plot', false, ...
            'Parallel', parallel, 'Progress', progress);

    case 'mvfuzzen'
        % Multivariate Fuzzy Entropy
        [entropy, ~, ~] = compute_mvFuzzEn(data, 'm', m, 'tau', tau, 'n', n, 'r', r, ...
            'Kernel', kernel_meth, 'BlockSize', blocksize, ...
            'Parallel', parallel, 'Progress', progress);

    case 'fracdim'
        % Box-counting Fractal Dimension — Esteller et al. (2001)
        [entropy, ~, ~] = compute_FracDim(data, ...
            'RejectBursts', true, 'WinFrac', 0.02, 'ZThresh', 6, ...
            'RobustFit', 'theilsen', 'ScaleTrimIQR', true, ...
            'Parallel', parallel, 'Progress', progress);

    case 'higfracdim'
        % Higuchi Fractal Dimension — Higuchi (1988)
        [entropy, SD, info] = compute_HigFracDim(data, ...
            'Kmax', 16, 'RobustFit', 'theilsen', 'MinScales', 4, ...
            'Parallel', true, 'Progress', true);

    case 'aperiodic'
        % Aperiodic (1/f) exponent and offset — Donoghue et al. (2020)

        % Step 1: Welch PSD (linear, µV²/Hz)
        [psd, freqs] = compute_psd(data, fs, ...
            'freqRange', freqRange,   ...
            'winSec',    winSec,      ...
            'overlap',   psdOverlap,  ...
            'window',    windowType,  ...
            'detrend',   true,        ...
            'Parallel',  parallel,    ...
            'Progress',  progress);

        % Step 2: Fit aperiodic model
        [exponent, offset, info] = compute_AperiodicFit(freqs, psd, ...
            'FreqRange',       fitFreqRange,    ...
            'AperiodicMode',   aperiodicMode,   ...
            'MaxPeaks',        maxPeaks,        ...
            'MinPeakHeight',   minPeakHeight,   ...
            'PeakThreshold',   peakThreshold,   ...
            'PeakWidthLimits', peakWidthLimits, ...
            'Parallel',        parallel,        ...
            'Progress',        progress);

        % Step 3: aperiodic model, band power, and corrected PSD (FIXED MODE ONLY)
        % Two different quantities, deliberately kept apart:
        %   psd_corrected — RATIO to the fit (dimensionless, ~1 where the fit is good)
        %   alpha_osc     — absolute oscillatory power in uV^2/Hz (alpha_raw minus fit)
        % They are not on a common scale; do not plot them against each other.
        [alpha_raw, alpha_osc, ap_model] = ...
            compute_AperiodicBandPower(freqs, psd, exponent, offset, alphaBand);

        psd_corrected = [];   % stays empty when correction is off or not applicable
        if strcmpi(aperiodicMode, 'fixed')
            if correctAperiodic
                psd_corrected = psd ./ ap_model;   % = 10.^(log10(psd) - (offset - exp*log10(f)))
                disp('Aperiodic component subtracted from PSD.');
            end
        else
            % The knee model is not a straight line in log-log space, so the
            % linear reconstruction the helper uses does not describe it.
            % alpha_raw is still valid (plain band-averaged PSD); the
            % oscillatory split is not.
            alpha_osc = nan(size(alpha_osc));
            if correctAperiodic
                warning('Aperiodic correction is only valid for fixed mode. Skipping correction.');
            end
        end

        % Step 4 (optional): sliding-window (time-resolved) aperiodic fit.
        if timeResolved
                [exp_slid, off_slid, times_slid, freqs_slid, psd_slid, psd_corr_slid, fiterr_slid, peaks_slid] = ...
                    compute_AperiodicFit_sliding(data, fs, ...
                    'freqRange',        freqRange,        ...
                    'aperiodicMode',    aperiodicMode,    ...
                    'maxPeaks',         maxPeaks,         ...
                    'minPeakHeight',    minPeakHeight,    ...
                    'peakThreshold',    peakThreshold,    ...
                    'peakWidthLimits',  slidPeakWidthLimits,  ...
                    'correctAperiodic', correctAperiodic, ...
                    'winSec',           slidWinSec,       ...
                    'winOverlap',       slidOverlap,      ...
                    'nAvg',             slidNAvg,         ...
                    'Parallel',         parallel,         ...
                    'Progress',         progress);
        end

    case 'mse'
        % % Multiscale Entropy — Costa et al. (2002)
        % entropy = nan(nChan, num_scales);
        % for iChan = 1:EEG.nbchan
        %     for iScale = 1:num_scales
        %         entropy(iChan,iScale) = compute_mse_costa(zscore(EEG.data(iChan,:)), m, r, tau, coarsing);
        %     end
        % end
        % scales = 1:num_scales;

        % Ascent improved MSE (BUG???)
        [entropy, scales] = compute_MSE(data, 'm', m, 'tau', tau, ...
            'coarsing', coarsing, 'num_scales', num_scales, ...
            'Parallel', parallel, 'Progress', progress);

    case 'mmse'
        % Modified MSE — Kosciessa et al. (2020)
        [entropy, scales, info] = compute_mMSE(data, 'Fs', fs, ...
            'm', m, 'tau', tau, 'r', r, 'num_scales', num_scales, 'coarsing', coarsing, ...
            'filter_mode', filter_mode, ...
            'TimeWin', TimeWin, 'TimeStep', TimeStep, 'TimeOnly', TimeOnly, ...
            'Parallel', parallel, 'Progress', progress);

    case 'mfe'
        % Multiscale Fuzzy Entropy — Azami & Escudero (2016)
        [entropy, scales] = compute_MFE(data, 'm', m, ...
            'tau', tau, 'r', r, 'coarsing', coarsing, 'num_scales', num_scales, ...
            'Parallel', parallel, 'Progress', progress);

    case 'cmfe'
        % Composite Multiscale Fuzzy Entropy (CMFE) — Azami & Escudero (2016)
        [entropy, scales] = compute_CMFE(data, 'm', m, 'tau', tau, ...
            'r', r, 'n', n, 'coarsing', coarsing, 'num_scales', num_scales, ...
            'Parallel', parallel, 'Progress', progress);

    case 'rcmfe'
        % Refined Composite Multiscale Fuzzy Entropy (RCMFE) — Azami & Escudero (2016)
        [entropy, scales] = compute_RCMFE(data, 'm', m, 'tau', tau, ...
            'r', r, 'n', n, 'coarsing', coarsing, 'num_scales', num_scales, ...
            'Parallel', parallel, 'Progress', progress);

    case 'rcmvmfe'
        % Refined Composite Multivariate Multiscale Fuzzy Entropy
        [entropy, scales] = compute_RCmvMFE(data, 'm', m, 'tau', tau, 'r', r, ...
            'coarsing', coarsing, 'num_scales', num_scales, ...
            'Parallel', parallel, 'Progress', progress);

    otherwise
        error('Unknown measure: %s. See help ascent_compute for valid options.', measure);
end


% ----------------------------
% Store outputs in EEG structure
% ----------------------------
switch lower(measure)
    case 'exsent'
        EEG.ascent.(measure).data.HD  = HD;
        EEG.ascent.(measure).data.HA  = HA;
        EEG.ascent.(measure).data.HDA = HDA;
    case 'aperiodic'
        EEG.ascent.(measure).data.exponent  = exponent;
        EEG.ascent.(measure).data.offset    = offset;
        EEG.ascent.(measure).data.psd       = psd;
        EEG.ascent.(measure).data.freqs     = freqs;
        % Band power saved so the figures plot what was computed, not their own
        % re-derivation. alpha_raw is total power (uV^2/Hz, includes the 1/f);
        % alpha_osc is the part above the fit (same unit). NaN in knee mode.
        EEG.ascent.(measure).data.alpha_raw = alpha_raw;
        EEG.ascent.(measure).data.alpha_osc = alpha_osc;
        EEG.ascent.(measure).params.alphaBand = alphaBand;
        EEG.ascent.(measure).params.freqRange       = freqRange;
        EEG.ascent.(measure).params.winSec          = winSec;
        EEG.ascent.(measure).params.overlap         = psdOverlap;
        EEG.ascent.(measure).params.window          = windowType;
        EEG.ascent.(measure).params.aperiodicMode   = aperiodicMode;
        EEG.ascent.(measure).params.fitFreqRange    = fitFreqRange;
        EEG.ascent.(measure).params.maxPeaks        = maxPeaks;
        EEG.ascent.(measure).params.minPeakHeight   = minPeakHeight;
        EEG.ascent.(measure).params.peakThreshold   = peakThreshold;
        EEG.ascent.(measure).params.peakWidthLimits = peakWidthLimits;
        EEG.ascent.(measure).params.timeResolved    = timeResolved;
        if ~isempty(psd_corrected)
            EEG.ascent.(measure).data.psd_corrected = psd_corrected;
        end
        if timeResolved
            EEG.ascent.(measure).data.exp_slid      = exp_slid;
            EEG.ascent.(measure).data.off_slid      = off_slid;
            EEG.ascent.(measure).data.times_slid    = times_slid;
            EEG.ascent.(measure).data.freqs_slid    = freqs_slid;
            EEG.ascent.(measure).data.psd_slid      = psd_slid;
            EEG.ascent.(measure).data.psd_corr_slid = psd_corr_slid;
            EEG.ascent.(measure).data.fiterr_slid = fiterr_slid;
            EEG.ascent.(measure).data.peaks_slid  = peaks_slid;
        end
    otherwise
        EEG.ascent.(measure).data = entropy;
end
EEG.ascent.(measure).domain              = domain;
EEG.ascent.(measure).electrode_labels    = chanLabels;   % IC1..ICn when domain = 'ica'
EEG.ascent.(measure).electrode_locations = chanlocs;     % scalp channels ICs project onto
if isICA
    % Saved so the back-projected figures can be reproduced from the .set alone
    EEG.ascent.(measure).icawinv     = icawinv;
    EEG.ascent.(measure).icachansind = icachansind;
end
if contains(lower(measure), {'mse','mmse','mfe','rcmfe','rcmvmfe'})
    EEG.ascent.(measure).scales = scales;
end
if exist('info','var')
    EEG.ascent.(measure).info = info;
end

% ----------------------------
% VISUALIZATIONS | PLOTS
% ----------------------------
if vis
    if nChan > 1
        % Remove scale 1 for std and var coarse-graining
        if ~isempty(coarsing) && contains(coarsing, {'std' 'sd' 'standard deviation' 'var' 'variance'})
            entropy(:,1) = [];
            scales(1) = [];
        end

        % Domain-specific arguments, appended to every ascent_plot call.
        % In channel mode this is empty and the calls are unchanged.
        icaArgs = {};
        if isICA
            icaArgs = {'ICA', true, 'icawinv', icawinv};
            if isfield(EEG,'dipfit') && ~isempty(EEG.dipfit)
                icaArgs = [icaArgs, {'dipfit', EEG.dipfit}];
            end
        end

        switch lower(measure)
            case 'exsent'
                ascent_plot(HD,  chanlocs, 'SampEn of durations',                    [], icaArgs{:});
                ascent_plot(HA,  chanlocs, 'SampEn of amplitudes',                   [], icaArgs{:});
                ascent_plot(HDA, chanlocs, 'Joint SampEn of durations & amplitudes', [], icaArgs{:});

            case 'aperiodic'
                ascent_plot(exponent, chanlocs, 'Aperiodic', [], offset, freqs, psd, psd_corrected, ...
                    icaArgs{:}, 'alphaPower', struct('raw',alpha_raw,'osc',alpha_osc,'band',alphaBand));
                if timeResolved
                    ascent_plot(exp_slid, chanlocs, 'Aperiodic', [], off_slid, ...
                        freqs_slid, psd_slid, psd_corr_slid, times_slid, icaArgs{:});
                end

            case 'mmse'
                hasTime = exist('info','var') && isstruct(info) && isfield(info,'mse_time') ...
                          && ~isempty(info.mse_time);
                if hasTime
                    plot_mMSE_timecourse(info.mse_time, info.time_sec, scales);
                end
                if ~TimeOnly && ~all(isnan(entropy(:)))
                    ascent_plot(entropy, chanlocs, measure, scales, icaArgs{:});
                end
            otherwise
                ascent_plot(entropy, chanlocs, measure, scales, icaArgs{:});
        end
    else
        disp("Ascent's visualizations require more than 1 signal.");
    end
end


% ----------------------------
% Command history (eegh)
% Build chanlist as a valid MATLAB cell-array literal so that replaying
% history from eegh reproduces the exact channel selection.
% ----------------------------
chanlist_str = ['{' strjoin(cellfun(@(c) sprintf('''%s''', c), cellstr(chanlist), ...
    'UniformOutput', false), ',') '}'];

if strcmpi(measure, 'aperiodic')
    com = sprintf( ...
        ['EEG = ascent_compute(EEG, ''measure'', ''%s'', ''domain'', ''%s'', '     ...
         '''chanlist'', %s, '                                                      ...
         '''vis'', %d, ''parallel'', %d, ''progress'', %d, '                       ...
         '''freqrange'', [%g %g], ''slidWinSec'', %g, ''slidOverlap'', %g, '               ...
         '''aperiodicmode'', ''%s'', ''fitfreqrange'', [%g %g], '                   ...
         '''maxpeaks'', %d, ''minpeakheight'', %g, ''peakthreshold'', %g, '        ...
         '''peakwidthlimits'', [%g %g], ''correctaperiodic'', %d, '                ...
         '''alphaband'', [%g %g], ''timeresolved'', %d;'],                          ...
        measure, domain, chanlist_str, vis, parallel, progress,                     ...
        freqRange(1), freqRange(2), slidWinSec, slidOverlap,                            ...
        aperiodicMode, fitFreqRange(1), fitFreqRange(2),                            ...
        maxPeaks, minPeakHeight, peakThreshold,                                     ...
        peakWidthLimits(1), peakWidthLimits(2),                                     ...
        correctAperiodic, alphaBand(1), alphaBand(2), timeResolved);
else
    com = sprintf( ...
        ['EEG = ascent_compute(EEG, ''measure'', ''%s'', ''domain'', ''%s'', ' ...
         '''chanlist'', %s, '                                                  ...
         '''tau'', %d, ''m'', %d, ''r'', %.3f, '                               ...
         '''vis'', %d, ''parallel'', %d, ''progress'', %d);'],                 ...
        measure, domain, chanlist_str, tau, m, r, vis, parallel, progress);
end

% Save dataset with new outputs
disp("Overwriting .set dataset to save the ASCENT outputs you computed. Note that this does not affect the rest of your EEG dataset.")
EEG = pop_saveset(EEG, 'filename', EEG.filename, 'filepath', EEG.filepath);

disp('Done!');
fprintf('Time to compute: %.2f minutes.\n', toc(tstart)/60);
disp('All outputs are stored in EEG.ascent.');