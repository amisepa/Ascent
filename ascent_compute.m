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
%     'chanlist'         - Channel labels, e.g. {'Fz','Cz', 'POz'} or 'all' (default: all)
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
%     'timeresolved'     - Compute sliding-window aperiodic fit, logical (default: false)
%                          slidWinSec (2 s) and slidOverlap (.5 for 50%).
%                          expose via GUI or name-value pairs in a future version.
%
% OUTPUTS:
%   EEG.ascent.(measure).data               - computed measure [channels x scales]
%   EEG.ascent.(measure).electrode_labels   - channel labels used
%   EEG.ascent.(measure).electrode_locations- channel locations used
%   EEG.ascent.(measure).scales             - scale vector (multiscale only)
%   EEG.ascent.aperiodic.data.exponent      - aperiodic exponent [channels x 1]
%   EEG.ascent.aperiodic.data.offset        - aperiodic offset [channels x 1]
%   EEG.ascent.aperiodic.data.psd           - raw PSD [channels x freqs]
%   EEG.ascent.aperiodic.data.freqs         - frequency vector
%   EEG.ascent.aperiodic.data.psd_corrected - aperiodic-corrected PSD (if requested)
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
addpath(genpath(plugin_path));

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
timeResolved     = p.timeResolved;
slidWinSec       = p.slidWinSec;
slidOverlap      = p.slidOverlap;
slidNAvg         = p.slidNAvg;
slidPeakWidthLimits = p.slidPeakWidthLimits;

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
% Channel indexing
% ----------------------------
nChan = length(chanlist);
if nChan > 1 && nChan < EEG.nbchan
    [~, chanIdx] = intersect({EEG.chanlocs.labels}, split(chanlist));
else
    chanIdx = 1:EEG.nbchan;
end

% Extract data — avoids struct broadcast overhead in parfor
data       = EEG.data(chanIdx, :);
fs         = EEG.srate;
nChan      = size(data, 1);
chanLabels = {EEG.chanlocs(chanIdx).labels};
chanlocs   = EEG.chanlocs(chanIdx);

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
        [psd, freqs] = compute_psd(EEG.data, EEG.srate, ...
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

        % Step 3: subtract aperiodic model from static PSD (ONLY VALID FOR FIXED MODE)
        % Model: L(f) = offset - exponent * log10(f)
        % Correction: log10(psd_corrected) = log10(psd) - L(f)
        if correctAperiodic && strcmpi(aperiodicMode, 'fixed')
            log_freqs = log10(freqs);
            psd_corrected = nan(size(psd));
            for iChan = 1:nChan
                ap_model = offset(iChan) - exponent(iChan) .* log_freqs;
                psd_corrected(iChan,:) = 10.^(log10(psd(iChan,:)) - ap_model);
            end
            disp('Aperiodic component subtracted from PSD.');
        else
            warning('Aperiodic correction is only valid for fixed mode. Skipping correction.');
        end

        % Step 4 (optional): sliding-window (time-resolved) aperiodic fit.
        if timeResolved
                [exp_slid, off_slid, times_slid, freqs_slid, psd_slid, psd_corr_slid, fiterr_slid, peaks_slid] = ...
                    compute_AperiodicFit_sliding(EEG.data, EEG.srate, ...
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
        % Multiscale Entropy — Costa et al. (2002)
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
        EEG.ascent.(measure).data.exponent = exponent;
        EEG.ascent.(measure).data.offset   = offset;
        EEG.ascent.(measure).data.psd      = psd;
        EEG.ascent.(measure).data.freqs    = freqs;
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
        if correctAperiodic
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
EEG.ascent.(measure).electrode_labels    = chanLabels;
EEG.ascent.(measure).electrode_locations = chanlocs;
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

        switch lower(measure)
            case 'exsent'
                ascent_plot(HD,  chanlocs, 'SampEn of durations',                    []);
                ascent_plot(HA,  chanlocs, 'SampEn of amplitudes',                   []);
                ascent_plot(HDA, chanlocs, 'Joint SampEn of durations & amplitudes', []);

            case 'aperiodic'
                if correctAperiodic
                    ascent_plot(exponent, chanlocs, 'Aperiodic', [], offset, freqs, psd, psd_corrected);
                else
                    ascent_plot(exponent, chanlocs, 'Aperiodic', [], offset, freqs, psd);
                end
                if timeResolved
                    ascent_plot(exp_slid, chanlocs, 'Aperiodic', [], off_slid, ...
                        freqs_slid, psd_slid, psd_corr_slid, times_slid);
                end

            case 'mmse'
                hasTime = exist('info','var') && isstruct(info) && isfield(info,'mse_time') ...
                          && ~isempty(info.mse_time);
                if hasTime
                    plot_mMSE_timecourse(info.mse_time, info.time_sec, scales);
                end
                if ~TimeOnly && ~all(isnan(entropy(:)))
                    ascent_plot(entropy, chanlocs, measure, scales);
                end
            otherwise
                ascent_plot(entropy, chanlocs, measure, scales);
        end
    else
        disp("Ascent's visualizations require more than 1 channel.");
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
        ['EEG = ascent_compute(EEG, ''measure'', ''%s'', ''chanlist'', %s, '       ...
         '''vis'', %d, ''parallel'', %d, ''progress'', %d, '                       ...
         '''freqrange'', [%g %g], ''slidWinSec'', %g, ''slidOverlap'', %g, '               ...
         '''aperiodicmode'', ''%s'', ''fitfreqrange'', [%g %g], '                   ...
         '''maxpeaks'', %d, ''minpeakheight'', %g, ''peakthreshold'', %g, '        ...
         '''peakwidthlimits'', [%g %g], ''correctaperiodic'', %d, '                ...
         '''timeresolved'', %d;'],                                                   ...
        measure, chanlist_str, vis, parallel, progress,                             ...
        freqRange(1), freqRange(2), slidWinSec, slidOverlap,                            ...
        aperiodicMode, fitFreqRange(1), fitFreqRange(2),                            ...
        maxPeaks, minPeakHeight, peakThreshold,                                     ...
        peakWidthLimits(1), peakWidthLimits(2),                                     ...
        correctAperiodic, timeResolved);
else
    com = sprintf( ...
        ['EEG = ascent_compute(EEG, ''measure'', ''%s'', ''chanlist'', %s, ' ...
         '''tau'', %d, ''m'', %d, ''r'', %.3f, '                             ...
         '''vis'', %d, ''parallel'', %d, ''progress'', %d);'],               ...
        measure, chanlist_str, tau, m, r, vis, parallel, progress);
end

% Save dataset with new outputs
disp("Overwriting .set dataset to save the ASCENT outputs you computed. Note that this does not affect the rest of your EEG dataset.")
EEG = pop_saveset(EEG, 'filename', EEG.filename, 'filepath', EEG.filepath);

disp('Done!');
fprintf('Time to compute: %.2f minutes.\n', toc(tstart)/60);
disp('All outputs are stored in EEG.ascent.');