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
%     'zNorm'            - Varying-tolerance rescaling across scales, integer 0-4
%                          (default: 0 = OFF, i.e. fixed tolerance / classic Costa)
%                            1 = std, 2 = var, 3 = mad(mean), 4 = mad(median)
%                          of each coarse-grained series. MSE, MFE, CMFE and
%                          RCMFE only; ignored by mMSE and RCmvMFE.
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
%     'maxpeaks'         - Max number of spectral peaks to fit (default: 3)
%     'minpeakheight'    - Min peak height above aperiodic (default: 0.05)
%     'peakthreshold'    - Peak detection threshold in SDs (default: 2.0)
%     'peakwidthlimits'  - [min max] peak width in Hz (default: [1 12])
%     'correctaperiodic' - Subtract aperiodic model from PSD, logical (default: true)
%                          In knee mode the full Lorentzian model is removed.
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
%   EEG.ascent.(measure).ic_info            - per IC: ICLabel class, dipole RV, MNI
%                                             position, inside brain, AAL region and
%                                             lobe (ICA mode, with ICLabel or DIPFIT;
%                                             see ascent_ic_info)
%   EEG.ascent.(measure).scales             - scale vector (multiscale only)
%   EEG.ascent.Aperiodic.data.exponent      - aperiodic exponent [signals x 1]
%   EEG.ascent.Aperiodic.data.offset        - aperiodic offset [signals x 1]
%   EEG.ascent.Aperiodic.data.psd           - raw PSD [signals x freqs]
%   EEG.ascent.Aperiodic.data.freqs         - frequency vector
%   EEG.ascent.Aperiodic.data.psd_corrected - aperiodic-corrected PSD (if requested).
%                                             RATIO to the 1/f fit (dimensionless).
%   EEG.ascent.Aperiodic.data.alpha_raw     - total alpha power [signals x 1], µV²/Hz
%                                             (includes the aperiodic background)
%   EEG.ascent.Aperiodic.data.alpha_osc     - alpha power above the 1/f fit
%                                             [signals x 1], µV²/Hz (both modes;
%                                             knee mode uses the Lorentzian model)
%   EEG.ascent.Aperiodic.data.exp_slid      - time-resolved exponent [channels x times]
%   EEG.ascent.Aperiodic.data.off_slid      - time-resolved offset   [channels x times]
%   EEG.ascent.Aperiodic.data.times_slid    - window centre times (s) [1 x times]
%   EEG.ascent.Aperiodic.data.freqs_slid    - frequency vector for sliding PSD
%   EEG.ascent.Aperiodic.data.psd_slid      - sliding PSD [channels x freqs x times]
%   EEG.ascent.Aperiodic.data.psd_corr_slid - sliding corrected PSD (if requested)
%
% REFERENCES:
%   Sample entropy:
%       Richman, J. S., & Moorman, J. R. (2000). Physiological time-series analysis using
%       approximate entropy and sample entropy. Am J Physiol Heart Circ Physiol, 278,
%       H2039-H2049.
%   Fuzzy entropy:
%       Chen, W., Zhuang, J., Yu, W., & Wang, Z. (2009). Measuring complexity using
%       FuzzyEn, ApEn, and SampEn. Med Eng Phys, 31, 61-68.
%   ExSEnt:
%       Kamali, S., Baroni, F., & Varona, P. (2026). ExSEnt: Extrema-segmented entropy
%       analysis of time series. Nonlinear Dyn, 114, 980.
%   Fractal dimension (box-counting):
%       Esteller, R., Vachtsevanos, G., Echauz, J., & Litt, B. (2001). A comparison of
%       waveform fractal dimension algorithms. IEEE Trans Circuits Syst I, 48, 177-183.
%   Higuchi fractal dimension:
%       Higuchi, T. (1988). Approach to an irregular time series on the basis of the
%       fractal theory. Physica D, 31, 277-283.
%   Multiscale entropy (MSE):
%       Costa, M., Goldberger, A. L., & Peng, C.-K. (2002). Multiscale entropy analysis
%       of complex physiologic time series. Phys Rev Lett, 89, 068102.
%       Costa, M. D., & Goldberger, A. L. (2015). Generalized multiscale entropy
%       analysis: Application to quantifying the complex volatility of human heartbeat
%       time series. Entropy, 17, 1197-1203.
%   Modified MSE (mMSE):
%       Kosciessa, J. Q., Kloosterman, N. A., & Garrett, D. D. (2020). Standard
%       multiscale entropy reflects neural dynamics at mismatched temporal scales: What's
%       signal irregularity got to do with it? PLoS Comput Biol, 16, e1007885.
%       Kloosterman, N. A., Kosciessa, J. Q., Lindenberger, U., Fahrenfort, J. J., &
%       Garrett, D. D. (2020). Boosts in brain signal variability track liberal shifts in
%       decision bias. eLife, 9, e54201.
%   MFE, CMFE, RCMFE:
%       Wu, S.-D., Wu, C.-W., Lin, S.-G., Wang, C.-C., & Lee, K.-Y. (2013). Time series
%       analysis using composite multiscale entropy. Entropy, 15, 1069-1084.
%       Azami, H., Fernandez, A., & Escudero, J. (2017). Refined multiscale fuzzy entropy
%       based on standard deviation for biomedical signal analysis. Med Biol Eng Comput,
%       55, 2037-2052.
%       Azami, H., & Escudero, J. (2018). Coarse-graining approaches in univariate
%       multiscale sample and dispersion entropy. Entropy, 20, 138.
%   mvFuzzEn, RCmvMFE:
%       Ahmed, M. U., & Mandic, D. P. (2011). Multivariate multiscale entropy: A tool for
%       complexity analysis of multichannel data. Phys Rev E, 84, 061918.
%       Azami, H., Smith, K., & Escudero, J. (2016). MEMD-enhanced multivariate fuzzy
%       entropy for the evaluation of complexity in biomedical signals. Proc IEEE EMBC
%       2016, 3761-3764.
%       Azami, H., & Escudero, J. (2017). Refined composite multivariate generalized
%       multiscale fuzzy entropy: A tool for complexity analysis of multichannel signals.
%       Physica A, 465, 261-276.
%   Aperiodic (specparam) and time-resolved (SPRiNT):
%       Donoghue, T., Haller, M., Peterson, E. J., et al. (2020). Parameterizing neural
%       power spectra into periodic and aperiodic components. Nat Neurosci, 23,
%       1655-1665.
%       Wilson, L. E., da Silva Castanheira, J., & Baillet, S. (2022). Time-resolved
%       parameterization of aperiodic and periodic brain activity. eLife, 11, e77348.
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).

com    = '';
tstart = tic;

% Add path to subfolders
plugin_path = fileparts(which('eegplugin_ascent.m'));
addpath(fullfile(plugin_path, 'functions'));

% ----------------------------
% Input checks
% ----------------------------
if nargin < 1, help ascent_compute; return; end
% With no name-value pairs the GUI opens, and it can load a dataset (or the
% sample dataset) when none is loaded, so the dataset is checked after it
guiMode = isempty(varargin) || (numel(varargin) == 1 && isempty(varargin{1}));
if ~guiMode, check_dataset(EEG); end

% ----------------------------
% Get and validate all parameters
% ----------------------------
[p, EEGin] = ascent_get_params(EEG, varargin{:});
if isempty(p)
    disp('User abort.');
    return          % the dataset passed in, not one loaded in the window
end
EEG = EEGin;
if guiMode, check_dataset(EEG); end

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
IncludeScale1    = p.IncludeScale1;
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
pool = gcp('nocreate');   % not 'p': p.save is read at the end
if parallel
    if isempty(pool)
        disp('Launching parallel pool...');
        parpool;
    end
else
    if ~isempty(pool)
        disp('Shutting down parallel pool (parallel computing OFF).');
        delete(pool);
    end
end

% ----------------------------
% Signal selection
% ----------------------------
% Extract data up front - avoids struct broadcast overhead in parfor.
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
        % Sample Entropy - Richman & Moorman (2000)
        entropy = compute_SampEn(data, 'm', m, 'r', r, ...
            'Parallel', parallel, 'Progress', progress);

    case 'fuzzen'
        % Fuzzy Entropy - Chen et al. (2009)
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
        % Box-counting Fractal Dimension - Esteller et al. (2001)
        [entropy, ~, ~] = compute_FracDim(data, ...
            'RejectBursts', true, 'WinFrac', 0.02, 'ZThresh', 6, ...
            'RobustFit', 'theilsen', 'ScaleTrimIQR', true, ...
            'Parallel', parallel, 'Progress', progress);

    case 'higfracdim'
        % Higuchi Fractal Dimension - Higuchi (1988)
        [entropy, SD, info] = compute_HigFracDim(data, ...
            'Kmax', 16, 'RobustFit', 'theilsen', 'MinScales', 4, ...
            'Parallel', parallel, 'Progress', progress);

    case 'aperiodic'
        % Aperiodic (1/f) exponent and offset - Donoghue et al. (2020)

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

        % Step 3: aperiodic model, band power, and corrected PSD (all modes)
        % Two different quantities, deliberately kept apart:
        %   psd_corrected - RATIO to the fit (dimensionless, ~1 where the fit is good)
        %   alpha_osc     - absolute oscillatory power in uV^2/Hz (alpha_raw minus fit)
        % They are not on a common scale; do not plot them against each other.
        % Knee mode uses the full Lorentzian model (offset - log10(knee+f^exp)),
        % matching upstream specparam/FOOOF, whose flattened spectrum is
        % computed on every fit regardless of aperiodic mode.
        if strcmpi(aperiodicMode, 'knee')
            knee_vec = cellfun(@(k) k(2), info.knee);   % knee is param 2 in knee mode
            [alpha_raw, alpha_osc, ap_model] = ...
                compute_AperiodicBandPower(freqs, psd, exponent, offset, alphaBand, knee_vec);
        else
            [alpha_raw, alpha_osc, ap_model] = ...
                compute_AperiodicBandPower(freqs, psd, exponent, offset, alphaBand);
        end

        psd_corrected = [];   % stays empty when correction is off
        if correctAperiodic
            psd_corrected = psd ./ ap_model;   % ratio to the fitted model
            disp('Aperiodic component subtracted from PSD.');
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
        % % Multiscale Entropy - Costa et al. (2002)
        % entropy = nan(nChan, num_scales);
        % for iChan = 1:EEG.nbchan
        %     for iScale = 1:num_scales
        %         entropy(iChan,iScale) = compute_mse_costa(zscore(EEG.data(iChan,:)), m, r, tau, coarsing);
        %     end
        % end
        % scales = 1:num_scales;

        % Ascent improved MSE (BUG???)
        [entropy, scales] = compute_MSE(data, 'm', m, 'tau', tau, 'r', r, ...
            'coarsing', coarsing, 'num_scales', num_scales, 'zNorm', zNorm, ...
            'IncludeScale1', IncludeScale1, ...
            'Parallel', parallel, 'Progress', progress);

    case 'mmse'
        % Modified MSE - Kosciessa et al. (2020)
        [entropy, scales, info] = compute_mMSE(data, 'Fs', fs, ...
            'm', m, 'tau', tau, 'r', r, 'num_scales', num_scales, 'coarsing', coarsing, ...
            'filter_mode', filter_mode, ...
            'TimeWin', TimeWin, 'TimeStep', TimeStep, 'TimeOnly', TimeOnly, ...
            'Parallel', parallel, 'Progress', progress);

    case 'mfe'
        % Multiscale Fuzzy Entropy - Azami & Escudero (2016)
        [entropy, scales] = compute_MFE(data, 'm', m, ...
            'tau', tau, 'r', r, 'coarsing', coarsing, 'num_scales', num_scales, ...
            'zNorm', zNorm, 'Parallel', parallel, 'Progress', progress);

    case 'cmfe'
        % Composite Multiscale Fuzzy Entropy (CMFE) - Azami & Escudero (2016)
        [entropy, scales] = compute_CMFE(data, 'm', m, 'tau', tau, ...
            'r', r, 'n', n, 'coarsing', coarsing, 'num_scales', num_scales, ...
            'zNorm', zNorm, 'Parallel', parallel, 'Progress', progress);

    case 'rcmfe'
        % Refined Composite Multiscale Fuzzy Entropy (RCMFE) - Azami & Escudero (2016)
        [entropy, scales] = compute_RCMFE(data, 'm', m, 'tau', tau, ...
            'r', r, 'n', n, 'coarsing', coarsing, 'num_scales', num_scales, ...
            'zNorm', zNorm, 'IncludeScale1', IncludeScale1, ...
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
        % alpha_osc is the part above the fit (same unit, both modes).
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
    % ICLabel class and equivalent dipole of each IC, when the dataset has them
    icl = [];
    if isfield(EEG,'etc') && isstruct(EEG.etc) && isfield(EEG.etc,'ic_classification') ...
            && isfield(EEG.etc.ic_classification,'ICLabel')
        icl = EEG.etc.ic_classification.ICLabel;
    end
    dip = [];
    if isfield(EEG,'dipfit') && ~isempty(EEG.dipfit), dip = EEG.dipfit; end
    if ~isempty(icl) || ~isempty(dip)
        EEG.ascent.(measure).ic_info = ascent_ic_info(size(icawinv,2), icl, dip, chanlocs, [], icawinv);
    end
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
  try   % a plotting error must not lose the results computed above
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
            if ~isempty(icl), icaArgs = [icaArgs, {'iclabel', icl}]; end
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
  catch ME
    warning('ascent_compute: plotting failed (%s); the results are kept.', ME.message);
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
% ⚠ GUARD (2026-10-04): this save once OVERWROTE sourced data files. Batch
% scripts that compute on a sliced copy of a dataset (e.g. tmp = EEG;
% tmp.data = EEG.data(:, s1:s2)) inherit EEG.filename/filepath, so the save
% silently replaced the ORIGINAL .set on disk with the one-block copy.
% Behaviour now:
%   'save', 'auto' (default): save ONLY when the dataset has not been
%       resampled/referenced/filtered/sliced since load, i.e. when saving is
%       actually safe; otherwise SKIP with a notice (outputs are in EEG.ascent
%       in memory either way - callers who need them on disk should save under
%       their own filename).
%   'save', 'never': never touch disk.
%   'save', 'always': old behaviour (explicitly requested overwrite).
if isfield(p, 'save') && ~isempty(p.save)
    saveOpt = lower(char(p.save));
else
    saveOpt = 'auto';
end
doSave = false;
switch saveOpt
    case 'never'
        doSave = false;
    case 'always'
        doSave = true;
    case 'auto'
        % Save ONLY when the in-memory dataset still matches the file it
        % claims to come from (channel count, sample count, event count).
        % The known destructive pattern - batch scripts computing on a
        % SLICED copy that inherited filename/filepath - always changes
        % pnts and/or events, and is thereby refused. A plain GUI flow
        % (load -> compute -> save, no processing in between) still saves
        % exactly as before. Anything that cannot be verified is not saved.
        try
            assert(~isempty(EEG.filename) && ~isempty(EEG.filepath));
            % never the sample dataset shipped with the plugin
            assert(~(strcmpi(EEG.filename, 'ascent_sample_data.set') && ...
                     exist(fullfile(EEG.filepath, 'eegplugin_ascent.m'), 'file') == 2));
            fpath = fullfile(EEG.filepath, EEG.filename);
            assert(exist(fpath, 'file') == 2);
            % '-mat': a .set file is read as ASCII otherwise, which errored and
            % made this check refuse every save since 2026-10-04
            if ismember('EEG', who('-file', fpath))          % older .set: one EEG variable
                S0 = load('-mat', fpath, 'EEG');  S0 = S0.EEG;
            else
                S0 = load('-mat', fpath, 'nbchan', 'pnts', 'event');
            end
            sigNow  = sprintf('%d:%d:%d', EEG.nbchan, EEG.pnts, numel(EEG.event));
            sigFile = sprintf('%d:%d:%d', S0.nbchan, S0.pnts, numel(S0.event));
            doSave = strcmp(sigNow, sigFile);
        catch
            doSave = false;   % cannot verify -> do not risk overwriting
        end
    otherwise
        warning('ascent_compute: unknown ''save'' option ''%s''; not saving.', saveOpt);
        doSave = false;
end
if doSave
    disp("Overwriting .set dataset to save the ASCENT outputs you computed. Note that this does not affect the rest of your EEG dataset.")
    EEG = pop_saveset(EEG, 'filename', EEG.filename, 'filepath', EEG.filepath);
else
    disp("ASCENT outputs kept in EEG.ascent (in memory only; dataset on disk NOT overwritten).")
end

disp('Done!');
fprintf('Time to compute: %.2f minutes.\n', toc(tstart)/60);
disp('All outputs are stored in EEG.ascent.');


%% Dataset checks (before computing)
function check_dataset(EEG)   % no closing END: the main function has none
if isempty(EEG) || ~isfield(EEG,'data') || isempty(EEG.data), error('Empty dataset.'); end
if isempty(EEG.chanlocs) || isempty(EEG.chanlocs(1).labels), error('No channel labels.'); end
if ~isfield(EEG.chanlocs,'X') || isempty(EEG.chanlocs(1).X)
    error("Electrode locations are required. Go to 'Edit > Channel locations' and import coordinates.");
end
if isempty(EEG.ref)
    warning('EEG data not referenced. Referencing is strongly recommended (e.g., average-reference).');
end
