%% Tutorial for the ASCENT EEGLAB plugin.
% 
% Aperiodic Signal Complexity Estimation for Neurophysiological Time series
% 
% Please cite when using this plugin or algorithms:
%   Cannard, C., & Delorme, A. (2022). An open-source EEGLAB plugin for 
%   computing entropy-based measures on MEEG signals.
% 
% Copyright (C) - ASCENT EEGLAB PLUGIN - Cedric Cannard, 2021-2025

clear; close all; clc

% launch eeglab
eeglab; close;
% pop_editoptions('option_parallel', 1); % turn parrallel computing on (1) or off (0)

% if you haven'installed the plugin yet, either go to File > Manage EEGLAB
% extensions > search ascent_compute > Install
% or clone the github directory in the EEGLAB plugins folder
% or donwload the github repo and unzip it in the EEGLAB plugins folder

% Load provided sample EEG data from the tutorial directory 
% (2 minutes of resting state eyes-closed, 64-channel Biosemi).
% Dataset source: 
%   Cannard, C., Wahbeh, H., & Delorme, A. (2021).
%   Validating the wearable MUSE headset for EEG spectral analysis and Frontal 
%   Alpha Asymmetry. IEEE International Conference on Bioinformatics
%   and Biomedicine (BIBM) (pp. 3603-3610). 

pluginPath = fileparts(which('eegplugin_ascent.m'));
addpath(genpath(pluginPath))
cd(pluginPath)
EEG = pop_loadset('filename','ascent_sample_data.set','filepath',fullfile(pluginPath));

disp("Are you ready for your 1st ascent?!")
disp("Run each cell below one by one by clicking in the cell and then pressing CTRL + ENTER")

%% Launch main GUI

EEG = ascent_compute(EEG);  % or Tools > Compute entropy

%% SampEn via command line with default parameters

% Compute Sample entropy with command line using default parameters
EEG = ascent_compute(EEG, 'measure', 'SampEn');
% print(gcf, 'SampEn_topo.png','-dpng','-r300');   % 300 dpi .png

% % Outputs can be found in:
% EEG.ascent.SampEn

% Fine-tuning parameters (for demonstration only)
[EEG, com] = ascent_compute(EEG, 'measure', 'SampEn', ...
    'chanlist', 'Cz Fz F3 F4 Pz Oz O1 O2', ...   % channel selection
    'tau', 2, ...
    'm', 3, ...
    'r', .3, ...
    'parallel', false, ...
    'progress', true, ...
    'vis', true);

%% Fuzzy entropy (FuzzEn)

% default parameters
EEG = ascent_compute(EEG, 'measure', 'FuzzEn');

% Fine-tuning parameters (for demonstration only)
EEG = ascent_compute(EEG, 'measure', 'FuzzEn', ...
    'n', 2, ...             % fuzy power (default = 2)
    'kernel','gaussian');    % default: 'exponential'


%% Extrema-Segmented Entropy (ExSEnt)

EEG = ascent_compute(EEG, 'measure', 'ExSEnt', 'r', .15);


%% Fractal dimension

EEG = ascent_compute(EEG, 'measure', 'HigFracDim');  % 'HigFracDim' = Higushi version; 'FracDim' = box-counting version

%% Aperiodic exponent (i.e. slope) and offset

EEG = ascent_compute(EEG, 'measure', 'Aperiodic');

%% Multiscale entropy (MSE)

EEG = ascent_compute(EEG, 'measure', 'MSE', ...
    'coarsing', 'mean', ...     % 'median' 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 20, ...       % number of scale factors to compute (default = 20; range = 5-100 depending on sample rate)
    'parallel', true, 'progress', true);


%% Modified MSE (mMSE)

EEG = ascent_compute(EEG, 'measure', 'mMSE', ...
    'coarsing', 'mean', ... % 'median' 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 30, ...
    'filter_mode', 'none', ...  %  'narrowband' (default), 'none' (normal MSE)
    'parallel', true, 'progress', true);

% ascent_plot(EEG.ascent.mMSE.data, EEG.chanlocs,'mMSE',EEG.ascent.mMSE.scales)


%% Modified MSE (mMSE) - Time-course mode 

% on contninuous resting-state data
EEG = ascent_compute(EEG, 'measure', 'mMSE', ...
    'num_scales', 9, ...
    'TimeWin',    8, ...   % 10 s: stable through scale ~25 at 256 Hz
    'TimeStep',   4, ...   % 50% overlap
    'TimeOnly',  true, ...
    'Parallel',  true, ...
    'Progress',  true);


%% Multiscale Fuzzy Entropy (MFE)

EEG = ascent_compute(EEG, 'measure', 'MFE', ...
    'coarsing', 'mean', ...     % 'mean' (default), 'median', 'trimmed mean', 'std', 'var'
    'num_scales', 30, ...       % number of scale factors to compute (default = 20; range = 5-100 depending on sample rate)
    'n', 2, ...                 % fuzzy power (default = 2)
    'parallel', true, 'progress', true);

%% Composite Multiscale Fuzzy Entropy (CMFE)

EEG = ascent_compute(EEG, 'measure', 'CMFE', ...
    'coarsing', 'mean', ...     %  'mean' (default) 'median'  'trimmed mean' 'std' 'var'
    'num_scales', 50, ...       % number of scale factors to compute (default = 20; range = 5-100 depending on sample rate)
    'n', 2, ...                 % fuzzy power (default = 2)
    'parallel', true, 'progress', true);


%% Refined Composite Multiscale Fuzzy Entropy (RCMFE)

EEG = ascent_compute(EEG, 'measure', 'RCMFE', ...
    'coarsing', 'std', ...     % 'median' 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 50, ...       % number of scale factors to compute (default = 20; range = 2-100 depending on sample rate)
    'n', 2, ...                 % fuzzy power (default = 2)
    'parallel', true, 'progress', true);

ascent_plot(EEG.ascent.RCMFE.data, EEG.chanlocs,'RCMFE',EEG.ascent.RCMFE.scales)


%% Multivariate Fuzzy entropy (mvFuzzEn)
% WARNING: use only for low-channel count (e.g. <10 EEG channels) or on few
% PCA components, ROIs, or very short segments, otherwise, extremely long to
% compute. 

% EEG = ascent_compute(EEG, 'measure', 'mvFuzzEn');  % all channels

EEG = ascent_compute(EEG, 'measure', 'mvFuzzEn', 'chanlist', 'Fpz Fz Cz Pz Iz Oz POz'); 

%% Refined Composite Multivariate Generalized Multiscale Fuzzy Entropy (RCmvMFE)
% WARNING: scale 1 is much longer than the rest of the scales.

EEG = ascent_compute(EEG, 'measure', 'RCmvMFE', 'chanlist', 'Fpz Fz Cz Pz Iz Oz POz', ...
    'coarsing', 'mean', ...     % 'median' 'mean' 'trimmed mean' 'std' (default) 'var'
    'num_scales', 10, ...       % number of scale factors to compute (default = 20; range = 5-100 depending on sample rate)
    'n', 2, ...                % fuzzy power (default = 2)
    'parallel', true, 'progress', true);

%% Multiscale entropy on indepent components (ICs) instead of 
% scalp channel signals. note: ICA was precomputed on the sample dataset.

if isempty(EEG.icaact)
    EEG.icaact = EEG.icaweights * EEG.icasphere * EEG.data; % sometimes not saved on disk to save space
end

% [entropy, scales] = compute_MSE(EEG.icaact, 'coarsing', 'mean', 'num_scales', 30);
[entropy, scales] = compute_RCMFE(EEG.icaact, 'coarsing', 'sd', 'num_scales', 30);

ascent_plot(entropy, EEG.chanlocs, 'RCMFE', scales, ...
    'ICA', true, 'icawinv', EEG.icawinv);
% xlim([2:size(entropy,2)])

%% Aperiodic parametrization on indepent components (ICs) instead of 
% scalp channel signals. note: ICA was precomputed on the sample dataset.

% Compute Welch PSD (linear)
[psd, freqs] = compute_psd(EEG.icaact, EEG.srate);

% Fit aperiodic model
%     'aperiodicmode'    - Aperiodic fit mode: 'fixed' (default) or 'knee'
%     'maxpeaks'         - Max number of spectral peaks to fit (default: 6)
%     'minpeakheight'    - Min peak height above aperiodic (default: 0.05)
%     'peakthreshold'    - Peak detection threshold in SDs (default: 2.0)
%     'peakwidthlimits'  - [min max] peak width in Hz (default: [1 12])
%     'correctaperiodic' - Subtract aperiodic model from PSD, logical (default: true)
[exponent, offset, info] = compute_AperiodicFit(freqs, psd);

% subtract aperiodic model from PSD
log_freqs = log10(freqs);
psd_corrected = nan(size(psd));
for iChan = 1:size(psd,1)
    ap_model = offset(iChan) - exponent(iChan) .* log_freqs;
    psd_corrected(iChan,:) = 10.^(log10(psd(iChan,:)) - ap_model);
end

% fallback plot with bar graph instead of topos (and EEG channels instead
% of components)
ascent_plot(exponent, EEG.chanlocs, 'mfe', scales);

% Plot with back-projected IC topographies
ascent_plot(exponent, EEG.chanlocs, 'aperiodic', [], ...
    offset, freqs, psd, psd_corrected, ...
    'ICA', true, 'icawinv', EEG.icawinv);
