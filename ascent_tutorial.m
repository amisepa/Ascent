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
% 
% or clone the github directory in the EEGLAB plugins folder
% or donwload the github repo and unzip it in the EEGLAB plugins folder

% Load provided sample EEG data from the tutorial directory 
% (2 minutes of resting state eyes-closed, 64-channel Biosemi).
% Dataset source: 
%   Cannard, C., Wahbeh, H., & Delorme, A. (2021, December).
%   Validating the wearable MUSE headset for EEG spectral analysis and Frontal 
%   Alpha Asymmetry. In 2021 IEEE International Conference on Bioinformatics
%   and Biomedicine (BIBM) (pp. 3603-3610). IEEE.
pluginPath = fileparts(which('eegplugin_ascent.m'));
addpath(genpath(pluginPath))
cd(pluginPath)
EEG = pop_loadset('filename','ascent_sample_data.set','filepath',fullfile(pluginPath));

disp("Are you ready for your 1st ascent?!")
disp("Run each cell below one by one by clicking in the cell and then pressing CTRL + ENTER")

%% Sample Entropy (SampEn) via GUI

EEG = ascent_compute(EEG);  % or Tools > Compute entrop y

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
    'r', .5, ...
    'parallel', false, ...
    'progress', true, ...
    'vis', true);

%% Fuzzy entropy (FuzzEn)

% default parameters
EEG = ascent_compute(EEG, 'measure', 'FuzzEn');

% Fine-tuning parameters (for demonstration only)
EEG = ascent_compute(EEG, 'measure', 'FuzzEn', ...
    'n', 2, ...
    'kernel','gaussian', ...    % default: exponential
    'blocksize', 128);          % default: 256


%% Extrema-Segmented Entropy (ExSEnt)

EEG = ascent_compute(EEG, 'measure', 'ExSEnt', 'r', .15);

%% Fractal dimension (votality)

EEG = ascent_compute(EEG, 'measure', 'FracDim');

%% Multiscale entropy (MSE)

EEG = ascent_compute(EEG, 'measure', 'MSE', ...
    'coarsing', 'median', ...     % 'median' (default) 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 50, ...       % number of scale factors to compute (default = 20; range = 5-100 depending on sample rate)
    'parallel', true, 'progress', true);


%% Modified MSE (mMSE)

EEG = ascent_compute(EEG, 'measure', 'mMSE', ...
    'coarsing', 'median', ... % 'median' (default) 'mean' 'std' 'variance'
    'num_scales', 50, ...
    'filter_mode', 'narrowband', ...  %  'narrowband' (default, annuli), 'none'
    'parallel', true, 'progress', true);

% ascent_plot(EEG.ascent.mMSE.data, EEG.chanlocs,'mMSE',EEG.ascent.mMSE.scales)


%% Modified MSE (mMSE) - Time-resolved version

EEG = ascent_compute(EEG, 'measure', 'mMSE', ...
    'coarsing', 'median', ...  % 'median' (default) 'mean' 'std' 'variance'
    'num_scales', 15, ...
    'filter_mode', 'narrowband', ...  %  'narrowband' (default, annuli), 'none'
    'TimeWin', 4, ...       %  window length (in s; default = 2 for continuous data)
    'TimeStep', 2, ...    % step between centers (s); default = TimeWin/2
    'parallel', true, 'progress', true);


%% Multiscale Fuzzy Entropy (MFE)

EEG = ascent_compute(EEG, 'measure', 'MFE', ...
    'coarsing', 'median', ...     % 'median' (default) 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 30, ...       % number of scale factors to compute (default = 20; range = 5-100 depending on sample rate)
    'n', 2, ...                 % fuzzy power (default = 2)
    'parallel', true, 'progress', true);

%% Refined Composite Multiscale Fuzzy Entropy (RCMFE)

EEG = ascent_compute(EEG, 'measure', 'RCMFE', ...
    'coarsing', 'std', ...     % 'median' (default) 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 10, ...       % number of scale factors to compute (default = 20; range = 5-100 depending on sample rate)
    'n', 2, ...                 % fuzzy power (default = 2)
    'parallel', true, 'progress', true);



%% Multivariate Fuzzy entropy (mvFuzzEn)
% WARNING: use only for low-channel count (e.g. <10 EEG channels) or on few
% PCA components, ROIs, or very short segments, otherwise, extremely long to
% compute. 

EEG = ascent_compute(EEG, 'measure', 'mvFuzzEn', 'chanlist', 'Fz Cz Pz Iz Oz'); 

%% Refined Composite Multivariate Generalized Multiscale Fuzzy Entropy (RCmvMFE)
% WARNING: use only for low-channel count (e.g. <10 EEG channels) or on few
% PCA components, or very short segments, otherwise, extremely long to
% compute. 

EEG = ascent_compute(EEG, 'measure', 'RCmvMFE', 'chanlist', 'Fz FCz FPz Cz Pz Iz Oz', ...
    'coarsing', 'var', ...     % 'median' (default) 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 50, ...       % number of scale factors to compute (default = 20; range = 5-100 depending on sample rate)
    'n', 2, ...                % fuzzy power (default = 2)
    'parallel', true, 'progress', true);

