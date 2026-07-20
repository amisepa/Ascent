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

% Launch EEGLAB
eeglab; close;

% if you haven'installed the plugin yet, either go to File > Manage EEGLAB
% extensions > search ascent_compute > Install
% or clone the github directory in the EEGLAB plugins folder
% or donwload the github repo and unzip it in the EEGLAB plugins folder

% Load provided sample 64-channel Biosemi dataset from the tutorial directory 
% (~1 minute of resting state eyes-closed).
% Source: 
%   Cannard, C., Wahbeh, H., & Delorme, A. (2021).
%   Validating the wearable MUSE headset for EEG spectral analysis and Frontal 
%   Alpha Asymmetry. IEEE International Conference on Bioinformatics
%   and Biomedicine (BIBM) (pp. 3603-3610). 

pluginPath = fileparts(which('eegplugin_ascent.m'));
cd(pluginPath)
EEG = pop_loadset('filename','ascent_sample_data.set','filepath',fullfile(pluginPath));

% eegplot(EEG.data, 'srate', EEG.srate, 'spacing', 80, 'winlength', 20);  % inspect the data (if needed)


fprintf("\n\n");
disp('-----------------------------------------------------------------------------------------')
disp("Are you ready for your 1st ascent?!")
disp("Run each cell below one by one by clicking in the cell and then pressing CTRL + ENTER")
disp('-----------------------------------------------------------------------------------------')
fprintf("\n");

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

EEG = ascent_compute(EEG, 'measure', 'ExSEnt', 'r', .15);

%% Fractal dimension

EEG = ascent_compute(EEG, 'measure', 'HigFracDim');  % 'HigFracDim' = Higushi version; 'FracDim' = box-counting version

%% Aperiodic exponent (i.e. slope) and offset

% normal mode
EEG = ascent_compute(EEG, 'measure', 'Aperiodic');

% time-resolved mode
EEG = ascent_compute(EEG, 'measure', 'Aperiodic', 'timeResolved', true,...
    'slidWinSec', 2, 'slidOverlap', .5, 'slidnAvg', 5);

%% Multiscale entropy (MSE)

EEG = ascent_compute(EEG, 'measure', 'MSE', ...
    'coarsing', 'mean', ...     % 'median' 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 20, ...       % number of scale factors to compute (default = 20; range = 5-100 depending on sample rate)
    'zNorm', 1, ...          % per-channel z-normalization
    'parallel', false, 'progress', true);


%% Modified MSE (mMSE)

EEG = ascent_compute(EEG, 'measure', 'mMSE', ...
    'coarsing', 'mean', ... % 'median' 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 30, ...
    'filter_mode', 'none', ...  %  'narrowband' (default), 'none' (normal MSE)
    'parallel', true, 'progress', true);

% ascent_plot(EEG.ascent.mMSE.data, EEG.chanlocs,'mMSE',EEG.ascent.mMSE.scales)


% Modified MSE (mMSE) - Time-resolved mode (on contninuous resting-state
% data)
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
    'coarsing', 'sd', ...     % 'median' 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 50, ...       % number of scale factors to compute (default = 20; range = 2-100 depending on sample rate)
    'n', 2, ...                 % fuzzy power (default = 2)
    'parallel', true, 'progress', true);

clustThresh = .5; % cluster percentile threshold for visualization (default 0.75)
ascent_plot(EEG.ascent.RCMFE.data, EEG.chanlocs,'RCMFE', EEG.ascent.RCMFE.scales, ...
    'ClusterThresh', clustThresh)

%% Multivariate Fuzzy entropy (mvFuzzEn)
% WARNING: use only for low-channel count (e.g. <10 EEG channels) or on few
% PCA components, ROIs, or very short segments, otherwise, can be extremely
% long to compute. 

% EEG = ascent_compute(EEG, 'measure', 'mvFuzzEn');  % all channels

EEG = ascent_compute(EEG, 'measure', 'mvFuzzEn', 'chanlist', 'Fpz Fz Cz Pz Iz Oz POz'); 

%% Refined Composite Multivariate Generalized Multiscale Fuzzy Entropy (RCmvMFE)
% WARNING: scale 1 is much longer than the rest of the scales.

EEG = ascent_compute(EEG, 'measure', 'RCmvMFE', 'chanlist', 'Fpz Fz Cz Pz Iz Oz POz', ...
    'coarsing', 'std', ...     % 'median' 'mean' 'trimmed mean' 'std' (default) 'var'
    'num_scales', 10, ...       % number of scale factors to compute (default = 20; range = 5-100 depending on sample rate)
    'n', 2, ...                % fuzzy power (default = 2)
    'parallel', true, 'progress', true);

%% Any measure on independent components (ICs) instead of scalp channels.
% Note: ICA was precomputed on the sample dataset using the Picard algorithm.
%
% Pass 'domain','ica' to run on EEG.icaact instead of EEG.data. Everything else
% is identical: results land in EEG.ascent.<measure> (with .domain = 'ica' and
% the mixing matrix alongside), and topographies are back-projected to the scalp
% via EEG.icawinv. The channel selection does not apply -- all ICs are used.

% % Compute ICA with Picard algorithm and effective rank reduction (if not
% % already done)
% dataRank = sum(eig(cov(double(EEG.data'))) > 1e-7);
% EEG = pop_runica(EEG, 'icatype', 'picard', 'mode', 'standard', ...
%     'maxiter', 500, 'pca', dataRank);
% % ascent_compute reconstructs EEG.icaact automatically when it is empty
% % (it is often not saved to disk, to keep .set files small).

EEG = ascent_compute(EEG, 'measure', 'RCMFE', 'domain', 'ica', ...
    'coarsing', 'sd', 'num_scales', 50);

% Outputs, now saved rather than only drawn:
% EEG.ascent.RCMFE.data     % [ICs x scales]
% EEG.ascent.RCMFE.domain   % 'ica'
% EEG.ascent.RCMFE.icawinv  % mixing matrix used for the back-projected topos

%% Run DIPFIT on ICA components, to examine dipole sources on the
% back-projected entropy topography.
%
% Prerequisites: EEGLAB open, ICA already run on EEG (EEG.icaweights etc.)
% Once EEG.dipfit exists, ascent_compute forwards it to the plots automatically
% and the region/dipole figure appears alongside the entropy figure.

% DIPFIT setup + fitting
dipfitPath = fileparts(which('dipfitdefs'));
EEG = pop_dipfit_settings(EEG, 'hdmfile', fullfile(dipfitPath,'standard_BEM','standard_vol.mat'), ...
    'coordformat','MNI','mrifile',fullfile(dipfitPath,'standard_BEM','standard_mri.mat'), ...
    'chanfile',fullfile(dipfitPath,'standard_BEM','elec','standard_1005.elc'),'chansel',1:EEG.nbchan);
EEG = pop_multifit(EEG, 1:size(EEG.icaweights,1), 'threshold',100, 'dipplot','off'); % include all ICs (the plot does the classic selection with RV<15)
EEG = pop_saveset(EEG, 'filepath', EEG.filepath, 'filename', EEG.filename);

% Same call as above; the dipole/region figure is added now that EEG.dipfit exists
EEG = ascent_compute(EEG, 'measure', 'RCMFE', 'domain', 'ica', ...
    'coarsing', 'sd', 'num_scales', 50);


%% Aperiodic parameterization on independent components (ICs) instead of
% scalp channel signals. note: ICA was precomputed on the sample dataset.

EEG = ascent_compute(EEG, 'measure', 'Aperiodic', 'domain', 'ica');

% The bottom row of the figure shows alpha power for the top IC (steepest 1/f)
% before and after removing the aperiodic component, both in µV²/Hz on a shared
% colour scale. Those values are saved, not re-derived by the plot:
% EEG.ascent.Aperiodic.data.alpha_raw   % total alpha power   [ICs x 1]
% EEG.ascent.Aperiodic.data.alpha_osc   % alpha above the 1/f [ICs x 1]

% Time-resolved mode
EEG = ascent_compute(EEG, 'measure', 'Aperiodic', 'domain', 'ica', ...
    'timeResolved', true, 'slidWinSec', 2, 'slidOverlap', .5, 'slidnAvg', 5);



