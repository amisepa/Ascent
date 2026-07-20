% ascent_group_analysis - Group-level statistical analysis script for the ASCENT EEGLAB plugin.
%
% This script preprocesses 64-channel EEG data (Biosemi) for two conditions 
%(eyes-closed and eyes-open), computes a broad range of entropy, complexity, 
% and aperiodic measures using the ASCENT plugin, and performs nonparametric 
% group-level statistical comparisons with multiple comparison correction (MCC). 
% Results are saved as figures and CSV summary tables.
%
% USAGE:
%   Run section by section or as a full script after setting data_path and
%   parameters at the top. Preprocessed .set files must already exist in
%   the eyes_closed/ and eyes_open/ subdirectories of data_path.
%
% WORKFLOW:
%   1. Compute all ASCENT measures per subject for each condition
%   2. Save outputs to .mat files per condition
%   3. Reload, reorganize electrode order (front-to-back), and align conditions
%   4. Run nonparametric paired permutation t-tests per measure
%   5. Apply MCC (cluster-based or TFCE) and visualize significant effects
%   6. Export figures (.fig, .png) and cluster summary tables (.csv)
%
% KEY PARAMETERS (set at top of script):
%   num_scales  - number of scale factors for multiscale measures (default: 30)
%   num_chan    - number of EEG channels (default: 64)
%   coarsing    - coarse-graining operator: 'sd' (standard deviation, default),
%                 'mean', 'median', 'var' (see ASCENT documentation)
%   nPerm       - number of permutations for null distribution (default: 5000)
%   alpha       - significance threshold (default: 0.01)
%   ct          - central tendency for test statistic: 'mean' (standard
%                 paired t-test) or 'trimmed mean' (Yuen t-test)
%   grp_type    - 'dpt' (dependent/paired) or 'idpt' (independent)
%   mcc_type    - MCC method: 0=uncorrected, 1=t-max, 2=cluster, 3=TFCE
%
% MEASURES COMPUTED:
%   Single-scale : SampEn, FuzzEn, ExSEnt (HD, HA, HDA), HigFracDim,
%                  Aperiodic exponent and offset, raw and corrected PSD
%   Multiscale   : MSE, mMSE, MFE, RCMFE
%   (Multivariate measures RCmvMFE are available but commented out to
%    reduce computation time and they are not directly comparable to the other measures)
%
% REQUIREMENTS:
%   - MATLAB R2019b or later (for tiledlayout, arguments blocks, pyenv)
%   - EEGLAB (Delorme & Makeig, 2004): https://sccn.ucsd.edu/eeglab
%       Plugins required: clean_rawdata, ICLabel, Picard (for preprocessing)
%   - ASCENT EEGLAB plugin (Cannard & Delorme, 2025): https://github.com/amisepa/ascent
%       Provides: ascent_compute, ascent_plot, run_stats_permutation,
%                 compute_mcc, pull_clusters, plot_results, plot_clusters
%   - EEG Robust Statistics toolbox (Cannard, 2025): https://github.com/amisepa/eeg_robust_statistics
%       Provides: run_stats_permutation, compute_mcc, pull_clusters,
%                 plot_results, plot_clusters, progressbar
%   - No additional MATLAB toolboxes required (Parallel Computing Toolbox
%     optional for faster entropy computation via ascent_compute)
%
% OUTPUTS:
%   ascent_outputs_EC_<coarsing>_new.mat  - EC condition measures (all subjects)
%   ascent_outputs_EO_<coarsing>_new.mat  - EO condition measures (all subjects)
%   figures_new2/group_results/           - figures (.fig, .png) and cluster
%                                           summary tables (.csv) per measure
%
% NOTES:
%   - Scale s=1 is excluded from all multiscale outputs (MSE, mMSE, MFE);
%     RCMFE(:,1,:) is additionally removed before statistics to align scales.
%   - mMSE scale labels reflect bandpass frequency bounds (Hz) derived from
%     the filt-skip scheme (Kosciessa et al., 2020).
%   - Statistics are run on linear PSD values; dB-normalized alternatives
%     are available as commented-out lines.
%   - The aperiodic correction is applied in fixed mode only (no knee).
%
% REFERENCES:
%   Cannard, C., & Delorme, A. (2025). Introducing the ASCENT EEGLAB plugin:
%     Aperiodic Signal Complexity Estimation for Neurophysiological Time series.
%     Entropy.
%
%   Delorme, A., & Makeig, S. (2004). EEGLAB: An open source toolbox for
%     analysis of single-trial EEG dynamics including independent component
%     analysis. Journal of Neuroscience Methods, 134, 9-21.
%
%   Groppe, D.M., Urbach, T.P., & Kutas, M. (2011). Mass univariate analysis
%     of event-related brain potentials/fields I: A critical tutorial review.
%     Psychophysiology, 48, 1711-1725.
%
%   Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of
%     EEG- and MEG-data. Journal of Neuroscience Methods, 164, 177-190.
%
%   Pernet, C.R., et al. (2015). Cluster-based computational methods for mass
%     univariate analyses of event-related brain potentials/fields: A simulation
%     study. Journal of Neuroscience Methods, 250, 85-93.
%
%   Smith, S.M., & Nichols, T.E. (2009). Threshold-free cluster enhancement:
%     Addressing problems of smoothing, threshold dependence and localisation
%     in cluster inference. NeuroImage, 44, 83-98.
%
% Author  : Cedric Cannard, 2021 - ASCENT EEGLAB Plugin
% Contact : ccannard@pm.me
% GitHub  : https://github.com/amisepa/ascent

clear; close all; clc
data_path = 'C:\Users\ccann\Documents\biosemi_data';
pluginPath = 'C:\Users\ccann\Documents\MATLAB\Ascent';
addpath(fullfile(pluginPath, 'functions'))
cd(pluginPath)
eeglab; close


%% Compute on whole group - Eyes closed (EC) condition

coarsing 	= 'sd'; 	% 'mean' and 'sd' were run for the paper

%------ PARAMETERS -------------------------
num_scales 	= 30; 		% 30 scales were run for the paper
num_chan 	= 64;
% ------------------------------------------

cd(fullfile(data_path, 'eyes_closed'))
filenames = {dir('*.bdf').name}';
num_files = length(filenames);

SampEn = nan(num_chan, num_files);
ExSEnt1 = nan(num_chan, num_files);
ExSEnt2 = nan(num_chan, num_files);
ExSEnt3 = nan(num_chan, num_files);
FuzzEn = nan(num_chan, num_files);
FracDim = nan(num_chan, num_files);
Exponent = nan(num_chan, num_files);
Offset = nan(num_chan, num_files);
PSD = nan(num_chan, 79, num_files);
PSD_corr = nan(num_chan, 79, num_files);
MSE = nan(num_chan, num_scales-1, num_files);
mMSE = nan(num_chan, num_scales-1, num_files);
MFE = nan(num_chan, num_scales-1, num_files);
RCMFE = nan(num_chan, num_scales, num_files);
% RCmvMFE = nan(num_chan, num_scales, num_files);
progressbar('Prosessing and computing measures for eyes-closed condition')
for iFile = 1:num_files

    disp('')
    fprintf('--------------------------------------------------------\n')
    fprintf('                        SUBJECT %g/%g \n', iFile, num_files)
    fprintf('--------------------------------------------------------\n')
    disp('')

    progressbar(iFile / num_files);

    % EEG = pop_biosig(filenames{iFile});
    % EEG = pop_resample(EEG, 256);
    % EEG = pop_select( EEG, 'chantype',{'EEG'});
    % EEG = pop_chanedit(EEG, {'lookup','standard_1005.elc'});
    % EEG = reorder_channels(EEG);
    % EEG = ref_infinity(EEG);
    % EEG = pop_eegfiltnew(EEG, 'locutoff',1);
    % % EEG = pop_eegfiltnew(EEG, 'locutoff',59,'hicutoff',61,'revfilt',1);
    % EEG = pop_eegfiltnew(EEG,'hicutoff',50);
    % oriEEG = EEG;
    % EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.75, ...
    %     'LineNoiseCriterion',5,'Highpass','off','BurstCriterion',40, ...
    %     'WindowCriterion','off','BurstRejection','off','Distance','Euclidian', ...
    %     'WindowCriterionTolerances','off');
    % % vis_artifacts(EEG,oriEEG);
    % EEG = pop_interp(EEG, oriEEG.chanlocs, 'spherical');
    % dataRank = sum(eig(cov(double(EEG.data'))) > 1e-7);
    % EEG = pop_runica(EEG, 'icatype', 'picard', 'mode', 'standard', ...
    %     'maxiter', 500, 'pca', dataRank);
    % EEG = pop_iclabel(EEG, 'default');
    % EEG = pop_icflag(EEG, [NaN NaN;0.9 1;0.9 1;0.95 1;0.95 1;0.9 1;NaN NaN]);
    % % pop_selectcomps(EEG, 1:24);
    % EEG = pop_subcomp(EEG, [], 0);
    % % pop_eegplot(EEG,1,1,1);
    
    % Load if already preprocessed
    EEG = pop_loadset('filepath', fullfile(data_path, 'eyes_closed'), 'filename', sprintf('%s.set', filenames{iFile}(1:end-4)));

    % Compute complexity measures
    EEG = ascent_compute(EEG, 'measure', 'SampEn', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'ExSEnt', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'FuzzEn', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'HigFracDim', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'Aperiodic', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'MSE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'mMSE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'MFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    % EEG = ascent_compute(EEG, 'measure', 'CMFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'RCMFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    % EEG = ascent_compute(EEG, 'measure', 'RCmvMFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);

    % ascent_plot(EEG.ascent.MSE.data, EEG.chanlocs, 'mse', EEG.ascent.MSE.scales);
    % ascent_plot(EEG.ascent.MFE.data, EEG.chanlocs, 'RCMFE', EEG.ascent.MFE.scales);
    % ascent_plot(EEG.ascent.RCMFE.data, EEG.chanlocs, 'RCMFE', EEG.ascent.RCMFE.scales);
    % ascent_plot(EEG.ascent.ExSEnt.data, EEG.chanlocs, 'exsent', []);
    
    % Save for running entropy on ICA time series later
    % EEG = pop_saveset(EEG, 'filepath', fullfile(data_path, 'eyes_closed'), 'filename', sprintf('%s.set', filenames{iFile}(1:end-4)));

    SampEn(:,iFile) = EEG.ascent.SampEn.data;
    ExSEnt1(:,iFile) = EEG.ascent.ExSEnt.data.HD;
    ExSEnt2(:,iFile) = EEG.ascent.ExSEnt.data.HA;
    ExSEnt3(:,iFile) = EEG.ascent.ExSEnt.data.HDA;
    FuzzEn(:,iFile) = EEG.ascent.FuzzEn.data;
    FracDim(:,iFile) = EEG.ascent.HigFracDim.data;
    Exponent(:,iFile) = EEG.ascent.Aperiodic.data.exponent;
    Offset(:,iFile) = EEG.ascent.Aperiodic.data.offset;
    PSD(:,:,iFile) = EEG.ascent.Aperiodic.data.psd;
    PSD_corr(:,:,iFile) = EEG.ascent.Aperiodic.data.psd_corrected;
    freqs = EEG.ascent.Aperiodic.data.freqs;
    MSE(:,:,iFile) = EEG.ascent.MSE.data;
    mMSE(:,:,iFile) = EEG.ascent.mMSE.data;
    MFE(:,:,iFile) = EEG.ascent.MFE.data;
    % CMFE(:,:,iFile) = EEG.ascent.CMFE.data;
    RCMFE(:,:,iFile) = EEG.ascent.RCMFE.data;
    % RCmvMFE(:,:,iFile) = EEG.ascent.RCmvMFE.data;

    chanlocs = EEG.chanlocs;
    scales = EEG.ascent.MSE.scales;
    scales_bounds = EEG.ascent.mMSE.scales;

    save(fullfile(data_path, sprintf('ascent_outputs_EC_%s_new.mat', coarsing)), ...
        "SampEn", "ExSEnt1", "ExSEnt2", "ExSEnt3", "FuzzEn", "FracDim", ...
        "Exponent", "Offset", "PSD", "PSD_corr", "freqs", ...
        "MSE", "MFE", "mMSE", "RCMFE", ...
        'chanlocs', 'scales', 'scales_bounds')

end
disp("Done computing on the whole group for eyes-closed condition!")

% Check there are no subjects with NaNs
nan_subject = any(squeeze(any(isnan(MFE(:,2:end,:)), 1)));  % [50 x 40] logical (scale x subject)
if any(nan_subject)
    error("some subjects have NaNs! Check data!")
    disp(nan_subject)
end


% Compute on whole group - Eyes opened (EO) condition
cd(fullfile(data_path, 'eyes_open'))
filenames = {dir('*.bdf').name}';
num_files = length(filenames);

SampEn = nan(num_chan, num_files);
ExSEnt1 = nan(num_chan, num_files);
ExSEnt2 = nan(num_chan, num_files);
ExSEnt3 = nan(num_chan, num_files);
FuzzEn = nan(num_chan, num_files);
FracDim = nan(num_chan, num_files);
Exponent = nan(num_chan, num_files);
Offset = nan(num_chan, num_files);
PSD = nan(num_chan, 79, num_files);
PSD_corr = nan(num_chan, 79, num_files);
MSE = nan(num_chan, num_scales-1, num_files);
mMSE = nan(num_chan, num_scales-1, num_files);
MFE = nan(num_chan, num_scales-1, num_files);
% CMFE = nan(num_chan, num_scales, num_files);
RCMFE = nan(num_chan, num_scales, num_files);
% RCmvMFE = nan(num_chan, num_scales, num_files);
progressbar('Prosessing and computing measures for eyes-open condition')
for iFile = 1:num_files
    disp('')
    fprintf('--------------------------------------------------------\n')
    fprintf('                        SUBJECT %g/%g \n', iFile, num_files)
    fprintf('--------------------------------------------------------\n')
    disp('')

    progressbar(iFile / num_files);

    % EEG = pop_biosig(filenames{iFile});
    % EEG = pop_resample(EEG, 256);
    % EEG = pop_select( EEG, 'chantype',{'EEG'});
    % EEG = pop_chanedit(EEG, {'lookup','standard_1005.elc'});
    % EEG = reorder_channels(EEG);
    % EEG = ref_infinity(EEG);
    % EEG = pop_eegfiltnew(EEG, 'locutoff',1);
    % EEG = pop_eegfiltnew(EEG, 'locutoff',59,'hicutoff',61,'revfilt',1);
    % % EEG = pop_eegfiltnew(EEG,'hicutoff',50);
    % oriEEG = EEG;
    % EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.75, ...
    %     'LineNoiseCriterion',5,'Highpass','off','BurstCriterion',40, ...
    %     'WindowCriterion','off','BurstRejection','off','Distance','Euclidian', ...
    %     'WindowCriterionTolerances','off');
    % % vis_artifacts(EEG,oriEEG);
    % EEG = pop_interp(EEG, oriEEG.chanlocs, 'spherical');
    % dataRank = sum(eig(cov(double(EEG.data'))) > 1e-7);
    % EEG = pop_runica(EEG, 'icatype', 'picard', 'mode', 'standard','pca', dataRank);
    % EEG = pop_iclabel(EEG, 'default');
    % EEG = pop_icflag(EEG, [NaN NaN;0.9 1;0.9 1;0.95 1;0.95 1;0.9 1;NaN NaN]);
    % % pop_selectcomps(EEG, 1:24);
    % EEG = pop_subcomp(EEG, [], 0);
    % % pop_eegplot(EEG,1,1,1);

    % % Load if already processed
    EEG = pop_loadset('filepath', fullfile(data_path, 'eyes_open'), 'filename', sprintf('%s.set', filenames{iFile}(1:end-4)));

    % Compute complexity measures
    EEG = ascent_compute(EEG, 'measure', 'SampEn', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'ExSEnt', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'FuzzEn', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'HigFracDim', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'Aperiodic', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'MSE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'mMSE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'MFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    % EEG = ascent_compute(EEG, 'measure', 'CMFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'RCMFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    % EEG = ascent_compute(EEG, 'measure', 'RCmvMFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    % 
    % % Save for running entropy on ICA time series later
    % EEG = pop_saveset(EEG, 'filepath', fullfile(data_path, 'eyes_open'), 'filename', sprintf('%s.set', filenames{iFile}(1:end-4)));

    SampEn(:,iFile) = EEG.ascent.SampEn.data;
    ExSEnt1(:,iFile) = EEG.ascent.ExSEnt.data.HD;
    ExSEnt2(:,iFile) = EEG.ascent.ExSEnt.data.HA;
    ExSEnt3(:,iFile) = EEG.ascent.ExSEnt.data.HDA;
    FuzzEn(:,iFile) = EEG.ascent.FuzzEn.data;
    FracDim(:,iFile) = EEG.ascent.HigFracDim.data;
    Exponent(:,iFile) = EEG.ascent.Aperiodic.data.exponent;
    Offset(:,iFile) = EEG.ascent.Aperiodic.data.offset;
    PSD(:,:,iFile) = EEG.ascent.Aperiodic.data.psd;
    PSD_corr(:,:,iFile) = EEG.ascent.Aperiodic.data.psd_corrected;
    freqs = EEG.ascent.Aperiodic.data.freqs;
    MSE(:,:,iFile) = EEG.ascent.MSE.data;
    mMSE(:,:,iFile) = EEG.ascent.mMSE.data;
    MFE(:,:,iFile) = EEG.ascent.MFE.data;
    % CMFE(:,:,iFile) = EEG.ascent.CMFE.data;
    RCMFE(:,:,iFile) = EEG.ascent.RCMFE.data;
    % RCmvMFE(:,:,iFile) = EEG.ascent.RCmvMFE.data;

    chanlocs = EEG.chanlocs;
    scales = EEG.ascent.MSE.scales;
    scales_bounds = EEG.ascent.mMSE.scales;

    save(fullfile(data_path, sprintf('ascent_outputs_EO_%s_new.mat', coarsing)), ...
        "SampEn", "ExSEnt1", "ExSEnt2", "ExSEnt3", "FuzzEn", "FracDim", ...
        "Exponent", "Offset", "PSD", "PSD_corr", "freqs", ...
        "MSE", "MFE", "mMSE", "RCMFE", ...
        'chanlocs', 'scales', 'scales_bounds')

end
disp("Done computing on the whole group for eyes-open condition!")

% Check there are no subjects with NaNs
nan_subject = any(squeeze(any(isnan(MFE(:,2:end,:)), 1)));  % [50 x 40] logical (scale x subject)
if any(nan_subject)
    error("some subjects have NaNs! Check data!")
    disp(nan_subject)
end

gong

%% Load, separate the data by condition, & reorganize electrodes

coarsing = 'sd';

% Load chanlocs only
load(fullfile(data_path, sprintf('ascent_outputs_EC_%s_new.mat', coarsing)), 'chanlocs')
cd(data_path)

% Front-to-back target order
desired_order = {
    'Fp1','FPz','FP2', ...
    'AF7','AF3','AFz','AF4','AF8', ...
    'F7','F5','F3','F1','Fz','F2','F4','F6','F8', ...
    'FT7','FT8', ...
    'FC5','FC3','FC1','FCz','FC2','FC4','FC6', ...
    'T7','T8', ...
    'C5','C3','C1','Cz','C2','C4','C6', ...
    'TP7','TP8', ...
    'CP5','CP3','CP1','CPz','CP2','CP4','CP6', ...
    'P9','P7','P5','P3','P1','Pz','P2','P4','P6','P8','P10', ...
    'PO7','PO3','POz','PO4','PO8', ...
    'O1','Oz','Iz','O2'
    };

% Map current labels -> desired positions, then sort
current_labels = {chanlocs.labels};
[found, loc_in_desired] = ismember(lower(current_labels), lower(desired_order));
if ~all(found)
    missing = current_labels(~found);
    error('Missing in desired_order: %s', strjoin(missing, ', '));
end
[~, order_idx] = sort(loc_in_desired, 'ascend');

% Reorder chanlocs
chanlocs = chanlocs(order_idx);
% save(fullfile(data_path, 'chanlocs_reorganized.mat'), 'chanlocs')
disp({chanlocs.labels}')

% Load EC metrics and reorder
load(fullfile(data_path, sprintf('ascent_outputs_EC_%s_new.mat', coarsing)), ...
    'SampEn','FuzzEn','ExSEnt1','ExSEnt2','ExSEnt3','FracDim','Exponent', ...
    'Offset','PSD','PSD_corr','MSE','MFE','mMSE','RCMFE','scales','scales_bounds')  % don't reload chanlocs here!!
load(fullfile(data_path, sprintf('ascent_outputs_EC_%s_new.mat', coarsing)))
SampEn1  = SampEn(order_idx, :);
FuzzEn1  = FuzzEn(order_idx, :);
ExSEnt1_1  = ExSEnt1(order_idx, :);
ExSEnt1_2  = ExSEnt2(order_idx, :);
ExSEnt1_3  = ExSEnt3(order_idx, :);
FracDim1 = FracDim(order_idx, :);
Exponent1 = Exponent(order_idx,:);
Offset1 = Offset(order_idx,:);
PSD1  = PSD(order_idx, :, :);
PSD_corr1 = PSD_corr(order_idx, :, :);
MSE1     = MSE(order_idx, :, :);
MFE1     = MFE(order_idx, :, :);
mMSE1    = mMSE(order_idx, :, :);
RCMFE1   = RCMFE(order_idx, :, :);
% RCmvMFE1   = RCmvMFE(order_idx, :, :);

% Load EO metrics and reorder
load(fullfile(data_path, sprintf('ascent_outputs_EO_%s_new.mat', coarsing)), ...
    'SampEn','FuzzEn','ExSEnt1','ExSEnt2','ExSEnt3','FracDim','Exponent', ...
    'Offset','PSD','PSD_corr','freqs','MSE','MFE','mMSE','RCMFE','scales', 'scales_bounds') % don't reload chanlocs here!!
% load(fullfile(data_path, sprintf('ascent_outputs_EO_%s_new.mat', coarsing)))
SampEn2  = SampEn(order_idx, :);
FuzzEn2  = FuzzEn(order_idx, :);
ExSEnt2_1  = ExSEnt1(order_idx, :);
ExSEnt2_2  = ExSEnt2(order_idx, :);
ExSEnt2_3  = ExSEnt3(order_idx, :);
FracDim2 = FracDim(order_idx, :);
Exponent2 = Exponent(order_idx,:);
Offset2 = Offset(order_idx,:);
PSD2 = PSD(order_idx, :, :);
PSD_corr2 = PSD_corr(order_idx, :, :);
MSE2     = MSE(order_idx, :, :);
MFE2     = MFE(order_idx, :, :);
mMSE2    = mMSE(order_idx, :, :);
RCMFE2   = RCMFE(order_idx, :, :);
% RCmvMFE2   = RCmvMFE(order_idx, :, :);

%% --------- PARAMETERS ---------------------------------------

nPerm       = 2000;       % number of permutations for H0
alpha       = 0.05;       % NaN to show maps of p-values
ct          = 'mean';        % central tendency method ('mean' for normal t-test, 'trimmed mean' for Yuen t-test)
grp_type    = 'dpt';   % groupe is dependent ('dpt') or independent ('idpt')
mcc_type    = 2;       % 0 = uncorrected; 1 = t-max; 2 = cluster; 3 = TFCE

outputs_path = fullfile(pluginPath, 'figures', 'group_results'); mkdir(outputs_path)

% load(fullfile(data_path,'chanlocs_reorganized.mat'), 'chanlocs')
% load colormap_bwr.mat; dmap(1,:) = [0.9 0.9 0.9]; % set NaNs to grey
rdbu_cmap   = interp1([1;128;256],[0.698 0.094 0.168;1 1 1;0.129 0.400 0.675],(1:256)','linear');
dmap        = flipud(max(0,min(1,rdbu_cmap)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% num_scales = 1:scales(end);

%% --------- UNISCALES ---------------------------------------

figure('Name', 'Results: Uniscale measures', 'Color','w','ToolBar','none','MenuBar','none');

% SampEn
disp("------------------------------------------------------")
disp("         MEASURE: SampEn")
% nexttile
subplot(2,4,1)
% figure('color','w');
% subplot(2,1,1); title('Mean across channels')
% load('colorbrewer.mat'); colormap = colorbrewer.qual.Set3{12}./255; colormap(8,:) = [];
% rm_raincloud({trimmean(SampEn1,20,1) trimmean(SampEn2,20,1)}, colormap([5 4],:));
% legend('EC','','EO',''); xlabel('SampEn') % title('SampEn')
% subplot(2,1,2); title('Mean across subjects')
% rm_raincloud({trimmean(SampEn1,20,2) trimmean(SampEn2,20,2)}, colormap([5 4],:));
% legend('EC','','EO',''); xlabel('SampEn') % title('SampEn')
% set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(SampEn1, SampEn2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(SampEn1, SampEn2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, [], chanlocs, ...
        'nonlinear', grp_type, {size(SampEn1,2) size(SampEn2,2)}, [], [], [], 'g');
    % tvals(~mask) = 0;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("SampEn (coarse: %s)", coarsing));
    title("SampEn")
end
disp("")


% FuzzEn
disp("------------------------------------------------------")
disp("         MEASURE: FuzzEn")
% nexttile
subplot(2,4,2)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(FuzzEn1, FuzzEn2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(FuzzEn1, FuzzEn2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, [], chanlocs, ...
        'nonlinear', grp_type, {size(SampEn1,2) size(SampEn2,2)}, [], [], [], 'g');
    % tvals(~mask) = 0;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("FuzzEn (coarse: %s)", coarsing));
    title("FuzzEn")
end
disp("")

% ExSEnt1
disp("------------------------------------------------------")
disp("         MEASURE: ExSent (Duration)")
% nexttile
subplot(2,4,3)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(ExSEnt1_1, ExSEnt2_1, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(ExSEnt1_1, ExSEnt2_1, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, [], chanlocs, ...
        'nonlinear', grp_type, {size(SampEn1,2) size(SampEn2,2)}, [], [], [], 'g');
    % tvals(~mask) = 0;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("ExSEnt (coarse: %s)", coarsing));
    title("ExSEnt (duration)")
end
disp("")

% ExSEnt2
disp("------------------------------------------------------")
disp("         MEASURE: ExSEnt (Amplitude)")
% nexttile
subplot(2,4,4)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(ExSEnt1_2, ExSEnt2_2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(ExSEnt1_2, ExSEnt2_2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, [], chanlocs, ...
        'nonlinear', grp_type, {size(SampEn1,2) size(SampEn2,2)}, [], [], [], 'g');
    % tvals(~mask) = NaN;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("ExSEnt (coarse: %s)", coarsing));
    title("ExSEnt (amplitude)")
end
disp("")

% ExSEnt3
disp("------------------------------------------------------")
disp("         MEASURE: ExSEnt (Duration + Amplitude)")
% nexttile
subplot(2,4,5)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(ExSEnt1_3, ExSEnt2_3, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(ExSEnt1_3, ExSEnt2_3, nPerm, ct, grp_type);
% mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% % [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, [], chanlocs, ...
        'nonlinear', grp_type, {size(SampEn1,2) size(SampEn2,2)}, [], [], [], 'g');
    % tvals(~mask) = NaN;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("ExSEnt (coarse: %s)", coarsing));
    title("ExSEnt (amp + dur)")
end
disp("")

% HigFracDim
disp("------------------------------------------------------")
disp("         MEASURE: Higuchi Fractal Dimension")
% nexttile
subplot(2,4,6)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(FracDim1, FracDim2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(FracDim1, FracDim2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, [], chanlocs, ...
        'nonlinear', grp_type, {size(SampEn1,2) size(SampEn2,2)}, [], [], [], 'g');
    % tvals(~mask) = 0;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("FracDim (coarse: %s)", coarsing));
    title("HigFracDim")
end
disp("")

% Aperiodic Exponent
disp("------------------------------------------------------")
disp("         MEASURE: Aperiodic Exponent")
% nexttile
subplot(2,4,7)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(Exponent1, Exponent2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(Exponent1, Exponent2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, [], chanlocs, ...
        'nonlinear', grp_type, {size(SampEn1,2) size(SampEn2,2)}, [], [], [], 'g');
    % tvals(~mask) = 0;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("FracDim (coarse: %s)", coarsing));
    title("Aperiodic Exponent")
end
disp("")

% Aperiodic Offset
disp("------------------------------------------------------")
disp("         MEASURE: Aperiodic Offset")
% nexttile
subplot(2,4,8)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(Offset1, Offset2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(Offset1, Offset2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, [], chanlocs, ...
        'nonlinear', grp_type, {size(SampEn1,2) size(SampEn2,2)}, [], [], [], 'g');
    % tvals(~mask) = 0;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("FracDim (coarse: %s)", coarsing));
    title("Aperiodic Offset")
end
disp("")


% set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');
set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold', 'CLim', [-10 10]);
saveas(gcf, fullfile(outputs_path, 'fig_uniscales.fig'));
print(gcf, fullfile(outputs_path, 'fig_uniscales.png'), '-dpng', '-r300');

%% ---------- CORRELATION ANALYSIS ON EC-EO DIFFERENCES ----------
% This section computes the difference (EC minus EO) for each single-scale
% measure per subject, then calculates:
%   1. A global correlation matrix (channel-averaged differences)
%   2. A spatial correlation matrix of the t-maps (EC vs EO)
% Both are plotted with the circle-based style (color + size + r-value),
% including the diagonal (r=1).

% Compute difference scores (EC - EO) for each measure
% All matrices are [nChans x nSubj]
SampEn_diff  = SampEn1  - SampEn2;
FuzzEn_diff  = FuzzEn1  - FuzzEn2;
ExSEnt1_diff = ExSEnt1_1 - ExSEnt2_1;
ExSEnt2_diff = ExSEnt1_2 - ExSEnt2_2;
ExSEnt3_diff = ExSEnt1_3 - ExSEnt2_3;
FracDim_diff = FracDim1 - FracDim2;
Exponent_diff= Exponent1 - Exponent2;
Offset_diff  = Offset1  - Offset2;

% List of measure names and corresponding difference matrices
measureNames = {'SampEn','FuzzEn','ExSEnt (HD)','ExSEnt (HA)','ExSEnt (HDA)',...
                'HigFracDim','Aperiodic Exponent','Aperiodic Offset'};
measureData  = {SampEn_diff, FuzzEn_diff, ExSEnt1_diff, ExSEnt2_diff, ...
                ExSEnt3_diff, FracDim_diff, Exponent_diff, Offset_diff};
nMeasures = numel(measureNames);
nChans = size(measureData{1}, 1);

% Correlation type: 'Spearman' (robust) or 'Pearson'
corrType = 'Spearman';

% ----- Global correlation matrix (channel-averaged differences) -----
% Average each difference measure across channels, then compute pairwise correlations
avgMeasures = zeros(size(measureData{1}, 2), nMeasures);  % subjects x measures
for m = 1:nMeasures
    avgMeasures(:, m) = mean(measureData{m}, 1, 'omitnan')';
end

[R_global, P_global] = corr(avgMeasures, 'Type', corrType, 'Rows', 'pairwise');

% BH-FDR correction across all unique pairs (lower triangle)
pvals_global = P_global(tril(true(nMeasures), -1));
pvals_global_fdr = bh_fdr(pvals_global);
% Rebuild FDR-corrected p-value matrix (for reference, not displayed)
P_global_fdr = ones(nMeasures);
idx = 1;
for i = 1:nMeasures
    for j = 1:i-1
        P_global_fdr(i,j) = pvals_global_fdr(idx);
        P_global_fdr(j,i) = pvals_global_fdr(idx);
        idx = idx + 1;
    end
end

% Plot global correlation matrix (circle style, with diagonal)
plot_corrmatrix_with_text(R_global, measureNames, ...
    sprintf('Global %s correlations of EC-EO differences (channel-averaged)', corrType), ...
    fullfile(outputs_path, sprintf('corr_diff_matrix_%s', corrType)));

% ----- Spatial correlation of t-maps (group-level) -----
% Recompute t-values for each measure (EC vs EO)
EC_data = {SampEn1, FuzzEn1, ExSEnt1_1, ExSEnt1_2, ExSEnt1_3, FracDim1, Exponent1, Offset1};
EO_data = {SampEn2, FuzzEn2, ExSEnt2_1, ExSEnt2_2, ExSEnt2_3, FracDim2, Exponent2, Offset2};

tvals_all = zeros(nChans, nMeasures);
for m = 1:nMeasures
    [tvals_all(:,m), ~, ~, ~] = run_stats_permutation(EC_data{m}, EO_data{m}, nPerm, ct, grp_type);
end

% Compute spatial correlation (across channels) for each pair
[R_spatial, P_spatial] = corr(tvals_all, 'Type', corrType, 'Rows', 'pairwise');

% BH-FDR correction
pvals_spatial = P_spatial(tril(true(nMeasures), -1));
pvals_spatial_fdr = bh_fdr(pvals_spatial);
P_spatial_fdr = ones(nMeasures);
idx = 1;
for i = 1:nMeasures
    for j = 1:i-1
        P_spatial_fdr(i,j) = pvals_spatial_fdr(idx);
        P_spatial_fdr(j,i) = pvals_spatial_fdr(idx);
        idx = idx + 1;
    end
end

% Plot spatial correlation matrix (circle style, with diagonal)
plot_corrmatrix_with_text(R_spatial, measureNames, ...
    'Spatial correlation of t-maps (EC vs EO)', ...
    fullfile(outputs_path, 'spatial_corr_tmaps'));

disp('Correlation analysis complete.');

%% ---------- Helper functions ----------

function plot_corrmatrix_with_text(C, labels, title_str, save_path)
% plot_corrmatrix_with_text  Visualize a correlation matrix with colored circles
% and the correlation coefficient displayed inside each circle.
% Includes the diagonal (r=1). Uses red-blue colormap, black text.
%
% Inputs
%   C         Square correlation matrix [N x N], values in [-1, 1]
%   labels    Cell array of variable names (length N)
%   title_str Figure title
%   save_path Base path for saving .fig and .png (without extension)

if size(C,1) ~= length(labels)
    error("Labels and correlation matrix must have the same number of variables: %g", size(C,1))
end

% Use the red-blue diverging colormap
cmap = rdbu_colormap();

% Keep lower triangle AND diagonal (set upper triangle to 0 to hide)
C_full = C;
C_full(triu(true(size(C)), 1)) = 0;

% Compute center of each circle
nVar = size(C_full, 1);
x = 1 : nVar;
y = 1 : nVar;
[xAll, yAll] = meshgrid(x, y);
hide = (C_full == 0) & ~eye(nVar);
xAll(hide) = nan;

% Scale colors from -1 to 1
clrLim = [-1, 1];
Cscaled = (C_full - clrLim(1)) / range(clrLim);
colIdx = discretize(Cscaled, linspace(0, 1, size(cmap,1)));

% Scale circle size
diamLim = [0.3, 0.9];
Cscaled_abs = abs(C_full);
diamSize = Cscaled_abs * range(diamLim) + diamLim(1);

% Create figure (wider to prevent cropping)
fh = figure('Color','w', 'Position', [100 100 1000 700]);  % increased width
ax = axes(fh);
hold(ax, 'on');
colormap(cmap);

% Draw axis labels (left and bottom) – BOLD
tickvalues = 1:nVar;
% Left labels: place further left to avoid cropping long names
x_left = -1 * ones(size(tickvalues));   % moved from -0.5 to -1
text(x_left, tickvalues, labels, 'HorizontalAlignment','right', ...
    'FontSize',12, 'FontWeight','bold', 'Color','k');
% Bottom labels: place at y = nVar+1, rotated
y_bottom = (nVar + 1) * ones(size(tickvalues));
text(tickvalues, y_bottom, labels, 'HorizontalAlignment','right', ...
    'Rotation',45, 'FontSize',12, 'FontWeight','bold', 'Color','k');

% Draw circles
theta = linspace(0, 2*pi, 100);
for i = 1:numel(xAll)
    if isnan(xAll(i)), continue; end
    fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
         diamSize(i)/2 * sin(theta) + yAll(i), ...
         cmap(colIdx(i),:), 'LineStyle','none');
    r_val = C_full(i);
    text(xAll(i), yAll(i), sprintf('%.2f', r_val), ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'FontSize', 9, 'Color', 'k', 'FontWeight','bold');
end

% Figure formatting
set(ax, 'YDir', 'reverse');
% Expand limits to show labels fully
xlim([-2, nVar+1]);   % left margin increased to -2
ylim([0.5, nVar+2]);  % bottom margin for rotated labels
cb = colorbar;
ylabel(cb, 'Correlation coefficient', 'FontSize', 12, 'FontWeight', 'bold');
clim(clrLim);
axis off;
set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');
title(title_str, 'FontSize', 13, 'FontWeight', 'bold');

% Save
saveas(gcf, [save_path '.fig']);
print(gcf, [save_path '.png'], '-dpng', '-r300');
end

function cmap = rdbu_colormap()
% Red-Blue diverging colormap (256 colors)
n = 256;
half = floor(n/2);
r1 = linspace(0.019, 0.97, half)';
g1 = linspace(0.188, 0.97, half)';
b1 = linspace(0.380, 0.97, half)';
r2 = linspace(0.97, 0.698, n-half)';
g2 = linspace(0.97, 0.094, n-half)';
b2 = linspace(0.97, 0.169, n-half)';
cmap = [r1 g1 b1; r2 g2 b2];
end



function p_adj = bh_fdr(p)
% Benjamini-Hochberg FDR correction
p = p(:); n = numel(p);
[ps, ord] = sort(p);
p_adj_sorted = min(1, ps .* n ./ (1:n)');
for k = n-1:-1:1
    p_adj_sorted(k) = min(p_adj_sorted(k), p_adj_sorted(k+1));
end
p_adj = nan(n,1);
p_adj(ord) = p_adj_sorted;
end



%% %%%%%%%%%%%%%%%%%%%% MULTISCALES %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coarsing = 'sd';
nPerm = 2000;          % number of permutations for H0
alpha = 0.05;           % NaN to show maps of p-values
% ct = 'mean';            % central tendency method ('mean' for normal t-test, 'trimmed mean' for Yuen t-test)
% grp_type = 'dpt';       % groupe is dependent ('dpt') or independent ('idpt')
mcc_type = 2;         % 0 = uncorrected; 1 = t-max; 2 = cluster; 3 = TFCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scales = scales(2:end);  % avoid frequencies filtered out by lowpass filter
num_scales = length(scales);

%% MSE

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(MSE1, MSE2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(MSE1, MSE2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
[mask_clusters, summary_tbl] = pull_clusters(mask, tvals, scales, chanlocs, ...
    'nonlinear', grp_type, {size(MSE1,3) size(MSE2,3)},  [], [], [], 'g');
plot_results('nonlinear', 'scalp', scales, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
title("MSE")
set(findall(gcf, 'type', 'axes'), 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_main.fig', coarsing)));
print(gcf, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_main.png', coarsing)), '-dpng', '-r300');
if ~isempty(mask_clusters)
    writetable(summary_tbl, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_summary.csv', coarsing)));
    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, scales, chanlocs, 'MSE', ...
        'DataType', 'scalp', 'Domain', 'nonlinear');
    for i = 1:numel(hs.curve)
        xlim(findobj(hs.curve{i}, 'Type', 'axes'), [2 scales(end)]);
        saveas(hs.topo{i}, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_cluster-%g_topo.fig', coarsing, i)));
        print(hs.topo{i}, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_cluster-%g_topo.png', coarsing, i)),'-dpng','-r300');
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_cluster-%g_curve.fig', coarsing, i)));
        print(hs.curve{i}, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_cluster-%g_curve.png', coarsing, i)),'-dpng','-r300');
    end
    % close([hs.topo{:} hs.curve{:}]) % to close figures
end


%% mMSE

% scales_bounds(1) = []; 

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(mMSE1, mMSE2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(mMSE1, mMSE2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
if any(mask, 'all')
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, scales, chanlocs, ...
        'nonlinear', grp_type, {size(mMSE1,3) size(mMSE2,3)}, 2, [], [], 'g');
    plot_results('nonlinear', 'scalp', scales, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
    title("mMSE"); set(findall(gcf, 'type', 'axes'), 'FontSize', 16, 'FontWeight', 'bold'); 
    ax = findobj(hs.curve{i}, 'Type', 'axes');
    xlim([1 scales(end)]);
    tick_idx = 1:2:scales(end);
    xticks(tick_idx);
    fmt_labels = cellfun(@(s) format_scale_label(s), scales_bounds(tick_idx), 'UniformOutput', false);
    xticklabels(fmt_labels);
    xtickangle(45);
    set(gca, 'FontSize', 14, 'LineWidth', 1.2);
    saveas(gcf, fullfile(outputs_path, sprintf('mMSE_%s_perm_main_uncorrected.fig', coarsing)));
    print(gcf, fullfile(outputs_path, sprintf('mMSE_%s_perm_main_uncorrected.png', coarsing)), '-dpng', '-r300');

    writetable(summary_tbl, fullfile(outputs_path, sprintf('mMSE_%s_perm_summary_uncorrected.csv', coarsing)));
    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, scales, chanlocs, 'mMSE', ...
        'DataType', 'scalp', 'Domain', 'nonlinear');
    for i = 1:numel(hs.curve)
        saveas(hs.topo{i}, fullfile(outputs_path, sprintf('mMSE_%s_perm_cluster-%g_topo_uncorrected.fig', coarsing, i)));
        print(hs.topo{i}, fullfile(outputs_path, sprintf('mMSE_%s_perm_cluster-%g_topo_uncorrected.png', coarsing, i)),'-dpng','-r300');
        ax = findobj(hs.curve{i}, 'Type', 'axes');
        xlim(ax, [1 scales(end)]);
        tick_idx = 1:2:scales(end);
        xticks(ax, tick_idx);
        fmt_labels = cellfun(@(s) format_scale_label(s), scales_bounds(tick_idx), 'UniformOutput', false);
        xticklabels(ax, fmt_labels);
        xtickangle(ax, 45);
        set(gca, 'FontSize', 14, 'LineWidth', 1.2);
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('mMSE_%s_perm_cluster-%g_curve_uncorrected.fig', coarsing, i)));
        print(hs.curve{i}, fullfile(outputs_path, sprintf('mMSE_%s_perm_cluster-%g_curve_uncorrected.png', coarsing, i)),'-dpng','-r300');
    end
    % close([hs.topo{:} hs.curve{:}]) % to close figures
end

%% MFE

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(MFE1, MFE2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(MFE1, MFE2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
if any(mask, 'all')
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, scales, chanlocs, ...
        'nonlinear', grp_type, {size(MFE1,3) size(MFE2,3)}, [], [], [], 'g');
    plot_results('nonlinear', 'scalp', scales, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
    title("MFE"); %set(get(gca,'Title'), 'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_main.fig', coarsing)));
    print(gcf, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_main.png', coarsing)), '-dpng', '-r300');

    writetable(summary_tbl, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_summary.csv', coarsing)));
    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, scales, chanlocs, 'MFE', ...
        'DataType', 'scalp', 'Domain', 'nonlinear');
    for i = 1:numel(hs.curve)
        xlim(findobj(hs.curve{i}, 'Type', 'axes'), [2 scales(end)]);
        saveas(hs.topo{i}, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_cluster-%g_topo.fig', coarsing, i)));
        print(hs.topo{i}, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_cluster-%g_topo.png', coarsing, i)),'-dpng','-r300');
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_cluster-%g_curve.fig', coarsing, i)));
        print(hs.curve{i}, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_cluster-%g_curve.png', coarsing, i)),'-dpng','-r300');
    end
    % close([hs.topo{:} hs.curve{:}]) % to close figures
end

%% RCMFE

if contains(coarsing, {'mean', 'median', 'trimmed mean'})
    scales = 1:scales(end);
    num_scales = length(scales);
    xlims_curve_plot = [1 scales(end)];
else
    RCMFE1(:,1,:) = [];
    RCMFE2(:,1,:) = [];
    xlims_curve_plot = [2 scales(end)];
end

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(RCMFE1, RCMFE2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(RCMFE1, RCMFE2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
if any(mask, 'all')
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, scales, chanlocs, ...
        'nonlinear', grp_type, {size(RCMFE1,3) size(RCMFE2,3)}, [], [], [], 'g');
    plot_results('nonlinear', 'scalp', scales, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
    title("RCMFE"); set(findall(gcf, 'type', 'axes'), 'FontSize', 16, 'FontWeight', 'bold');
    saveas(gcf, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_main.fig', coarsing)));
    print(gcf, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_main.png', coarsing)), '-dpng', '-r300');

    writetable(summary_tbl, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_summary.csv', coarsing)));
    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, scales, chanlocs, 'RCMFE', ...
        'DataType', 'scalp', 'Domain', 'nonlinear');
    for i = 1:numel(hs.curve)
        xlim(findobj(hs.curve{i}, 'Type', 'axes'), xlims_curve_plot);
        saveas(hs.topo{i}, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_cluster-%g_topo.fig', coarsing, i)));
        print(hs.topo{i}, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_cluster-%g_topo.png', coarsing, i)),'-dpng','-r300');
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_cluster-%g_curve.fig', coarsing, i)));
        print(hs.curve{i}, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_cluster-%g_curve.png', coarsing, i)),'-dpng','-r300');
    end
    % close([hs.topo{:} hs.curve{:}]) % to close figures
end


%% PSD (raw)

% PSD1 = log10(PSD1);
% PSD2 = log10(PSD2);
% fprintf('PSD1 negative/zero values: %d\n', sum(PSD1(:) <= 0));
% fprintf('PSD2 negative/zero values: %d\n', sum(PSD2(:) <= 0));

% figure('color','w'); hold on 
% % plot(freqs, mean(mean(PSD1,3),1),'LineWidth',2)
% % plot(freqs, mean(mean(PSD2,3),1),'LineWidth',2)
% % title("Linear")
% plot(freqs, mean(mean(10*log10(PSD1),3),1),'LineWidth',2)
% plot(freqs, mean(mean(10*log10(PSD2),3),1),'LineWidth',2)
% title("dB-norm")
% legend('eyes-closed','eyes-open')

y_unit = 'dB';   % '\muV^2/Hz'  or 'dB'

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(10*log10(PSD1), 10*log10(PSD2), nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(10*log10(PSD1), 10*log10(PSD2), nPerm, ct, grp_type);
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(PSD1, PSD2, nPerm, ct, grp_type);
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(PSD1, PSD2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
if any(mask, 'all')
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, freqs, chanlocs, ...
        'frequency', grp_type, {size(PSD1,3) size(PSD2,3)}, 1, [], [], 'g');
    writetable(summary_tbl, fullfile(outputs_path, 'PSD_raw_summary.csv'));

    plot_results('frequency', 'scalp', freqs, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
    title('Eyes closed vs. Eyes open - Classic PSD'); 
    set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf, fullfile(outputs_path, 'PSD_raw_main.fig'));
    print(gcf, fullfile(outputs_path, 'PSD_raw_main.png'), '-dpng', '-r300');

    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, freqs, chanlocs, ...
        sprintf('Power (%s)', y_unit), 'DataType', 'scalp','Domain','Frequency');
    for i = 1:numel(hs.curve)
        saveas(hs.topo{i},  fullfile(outputs_path, sprintf('PSD_raw_cluster-%g_topo.fig', i)));
        print(hs.topo{i},   fullfile(outputs_path, sprintf('PSD_raw_cluster-%g_topo.png', i)), '-dpng', '-r300');
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('PSD_raw_cluster-%g_curve.fig', i)));
        print(hs.curve{i},  fullfile(outputs_path, sprintf('PSD_raw_cluster-%g_curve.png', i)), '-dpng', '-r300');
    end
    % close([hs.topo{:} hs.curve{:}])
end

%% PSD (Aperiodic-corrected)

% figure('color','w'); hold on 
% plot(freqs, mean(mean(PSD_corr1,3),1),'LineWidth',2)
% plot(freqs, mean(mean(PSD_corr2,3),1),'LineWidth',2)
% title("Linear")
% % plot(freqs, mean(mean(10*log10(PSD_corr1),3),1),'LineWidth',2)
% % plot(freqs, mean(mean(10*log10(PSD_corr2),3),1),'LineWidth',2)
% % title("dB-norm")
% legend('eyes-closed','eyes-open')

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(10*log10(PSD_corr1), 10*log10(PSD_corr2), nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(10*log10(PSD_corr1), 10*log10(PSD_corr2), nPerm, ct, grp_type);
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(PSD_corr1, PSD_corr2, nPerm, ct, grp_type);
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(PSD_corr1, PSD_corr2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
if any(mask, 'all')
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, freqs, chanlocs, ...
        'frequency', grp_type, {size(PSD_corr1,3) size(PSD_corr2,3)}, 0.5, [], [], 'g');
    writetable(summary_tbl, fullfile(outputs_path, 'PSD_corrected_summary.csv'));

    plot_results('frequency', 'scalp', freqs, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
    title('Eyes closed vs. Eyes open - Aperiodic-corrected PSD')
    set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf, fullfile(outputs_path, 'PSD_corrected_main.fig'));
    print(gcf, fullfile(outputs_path, 'PSD_corrected_main.png'), '-dpng', '-r300');

    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, freqs, chanlocs, ...
        sprintf('Power (%s)', y_unit), 'DataType', 'scalp', 'Domain', 'frequency');
    for i = 1:numel(hs.curve)
        saveas(hs.topo{i},  fullfile(outputs_path, sprintf('PSD_corrected_cluster-%g_topo.fig', i)));
        print(hs.topo{i},   fullfile(outputs_path, sprintf('PSD_corrected_cluster-%g_topo.png', i)), '-dpng', '-r300');
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('PSD_corrected_cluster-%g_curve.fig', i)));
        print(hs.curve{i},  fullfile(outputs_path, sprintf('PSD_corrected_cluster-%g_curve.png', i)), '-dpng', '-r300');
    end
    close([hs.topo{:} hs.curve{:}])
end


%% try with eeglab/fieldtrip

% % --- Build FieldTrip freq structures (one per subject per condition) ---
% for s = 1:size(MSE1, 3)
%     freq1{s}.label     = {chanlocs.labels}';   % [64 x 1] cell
%     freq1{s}.freq      = 1:num_scales;          % scales as pseudo-frequencies
%     freq1{s}.powspctrm = MSE1(:,:,s);           % [64 x 50]
%     freq1{s}.dimord    = 'chan_freq';
%     freq1{s}.time      = 1;                     % dummy (required by some FT versions)
% 
%     freq2{s}.label     = {chanlocs.labels}';
%     freq2{s}.freq      = 1:num_scales;
%     freq2{s}.powspctrm = MSE2(:,:,s);
%     freq2{s}.dimord    = 'chan_freq';
%     freq2{s}.time      = 1;
% end
% 
% % --- Prepare neighbours from chanlocs ---
% cfg_neigh            = [];
% cfg_neigh.method     = 'triangulation';   % or 'distance'
% cfg_neigh.elec       = chanlocs2ft(chanlocs);  % convert EEGLAB -> FT format
% neighbours           = ft_prepare_neighbours(cfg_neigh);
% 
% % --- Run cluster permutation stats ---
% cfg                       = [];
% cfg.method                = 'montecarlo';
% cfg.statistic             = 'indepsamplesT';   % or 'depsamplesT' if within-subject
% cfg.correctm              = 'cluster';         % or 'tfce'
% cfg.clusteralpha          = 0.05;
% cfg.clusterstatistic      = 'maxsum';
% cfg.tail                  = 0;                 % two-tailed
% cfg.alpha                 = alpha;
% cfg.numrandomization      = nPerm;
% cfg.neighbours            = neighbours;
% cfg.avgoverfreq           = 'no';             % keep all scales
% cfg.design                = [ones(1,size(MSE1,3)) 2*ones(1,size(MSE2,3))];
% cfg.ivar                  = 1;
% 
% stat = ft_freqstatistics(cfg, freq1{:}, freq2{:});
% 
% % stat.stat         → [64 x 50] t-values
% % stat.mask         → [64 x 50] logical significant mask
% % stat.posclusters  → struct with cluster-level p-values
% % stat.negclusters  → same for negative clusters


%% Local helper


function lbl = format_scale_label(s)
    nums = regexp(s, '[\d.]+', 'match');
    if numel(nums) == 2
        lbl = sprintf('%.2f-%.2f Hz', str2double(nums{2}), str2double(nums{1}));  % flipped
    else
        lbl = sprintf('HP>%.2f Hz', str2double(nums{1}));
    end
end

% function elec = chanlocs2ft(chanlocs)
%     elec.label  = {chanlocs.labels}';
%     elec.elecpos = [[chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]'];
%     elec.chanpos = elec.elecpos;
%     elec.unit    = 'mm';
% end