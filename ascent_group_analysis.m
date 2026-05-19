% ascent_group_analysis
clear; close all; clc

data_path = 'C:\Users\ccann\Documents\biosemi_data';
pluginPath = fileparts(which('eegplugin_ascent.m'));
addpath(genpath(pluginPath))
cd(pluginPath)
eeglab; close

%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
num_scales = 30;
num_chan = 64;
coarsing = 'sd';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute on whole group - Eyes closed (EC) condition

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

% coarsing = 'sd';

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

% scales = scales(1:end);  % avoid frequencies filtered out by lowpass filter
% scales_bounds = scales_bounds(3:end);  % avoid frequencies filtered out by lowpass filter
% num_scales = numel(scales);  % update

%% Stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
nPerm = 1000;          % number of permutations for H0
alpha = 0.05;           % NaN to show maps of p-values
ct = 'mean';            % central tendency method ('mean' for normal t-test, 'trimmed mean' for Yuen t-test)
grp_type = 'dpt';       % groupe is dependent ('dpt') or independent ('idpt')
mcc_type = 2;   % 0 = uncorrected; 1 = t-max; 2 = cluster; 3 = TFCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load(fullfile(data_path,'chanlocs_reorganized.mat'), 'chanlocs')
outputs_path = fullfile(pluginPath, 'figures_new', 'group_results'); mkdir(outputs_path)
% load colormap_bwr.mat; dmap(1,:) = [0.9 0.9 0.9]; % set NaNs to grey
rdbu_cmap    = interp1([1;128;256],[0.698 0.094 0.168;1 1 1;0.129 0.400 0.675],(1:256)','linear');
dmap    = flipud(max(0,min(1,rdbu_cmap)));

%%%%%%%%%%%%%%%%%%%%% UNISCALES %%%%%%%%%%%%%%%%%%%%%

figure('Color','w','ToolBar','none','MenuBar','none');

% SampEn
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
    % tvals(~mask) = 0;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("SampEn (coarse: %s)", coarsing));
    title("SampEn")
else
    close(gcf)
end

% FuzzEn
% nexttile
subplot(2,4,2)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(FuzzEn1, FuzzEn2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(FuzzEn1, FuzzEn2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    % tvals(~mask) = 0;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("FuzzEn (coarse: %s)", coarsing));
    title("FuzzEn")
else
    close(gcf)
end

% ExSEnt1
% nexttile
subplot(2,4,3)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(ExSEnt1_1, ExSEnt2_1, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(ExSEnt1_1, ExSEnt2_1, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    % tvals(~mask) = 0;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("ExSEnt (coarse: %s)", coarsing));
    title("ExSEnt (duration)")
else
    close(gcf)
end

% ExSEnt2
% nexttile
subplot(2,4,4)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(ExSEnt1_2, ExSEnt2_2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(ExSEnt1_2, ExSEnt2_2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    % tvals(~mask) = NaN;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("ExSEnt (coarse: %s)", coarsing));
    title("ExSEnt (amplitude)")
else
    close(gcf)
end

% ExSEnt3
% nexttile
subplot(2,4,5)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(ExSEnt1_3, ExSEnt2_3, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(ExSEnt1_3, ExSEnt2_3, nPerm, ct, grp_type);
% mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% % [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    % tvals(~mask) = NaN;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("ExSEnt (coarse: %s)", coarsing));
    title("ExSEnt (amp + dur)")
else
    close(gcf)
% mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 0, 0.05, chanlocs);
% % [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
% topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
%     'verbose','off','whitebk','on');
% c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
% % title(sprintf("ExSEnt (coarse: %s)", coarsing));
% title("ExSEnt (amp + dur; uncorrected)")
end

% FracDim
% nexttile
subplot(2,4,6)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(FracDim1, FracDim2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(FracDim1, FracDim2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    tvals(mask) = 0;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("FracDim (coarse: %s)", coarsing));
    title("HigFracDim")
else
    close(gcf)
end

% Aperiodic Exponent
% nexttile
subplot(2,4,7)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(Exponent1, Exponent2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(Exponent1, Exponent2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    % tvals(~mask) = 0;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("FracDim (coarse: %s)", coarsing));
    title("Aperiodic Exponent")
else
    close(gcf)
end

% Aperiodic Ofset
% nexttile
subplot(2,4,8)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(Offset1, Offset2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(Offset1, Offset2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    % tvals(~mask) = 0;
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("FracDim (coarse: %s)", coarsing));
    title("Aperiodic Offset")
else
    close(gcf)
end

set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');

saveas(gcf, fullfile(outputs_path, 'fig_uniscales.fig'));
print(gcf, fullfile(outputs_path, 'fig_uniscales.png'), '-dpng', '-r300');


%% %%%%%%%%%%%%%%%%%%%% MULTISCALES %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
coarsing = 'mean';
nPerm = 1000;          % number of permutations for H0
alpha = 0.05;           % NaN to show maps of p-values
ct = 'mean';            % central tendency method ('mean' for normal t-test, 'trimmed mean' for Yuen t-test)
grp_type = 'dpt';       % groupe is dependent ('dpt') or independent ('idpt')
mcc_type = 3;         % 0 = uncorrected; 1 = t-max; 2 = cluster; 3 = TFCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
        saveas(hs.topo{i}, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_cluster-%g_topo.fig', coarsing, i)));
        print(hs.topo{i}, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_cluster-%g_topo.png', coarsing, i)),'-dpng','-r300');
        xlim(findobj(hs.curve{i}, 'Type', 'axes'), [1 num_scales]);
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_cluster-%g_curve.fig', coarsing, i)));
        print(hs.curve{i}, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_cluster-%g_curve.png', coarsing, i)),'-dpng','-r300');
    end
    % close([hs.topo{:} hs.curve{:}]) % to close figures
end

%% Same but with eeglab/fieldtrip

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

%% mMSE

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(mMSE1, mMSE2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(mMSE1, mMSE2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
if any(mask, 'all')
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, scales, chanlocs, ...
        'nonlinear', grp_type, {size(mMSE1,3) size(mMSE2,3)}, [], [], [], 'g');
    plot_results('nonlinear', 'scalp', scales, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
    title("mMSE"); set(findall(gcf, 'type', 'axes'), 'FontSize', 16, 'FontWeight', 'bold');
    ax = findobj(hs.curve{i}, 'Type', 'axes');
    xlim([1 num_scales]);
    tick_idx = 1:2:num_scales;
    xticks(tick_idx);
    fmt_labels = cellfun(@(s) format_scale_label(s), scales_bounds(tick_idx), 'UniformOutput', false);
    xticklabels(fmt_labels);
    xtickangle(45);
    set(gca, 'FontSize', 14, 'LineWidth', 1.2);
    saveas(gcf, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_main.fig', coarsing)));
    print(gcf, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_main.png', coarsing)), '-dpng', '-r300');

    writetable(summary_tbl, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_summary.csv', coarsing)));
    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, scales, chanlocs, 'mMSE', ...
        'DataType', 'scalp', 'Domain', 'nonlinear');
    for i = 1:numel(hs.curve)
        saveas(hs.topo{i}, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_cluster-%g_topo.fig', coarsing, i)));
        print(hs.topo{i}, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_cluster-%g_topo.png', coarsing, i)),'-dpng','-r300');
        ax = findobj(hs.curve{i}, 'Type', 'axes');
        xlim(ax, [1 num_scales]);
        tick_idx = 1:2:num_scales;
        xticks(ax, tick_idx);
        fmt_labels = cellfun(@(s) format_scale_label(s), scales_bounds(tick_idx), 'UniformOutput', false);
        xticklabels(ax, fmt_labels);
        xtickangle(ax, 45);
        set(gca, 'FontSize', 14, 'LineWidth', 1.2);
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_cluster-%g_curve.fig', coarsing, i)));
        print(hs.curve{i}, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_cluster-%g_curve.png', coarsing, i)),'-dpng','-r300');
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
    title("MFE"); set(findall(gcf, 'type', 'axes'), 'FontSize', 16, 'FontWeight', 'bold');
    saveas(gcf, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_main.fig', coarsing)));
    print(gcf, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_main.png', coarsing)), '-dpng', '-r300');

    writetable(summary_tbl, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_summary.csv', coarsing)));
    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, scales, chanlocs, 'MFE', ...
        'DataType', 'scalp', 'Domain', 'nonlinear');
    for i = 1:numel(hs.curve)
        saveas(hs.topo{i}, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_cluster-%g_topo.fig', coarsing, i)));
        print(hs.topo{i}, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_cluster-%g_topo.png', coarsing, i)),'-dpng','-r300');
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_cluster-%g_curve.fig', coarsing, i)));
        print(hs.curve{i}, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_cluster-%g_curve.png', coarsing, i)),'-dpng','-r300');
    end
    close([hs.topo{:} hs.curve{:}]) % to close figures
end

%% RCMFE

scales = 1:scales(end);
% num_scales = length(scales);

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

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(10*log10(PSD1), 10*log10(PSD2), nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(PSD1, PSD2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
if any(mask, 'all')
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, freqs, chanlocs, ...
        'frequency', grp_type, {size(PSD1,3) size(PSD2,3)}, [], [], [], 'g');
    plot_results('frequency', 'scalp', freqs, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
    title('PSD (raw)')
    set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf, fullfile(outputs_path, 'PSD_raw_perm_tfce_main.fig'));
    print(gcf, fullfile(outputs_path, 'PSD_raw_perm_tfce_main.png'), '-dpng', '-r300');

    writetable(summary_tbl, fullfile(outputs_path, 'PSD_raw_perm_tfce_summary.csv'));
    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, freqs, ...
        chanlocs, 'PSD (raw)', 'DataType', 'scalp','Domain','Frequency');
    for i = 1:numel(hs.curve)
        saveas(hs.topo{i},  fullfile(outputs_path, sprintf('PSD_raw_perm_tfce_cluster-%g_topo.fig', i)));
        print(hs.topo{i},   fullfile(outputs_path, sprintf('PSD_raw_perm_tfce_cluster-%g_topo.png', i)), '-dpng', '-r300');
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('PSD_raw_perm_tfce_cluster-%g_curve.fig', i)));
        print(hs.curve{i},  fullfile(outputs_path, sprintf('PSD_raw_perm_tfce_cluster-%g_curve.png', i)), '-dpng', '-r300');
    end
    % close([hs.topo{:} hs.curve{:}])
end

%% PSD (aperiodic-corrected)

% figure('color','w'); hold on 
% plot(freqs, mean(mean(PSD_corr1,3),1),'LineWidth',2)
% plot(freqs, mean(mean(PSD_corr2,3),1),'LineWidth',2)
% title("Linear")
% % plot(freqs, mean(mean(10*log10(PSD_corr1),3),1),'LineWidth',2)
% % plot(freqs, mean(mean(10*log10(PSD_corr2),3),1),'LineWidth',2)
% % title("dB-norm")
% legend('eyes-closed','eyes-open')

[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(10*log10(PSD_corr1), 10*log10(PSD_corr2), nPerm, ct, grp_type);
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(PSD_corr1, PSD_corr2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
if any(mask, 'all')
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, freqs, chanlocs, ...
        'frequency', grp_type, {size(PSD_corr1,3) size(PSD_corr2,3)}, [], [], [], 'g');
    plot_results('frequency', 'scalp', freqs, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
    title('PSD (aperiodic-corrected)')
    set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf, fullfile(outputs_path, 'PSD_corrected_perm_tfce_main.fig'));
    print(gcf, fullfile(outputs_path, 'PSD_corrected_perm_tfce_main.png'), '-dpng', '-r300');

    writetable(summary_tbl, fullfile(outputs_path, 'PSD_corrected_perm_tfce_summary.csv'));
    % hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, freqs, chanlocs, 'PSD (corrected)', 'DataType', 'scalp');
    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, freqs, chanlocs, 'Power', ...
        'DataType', 'scalp', 'Domain', 'frequency');
    for i = 1:numel(hs.curve)
        saveas(hs.topo{i},  fullfile(outputs_path, sprintf('PSD_corrected_perm_tfce_cluster-%g_topo.fig', i)));
        print(hs.topo{i},   fullfile(outputs_path, sprintf('PSD_corrected_perm_tfce_cluster-%g_topo.png', i)), '-dpng', '-r300');
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('PSD_corrected_perm_tfce_cluster-%g_curve.fig', i)));
        print(hs.curve{i},  fullfile(outputs_path, sprintf('PSD_corrected_perm_tfce_cluster-%g_curve.png', i)), '-dpng', '-r300');
    end
    % close([hs.topo{:} hs.curve{:}])
end


%% Local helper

function EEG = reorder_channels(EEG)

% Reorder to have Front-to-back of head order for better imagesc
% visualization
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

current_labels = {EEG.chanlocs.labels};
[found, loc_in_desired] = ismember(lower(current_labels), lower(desired_order));

if ~all(found)
    missing = current_labels(~found);
    error('Missing in desired_order: %s', strjoin(missing, ', '));
end

[~, order_idx] = sort(loc_in_desired, 'ascend');

EEG.chanlocs = EEG.chanlocs(order_idx);
EEG.data = EEG.data(order_idx,:);

end

function lbl = format_scale_label(s)
    nums = regexp(s, '[\d.]+', 'match');
    if numel(nums) == 2
        lbl = sprintf('%.2f-%.2f Hz', str2double(nums{2}), str2double(nums{1}));  % flipped
    else
        lbl = sprintf('HP>%.2f Hz', str2double(nums{1}));
    end
end



function elec = chanlocs2ft(chanlocs)
    elec.label  = {chanlocs.labels}';
    elec.elecpos = [[chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]'];
    elec.chanpos = elec.elecpos;
    elec.unit    = 'mm';
end