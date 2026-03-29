% ascent_group_analysis
clear; close all; clc

data_path = 'C:\Users\ccann\Documents\biosemi_data';
pluginPath = fileparts(which('eegplugin_ascent.m'));
addpath(genpath(pluginPath))
cd(pluginPath)
eeglab; close

%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
num_scales = 50;
num_chan = 64;
coarsing = 'mean';
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
MSE = nan(num_chan, num_scales, num_files);
mMSE = nan(num_chan, num_scales, num_files);
MFE = nan(num_chan, num_scales, num_files);
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

    EEG = pop_biosig(filenames{iFile});
    EEG = pop_resample(EEG, 256);
    EEG = pop_select( EEG, 'chantype',{'EEG'});
    EEG = pop_chanedit(EEG, {'lookup','standard_1005.elc'});
    EEG = reorder_channels(EEG);
    EEG = ref_infinity(EEG);
    EEG = pop_eegfiltnew(EEG, 'locutoff',1);
    % EEG = pop_eegfiltnew(EEG, 'locutoff',59,'hicutoff',61,'revfilt',1);
    EEG = pop_eegfiltnew(EEG,'hicutoff',50);
    oriEEG = EEG;
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.75, ...
        'LineNoiseCriterion',5,'Highpass','off','BurstCriterion',40, ...
        'WindowCriterion','off','BurstRejection','off','Distance','Euclidian', ...
        'WindowCriterionTolerances','off');
    % vis_artifacts(EEG,oriEEG);
    EEG = pop_interp(EEG, oriEEG.chanlocs, 'spherical');
    dataRank = sum(eig(cov(double(EEG.data'))) > 1e-7);
    EEG = pop_runica(EEG, 'icatype', 'picard', 'mode', 'standard', ...
        'maxiter', 500, 'pca', dataRank);
    EEG = pop_iclabel(EEG, 'default');
    EEG = pop_icflag(EEG, [NaN NaN;0.9 1;0.9 1;0.95 1;0.95 1;0.9 1;NaN NaN]);
    % pop_selectcomps(EEG, 1:24);
    EEG = pop_subcomp(EEG, [], 0);
    % pop_eegplot(EEG,1,1,1);
    
    % Load if already preprocessed
    % EEG = pop_loadset('filepath', fullfile(data_path, 'eyes_closed'), 'filename', sprintf('%s.set', filenames{iFile}(1:end-4)));

    % Compute complexity measures
    EEG = ascent_compute(EEG, 'measure', 'SampEn', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'FuzzEn', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'HigFracDim', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'Aperiodic', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'MSE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'mMSE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'ExSEnt', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'MFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'RCMFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    % EEG = ascent_compute(EEG, 'measure', 'RCmvMFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);

    % ascent_plot(EEG.ascent.MSE.data, EEG.chanlocs, 'mse', EEG.ascent.MSE.scales);

    % Save for running entropy on ICA time series later
    EEG = pop_saveset(EEG, 'filepath', fullfile(data_path, 'eyes_closed'), 'filename', sprintf('%s.set', filenames{iFile}(1:end-4)));

    SampEn(:,iFile) = EEG.ascent.SampEn.data;
    ExSEnt1(:,iFile) = EEG.ascent.ExSEnt.data.HD;
    ExSEnt2(:,iFile) = EEG.ascent.ExSEnt.data.HA;
    ExSEnt3(:,iFile) = EEG.ascent.ExSEnt.data.HDA;
    FuzzEn(:,iFile) = EEG.ascent.FuzzEn.data;
    FracDim(:,iFile) = EEG.ascent.HigFracDim.data;
    Exponent(:,iFile) = EEG.ascent.Aperiodic.data.exponent;
    Offset(:,iFile) = EEG.ascent.Aperiodic.data.offset;
    MSE(:,:,iFile) = EEG.ascent.MSE.data;
    mMSE(:,:,iFile) = EEG.ascent.mMSE.data;
    MFE(:,:,iFile) = EEG.ascent.MFE.data;
    RCMFE(:,:,iFile) = EEG.ascent.RCMFE.data;
    % RCmvMFE(:,:,iFile) = EEG.ascent.RCmvMFE.data;

    chanlocs = EEG.chanlocs;
    scales = EEG.ascent.MSE.scales;
    scales_bounds = EEG.ascent.mMSE.scales;

    save(fullfile(data_path, sprintf('ascent_outputs_EC_%s_new.mat', coarsing)), ...
        "SampEn", "ExSEnt1", "ExSEnt2", "ExSEnt3", "FuzzEn", "FracDim", "Exponent", "Offset", ...
        "MSE", "MFE", "mMSE", "RCMFE", ...
        'chanlocs', 'scales', 'scales_bounds')

end
disp("Done computing on the whole group for eyes-closed condition!")

% Check there are no subjects with NaNs
nan_subject = any(squeeze(any(isnan(MFE), 1)));  % [50 x 40] logical (scale x subject)
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
MSE = nan(num_chan, num_scales, num_files);
mMSE = nan(num_chan, num_scales, num_files);
MFE = nan(num_chan, num_scales, num_files);
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

    EEG = pop_biosig(filenames{iFile});
    EEG = pop_resample(EEG, 256);
    EEG = pop_select( EEG, 'chantype',{'EEG'});
    EEG = pop_chanedit(EEG, {'lookup','standard_1005.elc'});
    EEG = reorder_channels(EEG);
    EEG = ref_infinity(EEG);
    EEG = pop_eegfiltnew(EEG, 'locutoff',1);
    % EEG = pop_eegfiltnew(EEG, 'locutoff',59,'hicutoff',61,'revfilt',1);
    EEG = pop_eegfiltnew(EEG,'hicutoff',50);
    oriEEG = EEG;
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.75, ...
        'LineNoiseCriterion',5,'Highpass','off','BurstCriterion',40, ...
        'WindowCriterion','off','BurstRejection','off','Distance','Euclidian', ...
        'WindowCriterionTolerances','off');
    % vis_artifacts(EEG,oriEEG);
    EEG = pop_interp(EEG, oriEEG.chanlocs, 'spherical');
    dataRank = sum(eig(cov(double(EEG.data'))) > 1e-7);
    EEG = pop_runica(EEG, 'icatype', 'picard', 'mode', 'standard','pca', dataRank);
    EEG = pop_iclabel(EEG, 'default');
    EEG = pop_icflag(EEG, [NaN NaN;0.9 1;0.9 1;0.95 1;0.95 1;0.9 1;NaN NaN]);
    % pop_selectcomps(EEG, 1:24);
    EEG = pop_subcomp(EEG, [], 0);
    % pop_eegplot(EEG,1,1,1);

    % % Load if already processed
    % EEG = pop_loadset('filepath', fullfile(data_path, 'eyes_open'), 'filename', sprintf('%s.set', filenames{iFile}(1:end-4)));

    % Compute complexity measures
    EEG = ascent_compute(EEG, 'measure', 'SampEn', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'FuzzEn', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'HigFracDim', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'Aperiodic', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'MSE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'mMSE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'ExSEnt', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'MFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'RCMFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    % EEG = ascent_compute(EEG, 'measure', 'RCmvMFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);

    % Save for running entropy on ICA time series later
    EEG = pop_saveset(EEG, 'filepath', fullfile(data_path, 'eyes_open'), 'filename', sprintf('%s.set', filenames{iFile}(1:end-4)));

    SampEn(:,iFile) = EEG.ascent.SampEn.data;
    ExSEnt1(:,iFile) = EEG.ascent.ExSEnt.data.HD;
    ExSEnt2(:,iFile) = EEG.ascent.ExSEnt.data.HA;
    ExSEnt3(:,iFile) = EEG.ascent.ExSEnt.data.HDA;
    FuzzEn(:,iFile) = EEG.ascent.FuzzEn.data;
    FracDim(:,iFile) = EEG.ascent.HigFracDim.data;
    Exponent(:,iFile) = EEG.ascent.Aperiodic.data.exponent;
    Offset(:,iFile) = EEG.ascent.Aperiodic.data.offset;
    MSE(:,:,iFile) = EEG.ascent.MSE.data;
    mMSE(:,:,iFile) = EEG.ascent.mMSE.data;
    MFE(:,:,iFile) = EEG.ascent.MFE.data;
    RCMFE(:,:,iFile) = EEG.ascent.RCMFE.data;
    % RCmvMFE(:,:,iFile) = EEG.ascent.RCmvMFE.data;

    chanlocs = EEG.chanlocs;
    scales = EEG.ascent.MSE.scales;
    scales_bounds = EEG.ascent.mMSE.scales;

    save(fullfile(data_path, sprintf('ascent_outputs_EO_%s_new.mat', coarsing)), ...
        "SampEn", "ExSEnt1", "ExSEnt2", "ExSEnt3", "FuzzEn", "FracDim", "Exponent", "Offset", ...
        "MSE", "MFE", "mMSE", "RCMFE", ...
        'chanlocs', 'scales', 'scales_bounds')

end
disp("Done computing on the whole group for eyes-open condition!")

% Check there are no subjects with NaNs
nan_subject = any(squeeze(any(isnan(MFE), 1)));  % [50 x 40] logical (scale x subject)
if any(nan_subject)
    error("some subjects have NaNs! Check data!")
    disp(nan_subject)
end

gong

%% Load, separate the data by condition, & reorganize electrodes

% coarsing = 'mean';
% 
% % Load chanlocs only
% load(fullfile(data_path, sprintf('ascent_outputs_EC_%s_new.mat', coarsing)), 'chanlocs')
% cd(data_path)
% 
% % Front-to-back target order
% desired_order = {
%     'Fp1','FPz','FP2', ...
%     'AF7','AF3','AFz','AF4','AF8', ...
%     'F7','F5','F3','F1','Fz','F2','F4','F6','F8', ...
%     'FT7','FT8', ...
%     'FC5','FC3','FC1','FCz','FC2','FC4','FC6', ...
%     'T7','T8', ...
%     'C5','C3','C1','Cz','C2','C4','C6', ...
%     'TP7','TP8', ...
%     'CP5','CP3','CP1','CPz','CP2','CP4','CP6', ...
%     'P9','P7','P5','P3','P1','Pz','P2','P4','P6','P8','P10', ...
%     'PO7','PO3','POz','PO4','PO8', ...
%     'O1','Oz','Iz','O2'
%     };
% 
% % Map current labels -> desired positions, then sort
% current_labels = {chanlocs.labels};
% [found, loc_in_desired] = ismember(lower(current_labels), lower(desired_order));
% if ~all(found)
%     missing = current_labels(~found);
%     error('Missing in desired_order: %s', strjoin(missing, ', '));
% end
% [~, order_idx] = sort(loc_in_desired, 'ascend');
% 
% % Reorder chanlocs
% chanlocs = chanlocs(order_idx);
% % save(fullfile(data_path, 'chanlocs_reorganized.mat'), 'chanlocs')
% disp({chanlocs.labels}')
% 
% % Load EC metrics and reorder
% load(fullfile(data_path, sprintf('ascent_outputs_EC_%s_new.mat', coarsing)), ...
%     'SampEn','FuzzEn','ExSEnt1','ExSEnt2','ExSEnt3','FracDim','Exponent', ...
%     'Offset','MSE','MFE','mMSE','RCMFE','scales', 'scales_bounds')
% SampEn1  = SampEn(order_idx, :);
% FuzzEn1  = FuzzEn(order_idx, :);
% ExSEnt1_1  = ExSEnt1(order_idx, :);
% ExSEnt1_2  = ExSEnt2(order_idx, :);
% ExSEnt1_3  = ExSEnt3(order_idx, :);
% FracDim1 = FracDim(order_idx, :);
% Exponent1 = Exponent(order_idx,:);
% Offset1 = Offset(order_idx,:);
% MSE1     = MSE(order_idx, :, :);
% MFE1     = MFE(order_idx, :, :);
% mMSE1    = mMSE(order_idx, :, :);
% RCMFE1   = RCMFE(order_idx, :, :);
% % RCmvMFE1   = RCmvMFE(order_idx, :, :);
% 
% % Load EO metrics and reorder
% load(fullfile(data_path, sprintf('ascent_outputs_EO_%s_new.mat', coarsing)), ...
%     'SampEn','FuzzEn','ExSEnt1','ExSEnt2','ExSEnt3','FracDim','Exponent', ...
%     'Offset','MSE','MFE','mMSE','RCMFE','scales', 'scales_bounds')
% SampEn2  = SampEn(order_idx, :);
% FuzzEn2  = FuzzEn(order_idx, :);
% ExSEnt2_1  = ExSEnt1(order_idx, :);
% ExSEnt2_2  = ExSEnt2(order_idx, :);
% ExSEnt2_3  = ExSEnt3(order_idx, :);
% FracDim2 = FracDim(order_idx, :);
% Exponent2 = Exponent(order_idx,:);
% Offset2 = Offset(order_idx,:);
% MSE2     = MSE(order_idx, :, :);
% MFE2     = MFE(order_idx, :, :);
% mMSE2    = mMSE(order_idx, :, :);
% RCMFE2   = RCMFE(order_idx, :, :);
% % RCmvMFE2   = RCmvMFE(order_idx, :, :);


%% Stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
nPerm = 5000;          % number of permutations for H0
alpha = 0.01;           % NaN to show maps of p-values
ct = 'mean';            % central tendency method ('mean' for normal t-test, 'trimmed mean' for Yuen t-test)
grp_type = 'dpt';       % groupe is dependent ('dpt') or independent ('idpt')
mcc_type = 2;   % 0 = uncorrected; 1 = t-max; 2 = cluster; 3 = TFCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coarsing = 'mean';

% load(fullfile(data_path,'chanlocs_reorganized.mat'), 'chanlocs')
outputs_path = fullfile(pluginPath, 'figures_new', 'group_results'); mkdir(outputs_path)
load colormap_bwr.mat; dmap(1,:) = [0.9 0.9 0.9]; % set NaNs to grey


%%%%%%%%%%%%%%%%%%%%% UNISCALES %%%%%%%%%%%%%%%%%%%%%

figure('Color','w');

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
% if any(mask)
%     topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
%         'verbose','off','whitebk','on');
%     c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
%     % title(sprintf("ExSEnt (coarse: %s)", coarsing));
%     title("ExSEnt (amp + dur)")
% else
%     close(gcf)
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 0, 0.05, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
    'verbose','off','whitebk','on');
c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
% title(sprintf("ExSEnt (coarse: %s)", coarsing));
title("ExSEnt (amp + dur; uncorrected)")
% end

% FracDim
% nexttile
subplot(2,4,6)
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(FracDim1, FracDim2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(FracDim1, FracDim2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
% [mask, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,alpha,'pdep','yes');
if any(mask)
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("FracDim (coarse: %s)", coarsing));
    title("FracDim")
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
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, ...
        'verbose','off','whitebk','on');
    c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',12)
    % title(sprintf("FracDim (coarse: %s)", coarsing));
    title("Aperiodic Offset")
else
    close(gcf)
end

set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');

saveas(gcf, fullfile(outputs_path, 'fig3.fig'));
print(gcf, fullfile(outputs_path, 'fig3.png'), '-dpng', '-r300');

% % set(gcf, 'Renderer', 'painters');   % or 'opengl'
% exportgraphics(gcf, ...
%     fullfile(outputs_path,'single-scales.png'), ...
%     'Resolution', 300);
% fig = gcf;
% exportgraphics(fig, fullfile(outputs_path,'single-scales.png'), 'Resolution',300);
%
% fig = gcf;
% set(fig,'Color','w');
% opengl software
% print(fig, fullfile(outputs_path,'single-scales.png'), '-dpng', '-r300');

%% %%%%%%%%%%%%%%%%%%%% MULTISCALES %%%%%%%%%%%%%%%%%%%%%

coarsing = 'mean';
mcc_type = 3;   % 0 = uncorrected; 1 = t-max; 2 = cluster; 3 = TFCE

% MSE
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(MSE1, MSE2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(MSE1, MSE2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
[mask_clusters, summary_tbl] = pull_clusters(mask, tvals, scales, chanlocs, ...
    'nonlinear', grp_type, {size(MSE1,3) size(MSE2,3)},  2, 1, [], 'g');
plot_results('nonlinear', 'scalp', scales, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
title("MSE")
set(findall(gcf, 'type', 'axes'), 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_main.fig', coarsing)));
print(gcf, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_main.png', coarsing)), '-dpng', '-r300');
if ~isempty(mask_clusters)
    writetable(summary_tbl, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_summary.csv', coarsing)));
    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, scales, chanlocs, 'MSE', 'DataType', 'scalp');
    for i = 1:numel(hs.curve)
        saveas(hs.topo{i}, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_cluster-%g_topo.fig', coarsing, i)));
        print(hs.topo{i}, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_cluster-%g_topo.png', coarsing, i)),'-dpng','-r300');
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_cluster-%g_curve.fig', coarsing, i)));
        print(hs.curve{i}, fullfile(outputs_path, sprintf('MSE_%s_perm_tfce_cluster-%g_curve.png', coarsing, i)),'-dpng','-r300');
    end
    close([hs.topo{:} hs.curve{:}]) % to close figures
end

%% mMSE

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(mMSE1, mMSE2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(mMSE1, mMSE2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
[mask_clusters, summary_tbl] = pull_clusters(mask, tvals, scales, chanlocs, ...
    'nonlinear', grp_type, {size(mMSE1,3) size(mMSE2,3)}, 2, 1, [], 'g');
plot_results('nonlinear', 'scalp', scales, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
title("mMSE")
set(findall(gcf, 'type', 'axes'), 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_main.fig', coarsing)));
print(gcf, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_main.png', coarsing)), '-dpng', '-r300');
if ~isempty(mask_clusters)
    writetable(summary_tbl, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_summary.csv', coarsing)));
    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, scales, chanlocs, 'mMSE', 'DataType', 'scalp');
    xticklabels(scales_bounds);
    for i = 1:numel(hs.curve)
        saveas(hs.topo{i}, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_cluster-%g_topo.fig', coarsing, i)));
        print(hs.topo{i}, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_cluster-%g_topo.png', coarsing, i)),'-dpng','-r300');
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_cluster-%g_curve.fig', coarsing, i)));
        print(hs.curve{i}, fullfile(outputs_path, sprintf('mMSE_%s_perm_tfce_cluster-%g_curve.png', coarsing, i)),'-dpng','-r300');
    end
    close([hs.topo{:} hs.curve{:}]) % to close figures
end

%% MFE

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(MFE1, MFE2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(MFE1, MFE2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
[mask_clusters, summary_tbl] = pull_clusters(mask, tvals, scales, chanlocs, ...
    'nonlinear', grp_type, {size(MFE1,3) size(MFE2,3)}, [], [], [], 'g');
plot_results('nonlinear', 'scalp', scales, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
title("MFE")
set(findall(gcf, 'type', 'axes'), 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_main.fig', coarsing)));
print(gcf, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_main.png', coarsing)), '-dpng', '-r300');
if ~isempty(mask_clusters)
    writetable(summary_tbl, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_summary.csv', coarsing)));
    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, scales, chanlocs, 'MFE', 'DataType', 'scalp', 'LineNoiseHz', []);
    % xticklabels(scales_bounds);
    for i = 1:numel(hs.curve)
        saveas(hs.topo{i}, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_cluster-%g_topo.fig', coarsing, i)));
        print(hs.topo{i}, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_cluster-%g_topo.png', coarsing, i)),'-dpng','-r300');
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_cluster-%g_curve.fig', coarsing, i)));
        print(hs.curve{i}, fullfile(outputs_path, sprintf('MFE_%s_perm_tfce_cluster-%g_curve.png', coarsing, i)),'-dpng','-r300');
    end
    close([hs.topo{:} hs.curve{:}]) % to close figures
end

%% RCMFE

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(RCMFE1, RCMFE2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(RCMFE1, RCMFE2, nPerm, ct, grp_type);
mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
[mask_clusters, summary_tbl] = pull_clusters(mask, tvals, scales, chanlocs, ...
    'nonlinear', grp_type, {size(RCMFE1,3) size(RCMFE2,3)}, [], [], [], 'g');
plot_results('nonlinear', 'scalp', scales, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
title("RCMFE")
set(findall(gcf, 'type', 'axes'), 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_main.fig', coarsing)));
print(gcf, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_main.png', coarsing)), '-dpng', '-r300');
if ~isempty(mask_clusters)
    writetable(summary_tbl, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_summary.csv', coarsing)));
    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, scales, chanlocs, 'RCMFE', 'DataType', 'scalp', 'LineNoiseHz', []);
    xticklabels(scales_bounds);
    for i = 1:numel(hs.curve)
        saveas(hs.topo{i}, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_cluster-%g_topo.fig', coarsing, i)));
        print(hs.topo{i}, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_cluster-%g_topo.png', coarsing, i)),'-dpng','-r300');
        saveas(hs.curve{i}, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_cluster-%g_curve.fig', coarsing, i)));
        print(hs.curve{i}, fullfile(outputs_path, sprintf('RCMFE_%s_perm_tfce_cluster-%g_curve.png', coarsing, i)),'-dpng','-r300');
    end
    close([hs.topo{:} hs.curve{:}]) % to close figures
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
