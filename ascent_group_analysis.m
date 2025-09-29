% ascent_group_analysis
clear; close all; clc

data_path = '/Users/cedriccannard/Downloads/biosemi_data';
cd(data_path)
addpath(genpath('/Users/cedriccannard/Documents/MATLAB/Ascent'))
eeglab; close

%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
num_scales = 30;
num_chan = 64;
coarsing = 'median';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Eyes closed

cd(fullfile(data_path, 'eyes_closed'))
filenames = {dir('*.bdf').name}';
num_files = length(filenames);

SampEn = nan(num_chan, num_files);
ExSEnt = nan(num_chan, num_files);
FuzzEn = nan(num_chan, num_files);
FracDim = nan(num_chan, num_files);
MSE = nan(num_chan, num_scales, num_files);
mMSE = nan(num_chan, num_scales, num_files);
MFE = nan(num_chan, num_scales, num_files);
RCMFE = nan(num_chan, num_scales, num_files);
for iFile = 1:num_files
    EEG = pop_biosig(filenames{iFile});
    EEG = pop_resample(EEG, 256);
    EEG = pop_select( EEG, 'chantype',{'EEG'});
    EEG = pop_chanedit(EEG, {'lookup','standard_1005.elc'});
    EEG = ref_infinity(EEG);
    
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'usefftfilt',1);
    EEG = pop_eegfiltnew(EEG, 'locutoff',59,'hicutoff',61,'revfilt',1,'usefftfilt',1);

    EEG = ascent_compute(EEG, 'measure', 'SampEn', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'FuzzEn', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'FracDim', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'MSE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'mMSE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'ExSEnt', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'MFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'RCMFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);

    SampEn(:,iFile) = EEG.ascent.SampEn.data;
    ExSEnt(:,iFile) = EEG.ascent.ExSEnt.data.HDA;
    FuzzEn(:,iFile) = EEG.ascent.FuzzEn.data;
    FracDim(:,iFile) = EEG.ascent.FracDim.data;
    MSE(:,:,iFile) = EEG.ascent.MSE.data;
    mMSE(:,:,iFile) = EEG.ascent.mMSE.data;
    MFE(:,:,iFile) = EEG.ascent.MFE.data;
    RCMFE(:,:,iFile) = EEG.ascent.MSE.data;

    chanlocs = EEG.chanlocs;
    scales = EEG.ascent.MSE.scales;
    scales_bounds = EEG.ascent.mMSE.scales;

    save(fullfile(data_path, sprintf('ascent_outputs_EC_%s.mat', coarsing)), ...
        "SampEn", "ExSEnt", "FuzzEn", "FracDim", ...
        "MSE", "MFE", "mMSE", "RCMFE", ...
        'chanlocs', 'scales', 'scales_bounds')

end
disp("Done computing on the whole group for eyes-closed condition!")

%% Eyes opened

cd(fullfile(data_path, 'eyes_open'))
filenames = {dir('*.bdf').name}';
num_files = length(filenames);

SampEn = nan(num_chan, num_files);
ExSEnt = nan(num_chan, num_files);
FuzzEn = nan(num_chan, num_files);
FracDim = nan(num_chan, num_files);
MSE = nan(num_chan, num_scales, num_files);
mMSE = nan(num_chan, num_scales, num_files);
MFE = nan(num_chan, num_scales, num_files);
RCMFE = nan(num_chan, num_scales, num_files);
for iFile = 1:num_files

    fprintf('--------------------------------------------------------\n')
    fprintf('                        FILE %g/%g \n', iFile, num_files)
    fprintf('--------------------------------------------------------\n')

    EEG = pop_biosig(filenames{iFile});
    EEG = pop_resample(EEG, 256);
    EEG = pop_select( EEG, 'chantype',{'EEG'});
    EEG = pop_chanedit(EEG, {'lookup','standard_1005.elc'});
    EEG = ref_infinity(EEG);
    
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'usefftfilt',1);
    EEG = pop_eegfiltnew(EEG, 'locutoff',59,'hicutoff',61,'revfilt',1,'usefftfilt',1);

    EEG = ascent_compute(EEG, 'measure', 'SampEn', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'FuzzEn', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'FracDim', 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'MSE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'mMSE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'ExSEnt', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'MFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);
    EEG = ascent_compute(EEG, 'measure', 'RCMFE', 'num_scales', num_scales, 'coarsing', coarsing, 'vis', false);

    SampEn(:,iFile) = EEG.ascent.SampEn.data;
    ExSEnt(:,iFile) = EEG.ascent.ExSEnt.data.HDA;
    FuzzEn(:,iFile) = EEG.ascent.FuzzEn.data;
    FracDim(:,iFile) = EEG.ascent.FracDim.data;
    MSE(:,:,iFile) = EEG.ascent.MSE.data;
    mMSE(:,:,iFile) = EEG.ascent.mMSE.data;
    MFE(:,:,iFile) = EEG.ascent.MFE.data;
    RCMFE(:,:,iFile) = EEG.ascent.MSE.data;

    chanlocs = EEG.chanlocs;
    scales = EEG.ascent.MSE.scales;
    scales_bounds = EEG.ascent.mMSE.scales;

    save(fullfile(data_path, sprintf('ascent_outputs_EO_%s.mat', coarsing)), ...
        "SampEn", "ExSEnt", "FuzzEn", "FracDim", ...
        "MSE", "MFE", "mMSE", "RCMFE", ...
        'chanlocs', 'scales', 'scales_bounds')

end

disp("Done computing on the whole group for eyes-open condition!")


%% Stats

cd(data_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
coarsing = 'median';
nPerm = 1000;          % number of permutations for H0
alpha = 0.05;           % NaN to show maps of p-values
ct = 'mean';            % central tendency method (default = 'trimmed mean')
grp_type = 'dpt';       % groupe is dependent ('dpt') or independent ('idpt')
mcc_type = 2;           % multiple comparison correction method (0 = uncorrected; 1 = t-max; 2 = cluster; 3 = TFCE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('colorbrewer.mat'); colormap = colorbrewer.qual.Set3{12}./255; colormap(8,:) = [];

% [neighbors, neighbormatrix] = get_channelneighbors(chanlocs);

load(fullfile(data_path, sprintf('ascent_outputs_EC_%s.mat', coarsing)))
SampEn1 = SampEn;
FuzzEn1 = FuzzEn;
ExSEnt1 = ExSEnt;
FracDim1 = FracDim;
MSE1 = MSE;
MFE1 = MFE;
mMSE1 = mMSE;
RCMFE1 = RCMFE;

load(fullfile(data_path, sprintf('ascent_outputs_EO_%s.mat', coarsing)))
SampEn2 = SampEn;
FuzzEn2 = FuzzEn;
ExSEnt2 = ExSEnt;
FracDim2 = FracDim;
MSE2 = MSE;
MFE2 = MFE;
mMSE2 = mMSE;
RCMFE2 = RCMFE;

% SampEn
figure('color','w'); 
subplot(2,1,1); title('Mean across channels')
rm_raincloud({trimmean(SampEn1,20,1) trimmean(SampEn2,20,1)}, colormap([5 4],:)); 
legend('EC','','EO',''); xlabel('SampEn') % title('SampEn')
subplot(2,1,2); title('Mean across subjects')
rm_raincloud({trimmean(SampEn1,20,2) trimmean(SampEn2,20,2)}, colormap([5 4],:)); 
legend('EC','','EO',''); xlabel('SampEn') % title('SampEn')
set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(PSD_WHM_SCALP, PSD_BSL_SCALP, nboots, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(SampEn1, SampEn2, nPerm, ct, grp_type);
[mask, clust_summary] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, 1, 'time', chanlocs, grp_type, {size(SampEn1,2)});

plot_results('frequency', 'scalp', [], tvals, mask, chanlocs, 'main')

saveas(gcf, fullfile(outputs_path,'2vs1_novices.fig'));
print(gcf, fullfile(outputs_path,'2vs1_novices.png'),'-dpng','-r300');   % 300 dpi .png
writetable(summary_tbl, fullfile(outputs_path,'2vs1_novices_summary.csv'));


% RCMFE
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(RCMFE1, RCMFE2, nPerm, ct, grp_type);
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_permutation(RCMFE1, RCMFE2, nPerm, ct, grp_type);
% [mask, clust_summary] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, scales, 'time', chanlocs, grp_type, {size(SampEn1,2) size(SampEn1,2)});
[mask, pcorr, tfce_score, tfce_H0_score] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcc_type, alpha, chanlocs);
[mask_clusters, clust_summary] = pull_clusters(mask, tvals, scales, chanlocs, 'frequency', grp_type, {size(MSE1,3) size(MSE2,3)}, 2, [], true, 'cohen-d');
plot_results('frequency', 'scalp', scales, tvals, mask_clusters, chanlocs, 'all', clust_summary)
% saveas(gcf, fullfile(outputs_path,'2vs1_novices.fig'));
% print(gcf, fullfile(outputs_path,'2vs1_novices.png'),'-dpng','-r300');   % 300 dpi .png
% writetable(clust_summary, fullfile(outputs_path,'2vs1_novices_summary.csv'));


% Course plot + topo of each sifgnificant cluster
for iClust = 1:size(clust_summary,1)

    % Topoplot
    peak_scale = clust_summary.Peak(iClust);        
    [~, idx] = min(abs(scales - peak_scale));           % closest frequency index
    peak_channel = clust_summary.Channel{iClust};      
    chan_idx = find(strcmpi({chanlocs.labels}, peak_channel));  % index of peak channel
    topo_mask = mask_clusters{iClust}(:, idx); % Logical mask of significant electrodes at peak frequency 
    figure('Color','w');
    topoplot(tvals(:, idx), chanlocs, 'verbose','off','colormap',dmap, ...
        'whitebk','on', 'pmask', topo_mask);
    title(sprintf('Cluster %g (Scale factor: %g)', iClust, peak_scale), ...
        'FontSize', 15, 'FontWeight', 'bold');
    cb = colorbar(); ylabel(cb,'T-values','Rotation',270,'fontweight','bold','fontsize',15); set(gcf,'color','w'); 
    set(findall(gcf,'type','axes'),'fontSize',20,'fontweight','bold');
    % saveas(gcf, fullfile(outputs_path, sprintf('2vs1_all_cluster%g-topoplot.fig', iClust)));
    % print(gcf, fullfile(outputs_path, sprintf('2vs1_all_cluster%g-topoplot.png', iClust)), '-dpng', '-r300');
    
    % Course plot (mean + 95% HDI)
    data1 = squeeze(RCMFE1(chan_idx,:,:));  % condition 1
    data2 = squeeze(RCMFE2(chan_idx,:,:));  % condition 2
    mask_chan = mask_clusters{iClust}(chan_idx, :); % Logical mask of significant frequencies at peak channel 
    figure('Color','w');
    plotDiff(scales, data2, data1, 'mean', 'CI', mask_chan, 'EC', 'EO');  % now using correct 1D mask
    title(sprintf('Cluster %g, Peak: %s at Scale factor %.0f ', iClust, peak_channel, peak_scale), ...
        'FontSize', 15, 'FontWeight', 'bold');
    ylabel('Difference (entropy)'); xlabel('Scale factors');
    set(findall(gcf,'type','axes'),'FontSize',20,'FontWeight','bold');
    axis tight
    % saveas(gcf, fullfile(outputs_path, sprintf('2vs1_cluster%g-psdplot.fig', iClust)));
    % print(gcf, fullfile(outputs_path, sprintf('2vs1_cluster%g-psdplot.png', iClust)), '-dpng', '-r300');

end


