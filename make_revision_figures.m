% make_revision_figures.m - Regenerate the manuscript figures that changed in
% the Sept 2026 revision, with the current ASCENT code. Run from the MATLAB
% desktop (figure rendering hangs under matlab -batch on this machine).
%
%   Fig. 4 : channel-domain aperiodic, static + time-resolved
%   Fig. 5 : ICA-domain aperiodic, static + time-resolved
%   Fig. 6 : single-scale EC vs EO topographies (ExSEnt HDA now corrected
%            with its own cluster mask; the old figure reused the HA mask)
%   Fig. 7 : Spearman correlation matrices (EC-EO differences; t-maps)
%   Figs. 10-11 : raw and aperiodic-corrected PSD, EC vs EO
%   (group aperiodic outputs were recomputed with the current code and
%   written into the cached ascent_outputs_*_new.mat files on 2026-09-25)
%
% Output: manuscript/figures/regen/*.png (300 dpi). Panels are composed into
% the final multi-panel figures outside MATLAB.

root    = fileparts(mfilename('fullpath'));
eegRoot = fullfile(fileparts(root), 'eeglab');
statRoot = fullfile(fileparts(root), 'eeg_robust_statistics');
outdir  = fullfile(root, 'manuscript', 'figures', 'regen');
if ~exist(outdir, 'dir'), mkdir(outdir); end

addpath(genpath(eegRoot)); addpath(root); addpath(fullfile(root, 'functions'));
% EEGLAB plugin compat stubs shadow MATLAB built-ins (pwelch, mean, ...)
bad = {fullfile(eegRoot,'plugins','Biosig3.8.5'), fullfile(eegRoot,'plugins','Fieldtrip-lite250523')};
for k = 1:numel(bad), if exist(bad{k},'dir'), rmpath(genpath(bad{k})); end, end

%% ---- Figs. 4-5: aperiodic, sample dataset ----
scratch = fullfile(root, 'scratch_aperiodic');
if ~exist(scratch, 'dir'), mkdir(scratch); end
copyfile(fullfile(root, 'ascent_sample_data.set'), fullfile(scratch, 'ascent_sample_data.set'));
EEG = pop_loadset('filename', 'ascent_sample_data.set', 'filepath', scratch);

for dom = {'channel', 'ica'}
    close all
    ascent_compute(EEG, 'measure', 'Aperiodic', 'domain', dom{1}, 'timeResolved', true, ...
        'vis', true, 'progress', false, 'parallel', false);
    figs = findobj('Type', 'figure');
    for k = 1:numel(figs)
        nm = get(figs(k), 'Name');
        if contains(nm, 'Aperiodic time course'), suffix = '_timecourse';
        elseif contains(nm, 'Aperiodic visualization'), suffix = '_static';
        else, continue
        end
        print(figs(k), fullfile(outdir, sprintf('aperiodic_%s%s.png', dom{1}, suffix)), '-dpng', '-r300');
    end
end
close all

%% ---- Fig. 6: single-scale group topographies ----
addpath(genpath(statRoot));
data_path = 'C:\Users\ccann\Documents\biosemi_data';
S = load(fullfile(data_path, 'ascent_outputs_EC_sd_new.mat'), 'chanlocs');
chanlocs = S.chanlocs;
desired_order = {'Fp1','FPz','FP2','AF7','AF3','AFz','AF4','AF8','F7','F5','F3','F1','Fz','F2','F4','F6','F8', ...
    'FT7','FT8','FC5','FC3','FC1','FCz','FC2','FC4','FC6','T7','T8','C5','C3','C1','Cz','C2','C4','C6', ...
    'TP7','TP8','CP5','CP3','CP1','CPz','CP2','CP4','CP6','P9','P7','P5','P3','P1','Pz','P2','P4','P6','P8','P10', ...
    'PO7','PO3','POz','PO4','PO8','O1','Oz','Iz','O2'};
[~, loc] = ismember(lower({chanlocs.labels}), lower(desired_order));
[~, order_idx] = sort(loc, 'ascend');
chanlocs = chanlocs(order_idx);

vars   = {'SampEn','FuzzEn','ExSEnt1','ExSEnt2','ExSEnt3','FracDim','Exponent','Offset'};
titles = {'SampEn','FuzzEn','ExSEnt (duration)','ExSEnt (amplitude)','ExSEnt (amp + dur)', ...
          'HigFracDim','Aperiodic Exponent','Aperiodic Offset'};
EC = load(fullfile(data_path, 'ascent_outputs_EC_sd_new.mat'), vars{:});
EO = load(fullfile(data_path, 'ascent_outputs_EO_sd_new.mat'), vars{:});
rdbu = interp1([1;128;256], [0.698 0.094 0.168; 1 1 1; 0.129 0.400 0.675], (1:256)', 'linear');
dmap = flipud(max(0, min(1, rdbu)));

figure('Name', 'Fig 6', 'Color', 'w', 'Position', [50 50 1600 900]);
tvals_all = zeros(numel(chanlocs), numel(vars));
for v = 1:numel(vars)
    x1 = EC.(vars{v})(order_idx, :); x2 = EO.(vars{v})(order_idx, :);
    rng(1);
    [tvals, pvals, tvals_H0, pvals_H0] = run_stats_permutation(x1, x2, 2000, 'mean', 'dpt');
    tvals_all(:, v) = tvals(:);
    mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 2, 0.05, chanlocs);   % cluster-based
    subplot(2, 4, v)
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, 'verbose', 'off', 'whitebk', 'on');
    c = colorbar; ylabel(c, 't-values', 'FontWeight', 'bold', 'FontSize', 12)
    title(titles{v})
end
set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold', 'CLim', [-10 10]);
print(gcf, fullfile(outdir, 'fig_uniscales.png'), '-dpng', '-r300');

%% ---- Fig. 7: correlation matrices ----
names7 = {'SampEn','FuzzEn','ExSEnt (HD)','ExSEnt (HA)','ExSEnt (HDA)','HigFracDim','Aperiodic Exponent','Aperiodic Offset'};
avgD = zeros(size(EC.SampEn, 2), numel(vars));
for v = 1:numel(vars), avgD(:, v) = mean(EC.(vars{v}) - EO.(vars{v}), 1, 'omitnan')'; end
R_global  = corr(avgD, 'Type', 'Spearman', 'Rows', 'pairwise');
R_spatial = corr(tvals_all, 'Type', 'Spearman', 'Rows', 'pairwise');
plot_corrmatrix_with_text(R_global, names7, 'Global Spearman correlations of EC-EO differences (channel-averaged)', fullfile(outdir, 'fig7A_corr_diff'));
plot_corrmatrix_with_text(R_spatial, names7, 'Spatial correlation of t-maps (EC vs EO)', fullfile(outdir, 'fig7B_corr_tmaps'));
writematrix(round(R_global, 2), fullfile(outdir, 'fig7A_values.csv'));
writematrix(round(R_spatial, 2), fullfile(outdir, 'fig7B_values.csv'));

%% ---- Figs. 10-11: raw and aperiodic-corrected PSD ----
P1 = load(fullfile(data_path, 'ascent_outputs_EC_sd_new.mat'), 'PSD', 'PSD_corr', 'freqs');
P2 = load(fullfile(data_path, 'ascent_outputs_EO_sd_new.mat'), 'PSD', 'PSD_corr');
freqs = P1.freqs;
psdSets = {'PSD', 'PSD_raw', 'Eyes closed vs. Eyes open - Classic PSD', 1; ...
           'PSD_corr', 'PSD_corrected', 'Eyes closed vs. Eyes open - Aperiodic-corrected PSD', 0.5};
for k = 1:2
    A = 10*log10(P1.(psdSets{k,1})(order_idx,:,:)); B = 10*log10(P2.(psdSets{k,1})(order_idx,:,:));
    rng(1);
    [tvals, pvals, tvals_H0, pvals_H0] = run_stats_permutation(A, B, 2000, 'mean', 'dpt');
    mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 2, 0.05, chanlocs);
    [mask_clusters, summary_tbl] = pull_clusters(mask, tvals, freqs, chanlocs, ...
        'frequency', 'dpt', {size(A,3) size(B,3)}, psdSets{k,4}, [], [], 'g');
    writetable(summary_tbl, fullfile(outdir, [psdSets{k,2} '_summary.csv']));
    plot_results('frequency', 'scalp', freqs, tvals, mask_clusters, chanlocs, 'main', summary_tbl);
    title(psdSets{k,3});
    set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');
    print(gcf, fullfile(outdir, [psdSets{k,2} '_main.png']), '-dpng', '-r300');
    hs = plot_clusters(summary_tbl, mask_clusters, tvals, tvals, freqs, chanlocs, ...
        'Power (dB)', 'DataType', 'scalp', 'Domain', 'Frequency');
    for i = 1:numel(hs.curve)
        print(hs.topo{i},  fullfile(outdir, sprintf('%s_cluster-%g_topo.png', psdSets{k,2}, i)), '-dpng', '-r300');
        print(hs.curve{i}, fullfile(outdir, sprintf('%s_cluster-%g_curve.png', psdSets{k,2}, i)), '-dpng', '-r300');
    end
end

fid = fopen(fullfile(outdir, 'DONE.txt'), 'w'); fprintf(fid, '%s\n', datestr(now)); fclose(fid);
disp('All revision figures written to manuscript/figures/regen')

%% ---- local functions (from ascent_group_analysis.m) ----
function plot_corrmatrix_with_text(C, labels, title_str, save_path)
cmap = rdbu_colormap();
C_full = C; C_full(triu(true(size(C)), 1)) = 0;
nVar = size(C_full, 1);
[xAll, yAll] = meshgrid(1:nVar, 1:nVar);
xAll((C_full == 0) & ~eye(nVar)) = nan;
clrLim = [-1, 1];
colIdx = discretize((C_full - clrLim(1)) / range(clrLim), linspace(0, 1, size(cmap,1)));
diamLim = [0.3, 0.9];
diamSize = abs(C_full) * range(diamLim) + diamLim(1);
fh = figure('Color','w', 'Position', [100 100 1000 700]);
ax = axes(fh); hold(ax, 'on'); colormap(cmap);
tv = 1:nVar;
text(-1*ones(size(tv)), tv, labels, 'HorizontalAlignment','right', 'FontSize',12, 'FontWeight','bold', 'Color','k');
text(tv, (nVar+1)*ones(size(tv)), labels, 'HorizontalAlignment','right', 'Rotation',45, 'FontSize',12, 'FontWeight','bold', 'Color','k');
theta = linspace(0, 2*pi, 100);
for i = 1:numel(xAll)
    if isnan(xAll(i)), continue; end
    fill(diamSize(i)/2*cos(theta) + xAll(i), diamSize(i)/2*sin(theta) + yAll(i), cmap(colIdx(i),:), 'LineStyle','none');
    text(xAll(i), yAll(i), sprintf('%.2f', C_full(i)), 'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', 'FontSize', 9, 'Color', 'k', 'FontWeight','bold');
end
set(ax, 'YDir', 'reverse'); xlim([-2, nVar+1]); ylim([0.5, nVar+2]);
cb = colorbar; ylabel(cb, 'Correlation coefficient', 'FontSize', 12, 'FontWeight', 'bold');
clim(clrLim); axis off;
set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');
title(title_str, 'FontSize', 13, 'FontWeight', 'bold');
print(gcf, [save_path '.png'], '-dpng', '-r300');
end

function cmap = rdbu_colormap()
n = 256; half = floor(n/2);
cmap = [linspace(0.019,0.97,half)' linspace(0.188,0.97,half)' linspace(0.380,0.97,half)'; ...
        linspace(0.97,0.698,n-half)' linspace(0.97,0.094,n-half)' linspace(0.97,0.169,n-half)'];
end
