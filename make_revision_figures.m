% make_revision_figures.m - Regenerate the manuscript figures that changed in
% the Sept 2026 revision, with the current ASCENT code. Run from the MATLAB
% desktop (figure rendering hangs under matlab -batch on this machine).
%
%   Fig. 4 : channel-domain aperiodic, static + time-resolved
%   Fig. 5 : ICA-domain aperiodic, static + time-resolved
%   Fig. 6 : single-scale EC vs EO topographies (ExSEnt HDA now corrected
%            with its own cluster mask; the old figure reused the HA mask)
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
for v = 1:numel(vars)
    x1 = EC.(vars{v})(order_idx, :); x2 = EO.(vars{v})(order_idx, :);
    [tvals, pvals, tvals_H0, pvals_H0] = run_stats_permutation(x1, x2, 2000, 'mean', 'dpt');
    mask = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 2, 0.05, chanlocs);   % cluster-based
    subplot(2, 4, v)
    topoplot(tvals, chanlocs, 'colormap', dmap, 'pmask', mask, 'verbose', 'off', 'whitebk', 'on');
    c = colorbar; ylabel(c, 't-values', 'FontWeight', 'bold', 'FontSize', 12)
    title(titles{v})
end
set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold', 'CLim', [-10 10]);
print(gcf, fullfile(outdir, 'fig_uniscales.png'), '-dpng', '-r300');

fid = fopen(fullfile(outdir, 'DONE.txt'), 'w'); fprintf(fid, '%s\n', datestr(now)); fclose(fid);
disp('All revision figures written to manuscript/figures/regen')
