function ascent_plot_multiscale(entropyData, chanlocs, entropyType, scales, opt)
% ascent_plot_multiscale  Entropy / complexity figures for Ascent.
%
% Called by ascent_plot; not intended for direct use. See ascent_plot for the
% public interface and for what opt carries. Handles every measure except the
% aperiodic ones (those go to ascent_plot_aperiodic).
%
% Four layouts, chosen from the shape of entropyData:
%   multiscale + multi-signal : heatmap with cluster overlay | per-cluster topo+curve
%   uniscale   + multi-signal : single topo
%   multiscale + one signal   : single curve
%   uniscale   + one signal   : scalar bar
%
% In ICA mode, topos are the cluster-mean entropy back-projected to the scalp
% through icawinv, and chanlocs are the channels the decomposition ran on.
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).

isICA         = opt.isICA;
icawinv       = opt.icawinv;
clusterThresh = opt.clusterThresh;
dipfit        = opt.dipfit;
iclabel       = [];
if isfield(opt,'iclabel'), iclabel = opt.iclabel; end
if isfield(opt,'fontSize') && ~isempty(opt.fontSize), fs = opt.fontSize; else, fs = 9; end
lw            = max(1, fs/9);   % line widths grow with the font for small prints

isTimeResolved = ndims(entropyData)==3;
if isTimeResolved
    entropyData = mean(entropyData,3,'omitnan');
end
multiscale = size(entropyData,2)>1;
if ~multiscale, entropyData(entropyData==0)=NaN; end

nChan     = size(entropyData,1);
multiChan = nChan>1;

if multiscale && multiChan
    % -------------------------------------------------------
    % Multiscale + multichannel: heatmap with the cluster outlines, plus
    % panels for each cluster (up to 3). Layout:
    %   one cluster, channels : heatmap left; topo top right, curve below it
    %   one cluster, ICs      : heatmap left; IC bar top right, curve | topo below it
    %   2-3 clusters          : heatmap left, one topo | curve row per cluster
    % -------------------------------------------------------

    finiteMask=isfinite(entropyData);
    if ~any(finiteMask(:))
        warning('ascent_plot: all values are NaN; nothing to plot.'); return
    end

    % Cluster detection (up to 3)
    [clusterMasks,nClusters,centroids,usedCluster] = ...
        find_entropy_clusters(entropyData, clusterThresh, 3);
    nShow   = min(nClusters,3);
    compact = nShow == 1;

    % ICLabel class and equivalent dipole of each IC (ICA only)
    icinfo = [];
    if isICA && (~isempty(iclabel) || ~isempty(dipfit))
        icinfo = ascent_ic_info(nChan, iclabel, dipfit, chanlocs);
    end
    isBrainIC = [];
    if ~isempty(icinfo) && ~all(cellfun(@isempty, {icinfo.iclabel}))
        isBrainIC = strcmp({icinfo.iclabel}, 'Brain')';
    end

    % Shared cluster colours: used for both heatmap contours and curve plots.
    % Chosen to be visible on the dark parula heatmap AND on white axes.
    % No yellow (invisible on white) or blue (clashes with parula).
    clrCluster = {[1.00 0.35 0.10], ...   % C1: vermillion-orange
                  [0.12 0.78 0.35], ...   % C2: bright green
                  [0.82 0.15 0.82]};      % C3: magenta

    hFig2=figure('Color','w','InvertHardCopy','off',...
                 'Name','Multiscale entropy visualization',...
                 'Toolbar','none','Menu','none','NumberTitle','Off');
    try icadefs; set(hFig2,'color',BACKCOLOR); catch; end %#ok<NODEF>
    if compact
        % One cluster: a smaller figure, so the cluster panels do not leave
        % two empty rows
        if isICA, sz = [1000 600]; else, sz = [1000 540]; end
        set(hFig2,'Units','pixels'); p0 = get(hFig2,'Position');
        set(hFig2,'Position',[p0(1) max(40, p0(2)+p0(4)-sz(2)) sz]);
    end

    nScales = size(entropyData,2);
    xvals   = ifelse(isnumeric(scales), scales, 1:nScales);

    % Compact layouts use OuterPosition so that MATLAB fits tick labels, axis
    % labels and the colorbar inside each panel (subplot slots let a large font
    % run into the neighbouring panel)
    if compact
        ax_heat = axes(hFig2,'OuterPosition',[0 0 0.56 1]);
    else
        ax_heat = subplot(3,4,[1 2 5 6 9 10]);
    end
    imagesc(ax_heat,1:nScales,1:nChan,entropyData); axis(ax_heat,'tight');
    set(ax_heat,'TickDir','out','Box','off');   % box ticks would poke into the title
    ax_heat.TickLabelInterpreter='none';
    colormap(ax_heat,'parula');
    ch=colorbar(ax_heat); ch.Label.String='Entropy';
    xlabel(ax_heat,'Scales');
    ylabel(ax_heat,ifelse(isICA,'ICs','EEG channels'));
    title(ax_heat,entropyType,'Interpreter','none');
    set(ax_heat,'FontSize',fs,'FontWeight','bold'); ch.FontSize=fs;
    set_heat_ticks(ax_heat, scales, nChan, isICA, chanlocs, fs);

    % Overlay cluster boundaries using the shared cluster colours. The outline
    % follows the mask exactly: the member signals over the cluster's scales.
    % Rows are signals in list order, so members need not be adjacent and one
    % cluster can be several boxes; the old smoothing (sigma 1.5 cells, level
    % 0.35) dropped isolated rows and drew a far smaller cluster than reported.
    % The mask is upsampled f times and smoothed by a quarter cell, which only
    % rounds the corners.
    if usedCluster
        f  = 8;
        g  = exp(-((-f:f).^2) / (2*(f/4)^2));  g = g / sum(g);
        xs = ((1:(nScales+2)*f) - 0.5)/f - 0.5;   % cell centres at 1..n after a
        ys = ((1:(nChan+2)*f)   - 0.5)/f - 0.5;   % one-cell zero border
        hold(ax_heat,'on');
        for k = 1:nShow
            U = kron(padarray0(double(clusterMasks(:,:,k))), ones(f));
            contour(ax_heat, xs, ys, conv2(g, g, U, 'same'), [0.5 0.5], ...
                    'Color', clrCluster{k}, 'LineWidth', 1.2*lw);
        end
        hold(ax_heat,'off');
    end

    % Shade scale extent on curves only when clusters differ meaningfully in
    % scale coverage -- compare first/last occupied scale per cluster,
    % tolerating up to 1-step differences from boundary effects.
    doShading = false;
    if usedCluster && nShow > 1
        sc_lo = zeros(1, nShow);
        sc_hi = zeros(1, nShow);
        for k = 1:nShow
            ii = find(any(clusterMasks(:,:,k), 1));
            sc_lo(k) = ii(1); sc_hi(k) = ii(end);
        end
        doShading = (range(sc_lo) > 1) || (range(sc_hi) > 1);
    end

    % Right-side subplot slot indices in the 3x4 grid (2-3 clusters)
    % Row 1: 3,4 | Row 2: 7,8 | Row 3: 11,12
    topoSlots  = [3  7  11];
    curveSlots = [4  8  12];
    unit = ifelse(isICA, 'ICs', 'ch');
    ic_idx_1 = [];   % cluster 1 members, for the ICLabel count

    for k = 1:ifelse(compact, 1, 3)
        if compact && ~isICA
            ax_topo  = axes(hFig2,'OuterPosition',[0.56 0.5 0.44 0.5]);
            ax_curve = axes(hFig2,'OuterPosition',[0.56 0   0.44 0.5]);
        elseif compact
            ax_bar   = axes(hFig2,'OuterPosition',[0.56 0.5 0.44 0.5]);
            ax_curve = axes(hFig2,'OuterPosition',[0.56 0   0.22 0.5]);
            ax_topo  = axes(hFig2,'OuterPosition',[0.78 0   0.22 0.5]);
        else
            ax_topo  = subplot(3,4, topoSlots(k));
            ax_curve = subplot(3,4, curveSlots(k));
            if k > nShow
                axis(ax_topo,'off'); axis(ax_curve,'off'); continue
            end
        end

        mask_k      = clusterMasks(:,:,k);
        scaleMask_k = any(mask_k,1);
        chanMask_k  = any(mask_k,2);
        topo_vals_k = mean(entropyData(:, scaleMask_k), 2, 'omitnan');
        curve_k     = mean(entropyData(chanMask_k, :), 1, 'omitnan');
        nCh_k = sum(chanMask_k);
        nSc_k = sum(scaleMask_k);
        pc_k  = centroids(k,1);
        ps_k  = centroids(k,2);
        if iscell(scales), scLabel_k = scales{ps_k}; else, scLabel_k = num2str(scales(ps_k)); end
        if isICA, sigName = sprintf('IC%d', pc_k); else, sigName = chanlocs(pc_k).labels; end
        sc_idx_k = find(scaleMask_k);
        if isempty(sc_idx_k)
            scRange = scLabel_k;
        elseif isnumeric(scales)
            scRange = sprintf('%g-%g', scales(sc_idx_k(1)), scales(sc_idx_k(end)));
        else
            scRange = sprintf('%s-%s', scales{sc_idx_k(1)}, scales{sc_idx_k(end)});
        end
        if k == 1, ic_idx_1 = find(chanMask_k); end

        % Console summary
        if usedCluster
            fprintf('Cluster %d: %d %s x %d scales  (centroid %s, scale %s)\n', ...
                    k, nCh_k, unit, nSc_k, sigName, scLabel_k);
        elseif k == 1
            fprintf('No cluster found - peak: %s, scale %s\n', sigName, scLabel_k);
        end

        % Panel titles: short in the compact layouts; two lines in the narrow
        % IC curve and topo panels
        if compact && isICA
            topoTitle = {'Back-projected', ['scales ' scRange]};
            if usedCluster
                curveTitle = {'Cluster mean', sprintf('(%d ICs)', nCh_k)};
            else
                curveTitle = {'Peak', sigName};
            end
        elseif compact
            topoTitle = sprintf('Scales %s, peak %s sc %s', scRange, sigName, scLabel_k);
            curveTitle = ifelse(usedCluster, sprintf('Cluster mean (%d ch)', nCh_k), ['Peak ' sigName]);
        else
            if usedCluster
                topoTitle = sprintf('Cluster %d  (%d %s x %d sc, %s sc%s)', k, nCh_k, unit, ...
                                    nSc_k, sigName, scLabel_k);
                curveTitle = sprintf('Cluster %d mean  (%d %s)', k, nCh_k, unit);
            else
                topoTitle  = sprintf('Peak  %s, sc %s', sigName, scLabel_k);
                curveTitle = 'Peak curve';
            end
        end

        % IC bar (one cluster, ICs)
        if compact && isICA
            bar(ax_bar, topo_vals_k, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none');
            hold(ax_bar,'on');
            if usedCluster
                bar(ax_bar, ic_idx_1, topo_vals_k(ic_idx_1), ...
                    'FaceColor', clrCluster{k}, 'EdgeColor', 'none');
            end
            bar(ax_bar, pc_k, topo_vals_k(pc_k), 'FaceColor', 'r', 'EdgeColor', 'none');
            if ~isempty(isBrainIC)   % dot above the ICs ICLabel classifies as brain
                ib = find(isBrainIC);
                plot(ax_bar, ib, topo_vals_k(ib) + 0.04*max(topo_vals_k), 'o', 'LineStyle', 'none', ...
                     'MarkerSize', 2.5*lw, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
            end
            xlim(ax_bar, [0 nChan+1]);
            xlabel(ax_bar,'IC'); ylabel(ax_bar,'Entropy'); box(ax_bar,'on');
            if usedCluster
                barTitle = sprintf('%d ICs x %d scales, peak %s sc %s', nCh_k, nSc_k, sigName, scLabel_k);
            else
                barTitle = sprintf('Peak %s, sc %s', sigName, scLabel_k);
            end
            if ~isempty(isBrainIC)
                if usedCluster
                    barTitle = {barTitle, sprintf('%s ICLabel brain: %d of %d', char(8226), ...
                                sum(isBrainIC(ic_idx_1)), nCh_k)};
                else
                    barTitle = {barTitle, sprintf('%s ICLabel brain: %d of %d ICs', char(8226), ...
                                sum(isBrainIC), nChan)};
                end
            end
            title(ax_bar, barTitle, 'Interpreter','none');
        end

        % Topo: channel entropy, or for ICs the cluster mean back-projected
        % to the scalp through icawinv
        if isICA, tv = icawinv * topo_vals_k; else, tv = topo_vals_k; end
        plot_topo(ax_topo, tv, chanlocs, max(7, round(7*lw)));
        title(ax_topo, topoTitle, 'Interpreter','none');

        % Curve
        hold(ax_curve,'on'); box(ax_curve,'on');
        if all(~isfinite(curve_k)), curve_k = nan(1,nScales); end
        if doShading && isnumeric(scales) && any(scaleMask_k)
            x_lo = xvals(find(scaleMask_k,1,'first'));
            x_hi = xvals(find(scaleMask_k,1,'last'));
            finC = curve_k(isfinite(curve_k));
            if ~isempty(finC)
                yl_pre = [min(finC)*0.95 max(finC)*1.05];
                if yl_pre(1)==yl_pre(2), yl_pre = yl_pre+[-0.1 0.1]; end
                patch(ax_curve,[x_lo x_hi x_hi x_lo],...
                      [yl_pre(1) yl_pre(1) yl_pre(2) yl_pre(2)],...
                      clrCluster{k},'FaceAlpha',0.15,'EdgeColor','none');
            end
        end
        plot(ax_curve, xvals, curve_k, 'LineWidth', 2*lw, 'Color', clrCluster{k});
        xlim(ax_curve,[xvals(1) xvals(end)]);
        xlabel(ax_curve,'Scale'); ylabel(ax_curve,'Entropy');
        title(ax_curve, curveTitle, 'Interpreter','none');
    end

    colormap(hFig2,parula);
    set(findall(hFig2,'type','axes'),'FontSize',fs,'FontWeight','bold', ...
        'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1);
    set(findall(hFig2,'type','colorbar'),'FontSize',fs);
    ch.Label.FontSize = fs; ch.Label.FontWeight = 'bold';

    % Region entropy + dipole sources figure (ICA + dipfit only)
    if isICA && ~isempty(dipfit)
        plot_regions(entropyData, xvals, icinfo, dipfit, entropyType, fs);
    end

elseif ~multiscale && multiChan
    % Uniscale topo
    hFig3=figure('Color','w','InvertHardCopy','off',...
                 'Name','Uniscale entropy visualization',...
                 'Toolbar','none','Menu','none','NumberTitle','Off');
    try icadefs; set(hFig3,'color',BACKCOLOR); catch; end %#ok<NODEF>
    vals=entropyData(:); finiteVals=vals(isfinite(vals));
    ax_uni=axes(hFig3);
    if isempty(finiteVals)
        axis(ax_uni,'off');
        text(0.5,0.5,'No finite values to plot','Parent',ax_uni,...
             'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold');
    else
        if isICA
            bp_vals=icawinv*vals(:); finBP=bp_vals(isfinite(bp_vals));
            axes(ax_uni); topoplot(bp_vals,chanlocs,'emarker',{'.','k',15,1},'electrodes','labels'); %#ok<LAXES>
            if min(finBP)<max(finBP), clim(ax_uni,[min(finBP)*0.95 max(finBP)*1.05]); end
            title(ax_uni,[entropyType ' (back-projected)'],'Interpreter','none');
        else
            axes(ax_uni); topoplot(vals,chanlocs,'emarker',{'.','k',15,1},'electrodes','labels'); %#ok<LAXES>
            if min(finiteVals)<max(finiteVals)
                clim(ax_uni,[min(finiteVals)*0.95 max(finiteVals)*1.05]);
            end
            title(ax_uni,entropyType,'Interpreter','none');
        end
        colormap(hFig3,parula);
        c=colorbar(ax_uni);
        c.Label.String='Entropy'; c.Label.FontSize=fs+2; c.Label.FontWeight='bold';
    end
    set(findall(hFig3,'type','axes'),'FontSize',fs+1,'FontWeight','bold');

elseif multiscale && ~multiChan
    % Single channel/IC curve
    hFig4=figure('Color','w','InvertHardCopy','off',...
                 'Name',entropyType,'Toolbar','none','Menu','none','NumberTitle','Off');
    try icadefs; set(hFig4,'color',BACKCOLOR); catch; end %#ok<NODEF>
    ax_s=axes(hFig4);
    xvals = ifelse(isnumeric(scales), scales, 1:numel(scales));
    plot(ax_s,xvals,entropyData,'LineWidth',2*lw);
    xlabel(ax_s,'Scale'); ylabel(ax_s,'Entropy');
    title(ax_s,entropyType,'Interpreter','none');
    axis(ax_s,'tight');
    set(findall(hFig4,'type','axes'),'FontSize',fs+1,'FontWeight','bold');

elseif ~multiscale && ~multiChan
    % Scalar
    hFig5=figure('Color','w','InvertHardCopy','off','Name',entropyType,...
                 'Toolbar','none','Menu','none','NumberTitle','Off','Position',[200 200 400 300]);
    try icadefs; set(hFig5,'color',BACKCOLOR); catch; end %#ok<NODEF>
    val=entropyData(isfinite(entropyData)); ax_sc=axes(hFig5);
    if isempty(val)
        axis(ax_sc,'off');
        text(0.5,0.5,sprintf('%s = NaN',entropyType),'Parent',ax_sc,...
             'HorizontalAlignment','center','VerticalAlignment','middle',...
             'FontSize',14,'FontWeight','bold','Interpreter','none');
    else
        bar(ax_sc,1,val,0.4,'FaceColor',[0.18 0.45 0.87],'EdgeColor','none');
        set(ax_sc,'XTick',1,'XTickLabel',{entropyType},'TickDir','out','TickLabelInterpreter','none');
        ylabel(ax_sc,'Entropy');
        title(ax_sc,sprintf('%s = %.4f',entropyType,val),'Interpreter','none','FontSize',13,'FontWeight','bold');
        ylim(ax_sc,[min(0,val*1.2) max(0,val*1.2)]);
        box(ax_sc,'on');
        if val<0
            text(1,val*0.5,'negative entropy is valid','Parent',ax_sc,...
                 'HorizontalAlignment','center','FontSize',9,'Color',[0.5 0.5 0.5],'Interpreter','none');
        end
    end
    set(findall(hFig5,'type','axes'),'FontSize',fs+1,'FontWeight','bold');

else
    error('ascent_plot: data format not recognized.')
end
end


%% Cluster detection (scale-first, channel-order agnostic)
function [clusterMasks, nClusters, centroids, usedCluster] = ...
        find_entropy_clusters(data, thresh_pct, maxClusters)
% find_entropy_clusters  Scale-first, channel-order-agnostic clustering.
%
% Strategy (avoids 2-D row-adjacency bias from arbitrary channel ordering):
%   1. Compute per-scale mean entropy across channels -> [1 x nScale] summary.
%   2. Threshold the summary at thresh_pct of all finite data values.
%   3. Find contiguous hot scale-ranges by 1-D connected components
%      (no channel dimension involved => no ordering bias).
%   4. For each scale-range: include channels whose mean entropy in that
%      range exceeds the threshold, evaluated independently per channel.
%   5. Sort resulting blobs by summed entropy mass; return top maxClusters.
%
% Falls back to the global peak cell when no valid blob is found
% (usedCluster = false).
%
% Outputs
%   clusterMasks  [nChan x nScale x nFound] logical
%   nClusters     number of blobs returned (always >= 1)
%   centroids     [nClusters x 2]  [weighted_channel, weighted_scale]
%   usedCluster   true when at least one real cluster was found

if nargin < 3, maxClusters = 3; end
[nChan, nScale] = size(data);
finiteMask      = isfinite(data);

% --- Fallback: global peak cell ---
tmp = data; tmp(~finiteMask) = -Inf;
[~, li] = max(tmp(:));
[peak_ch, peak_sc] = ind2sub([nChan nScale], li);
fallback = false(nChan, nScale);
fallback(peak_ch, peak_sc) = true;
clusterMasks = false(nChan, nScale, 1);
clusterMasks(:,:,1) = fallback;
nClusters = 1; centroids = [peak_ch peak_sc]; usedCluster = false;

if sum(finiteMask(:)) < 3, return; end

% Two separate thresholds, each relative to its own distribution:
%   thresh_scale  applied to the per-scale channel mean  -> selects hot scales
%   thresh_chan   applied to per-channel scale-range mean -> selects channels
% Using the mean's own distribution prevents the aggregation from always
% falling below a cell-level percentile (which would leave hot_scales empty).
scale_mean   = mean(data, 1, 'omitnan');                          % [1 x nScale]
thresh_scale = prctile(scale_mean(isfinite(scale_mean)), thresh_pct * 100);
thresh_chan  = prctile(data(finiteMask),                  thresh_pct * 100);

% --- Step 1: per-scale summary ---
hot_scales = scale_mean >= thresh_scale;   % [1 x nScale] logical

if ~any(hot_scales), return; end

% --- Step 2: 1-D connected scale-range segments ---
changes = diff([0 hot_scales 0]);
starts  = find(changes ==  1);
ends    = find(changes == -1) - 1;   % inclusive end
nRuns   = numel(starts);
if nRuns == 0, return; end

% --- Step 3: build per-segment channel masks and compute mass ---
blobMasks  = {};
blobMasses = [];
for r = 1:nRuns
    sc_idx = starts(r):ends(r);

    % Channel membership: mean entropy over this scale range >= cell threshold,
    % evaluated independently per channel (no channel-ordering bias)
    chan_mean = mean(data(:, sc_idx), 2, 'omitnan');   % [nChan x 1]
    chan_in   = chan_mean >= thresh_chan;
    if ~any(chan_in)
        % Fallback: include any channel with finite values in this range
        chan_in = any(finiteMask(:, sc_idx), 2);
    end
    if ~any(chan_in), continue; end

    mk = false(nChan, nScale);
    mk(chan_in, sc_idx) = true;
    mk = mk & finiteMask;
    if sum(mk(:)) < 3, continue; end

    blobMasks{end+1}  = mk; %#ok<AGROW>
    blobMasses(end+1) = sum(data(mk), 'omitnan'); %#ok<AGROW>
end

if isempty(blobMasks), return; end

% --- Step 4: sort by mass, keep top maxClusters ---
[~, order] = sort(blobMasses, 'descend');
nFound      = min(numel(order), maxClusters);
clusterMasks = false(nChan, nScale, nFound);
centroids    = zeros(nFound, 2);
usedCluster  = true;

for k = 1:nFound
    mk = blobMasks{order(k)};
    clusterMasks(:,:,k) = mk;

    % Centroid channel = IC/channel with highest mean entropy over the
    % cluster's scale extent (peak, not weighted mean).
    % Centroid scale   = scale with highest mean entropy across cluster channels.
    chan_in  = any(mk, 2);          % [nChan x 1]
    scale_in = any(mk, 1);          % [1 x nScale]
    chan_mean_k  = mean(data(:, scale_in), 2, 'omitnan');   % mean over cluster scales
    scale_mean_k = mean(data(chan_in, :),  1, 'omitnan');   % mean over cluster chans
    [~, pc] = max(chan_mean_k);
    [~, ps] = max(scale_mean_k);

    centroids(k,:) = [max(1,min(nChan,pc)), max(1,min(nScale,ps))];
end
nClusters = nFound;
end


%% Dipole sources + region entropy (ICA + dipfit)
function plot_regions(ent, scales, icinfo, dipfit, entropyType, fs, rv_thresh)
% Dipolar ICs (RV < rv_thresh) inside the brain, from ascent_ic_info: an MRI
% panel with the dipoles (colour: mean entropy across scales; size: fit
% quality, 1 - RV; circle: ICLabel brain, square: other ICLabel classes) and,
% below, entropy by AAL lobe for lobes with at least 3 ICs. Dipoles outside the
% brain are artifact ICs (eyes, muscle) and are left out.
if nargin < 7, rv_thresh = 0.15; end
lw = max(1, fs/9);
n_ics = size(ent, 1);
if isempty(icinfo), return; end
sel = icinfo([icinfo.dipolar] & [icinfo.inbrain] & [icinfo.ic] <= n_ics);
if isempty(sel)
    warning('ascent_plot: no dipolar IC (RV < %.0f%%) inside the brain; no dipole figure.', rv_thresh*100);
    return
end
ent_sel   = ent([sel.ic],:);
mean_ent  = mean(ent_sel, 2, 'omitnan');
hasICL    = ~all(cellfun(@isempty, {sel.iclabel}));
isBrainIC = strcmp({sel.iclabel}, 'Brain');

% Table in the command window
fprintf('Dipolar ICs inside the brain (RV < %.0f%%): %d\n', rv_thresh*100, numel(sel));
fprintf('   IC     RV   MNI x, y, z (mm)   AAL region (lobe)                      ICLabel (p)          Mean %s\n', entropyType);
for i = 1:numel(sel)
    s = sel(i);
    if isempty(s.aal), reg = 'unlabeled'; else, reg = sprintf('%s (%s)', s.aal, s.lobe); end
    if isempty(s.iclabel), icl = '-'; else, icl = sprintf('%s (%.2f)', s.iclabel, s.p_iclabel); end
    fprintf('  %3d  %4.1f%%  %5.0f %5.0f %5.0f    %-38s %-20s %.3f\n', s.ic, 100*s.rv, s.mni, reg, icl, mean_ent(i));
end

% Lobes with at least min_n ICs
lobe_names = {'Frontal','Central','Parietal','Temporal','Occipital','Cingulate','Insula','Subcortical','Cerebellum'};
lobe_clrs  = {[0.00 0.45 0.74],[0.90 0.55 0.00],[0.85 0.20 0.55],[0.45 0.30 0.85],[0.45 0.45 0.45], ...
              [0.15 0.15 0.55],[0.30 0.75 0.93],[0.00 0.00 0.00],[0.60 0.40 0.20]};
min_n = 3;
ic_lobes = {sel.lobe};
has_reg  = any(cellfun(@(r) sum(strcmp(ic_lobes, r)) >= min_n, lobe_names));

% ---- Run dipplot into temp figure, capture as image --------------------
dip_img = [];
vals = mean_ent';
cmap = plasma_map(256);   % purple to yellow, no green, readable on the grey MRI
cl = [min(vals) max(vals)]; if ~(cl(2) > cl(1)), cl = cl + [-1 1]*eps(cl(1)+1); end
if isfield(dipfit,'mrifile')
    try
        % Cluster membership is not marked: few high-entropy ICs are dipolar,
        % and a dipole position only means something for those that are.
        sources = struct('posxyz', {sel.mni}, 'momxyz', [0 0 0], 'rv', {sel.rv});
        ci = 1 + round((vals - cl(1)) / (cl(2) - cl(1)) * 255);
        colors = arrayfun(@(i) cmap(i,:), ci, 'UniformOutput', false);
        msz = 5 + 7 * (1 - min([sel.rv], rv_thresh) / rv_thresh);   % RV 0: 12 pt; RV at threshold: 5 pt
        mrk = repmat({'o'}, 1, numel(sel));
        if hasICL, mrk(~isBrainIC) = {'s'}; end
        pre_figs = findall(0,'Type','figure');
        dipplot(sources,'mri',template_mri_without_slab(dipfit),'coordformat','MNI','color',colors,...
            'summary','3d','dipolesize',15,'dipolelength',0,'gui','off','mesh','off');
        drawnow; pause(0.3);
        new_figs = setdiff(findall(0,'Type','figure'),pre_figs);
        if ~isempty(new_figs)
            tmp_fig = new_figs(end);
            hide_edge_on_slices(tmp_fig);
            outline_dipoles(tmp_fig, msz, mrk);
            set(tmp_fig,'Color','k','InvertHardcopy','off');
            try
                dip_img = print(tmp_fig,'-RGBImage','-r300');  % sharper than getframe
            catch
                dip_img = getframe(tmp_fig).cdata;
            end
            close(tmp_fig);
        end
    catch ME
        warning('ascent_plot: dipplot capture failed - %s', ME.message);
    end
end
has_dip = ~isempty(dip_img);
if ~has_dip && ~has_reg, return; end

% ---- Combined figure layout --------------------------------------------
fig = figure('Color','w','InvertHardCopy','off','NumberTitle','off',...
             'Name','Dipole sources + entropy by region');
try icadefs; set(fig,'color',BACKCOLOR); catch; end %#ok<NODEF>
set(fig,'Units','pixels'); p0 = get(fig,'Position');
if has_dip && has_reg
    sz = [1000 720]; rowH = 0.50; dipPos = [0.01 0.52 0.88 0.38];
elseif has_dip
    sz = [1000 380]; dipPos = [0.01 0.03 0.88 0.72];
else
    sz = [1000 420]; rowH = 1;
end
set(fig,'Position',[p0(1) max(40, p0(2)+p0(4)-sz(2)) sz]);
if has_dip
    ax_dip = axes(fig,'Position',dipPos);   % room on the right for the colour bar and its label
    image(ax_dip, dip_img); axis(ax_dip,'image','off');
    colormap(ax_dip, cmap); clim(ax_dip, cl);
    cb = colorbar(ax_dip, 'eastoutside'); cb.FontSize = fs; cb.FontWeight = 'bold';
    cb.Label.String = ['Mean ' entropyType]; cb.Label.FontSize = fs; cb.Label.FontWeight = 'bold';
    t2 = 'Size: fit quality (1 - RV)';
    if hasICL
        t2 = sprintf('%s; circles: ICLabel brain (%d), squares: other classes (%d)', t2, ...
                     sum(isBrainIC), sum(~isBrainIC));
    end
    title(ax_dip, {sprintf('Equivalent dipoles inside the brain (RV < %.0f%%, %d ICs), color: mean %s', ...
          rv_thresh*100, numel(sel), entropyType), t2}, 'FontSize',fs,'FontWeight','bold','Interpreter','none');
end
if ~has_reg
    fprintf('No AAL lobe has %d or more dipolar ICs: no region plot.\n', min_n);
    return
end
% OuterPosition: MATLAB fits the tick labels, axis labels and titles inside
ax1 = axes(fig,'OuterPosition',[0    0 0.64 rowH]);
ax2 = axes(fig,'OuterPosition',[0.64 0 0.36 rowH]);
hold(ax1,'on'); box(ax1,'on');
hold(ax2,'on'); box(ax2,'on');

xpos=0; xt_pos=[]; xt_lbl={}; leg_h=gobjects(0); allv=[];
for rr = 1:numel(lobe_names)
    rname = lobe_names{rr};
    idx   = strcmp(ic_lobes, rname);
    n_r   = sum(idx);
    if n_r < min_n, continue; end
    clr = lobe_clrs{rr};

    mu  = mean(ent_sel(idx,:),1,'omitnan');
    sem = std( ent_sel(idx,:),0,1,'omitnan') / sqrt(n_r);
    fill(ax1,[scales fliplr(scales)],[mu+sem fliplr(mu-sem)],...
         clr,'FaceAlpha',0.15,'EdgeColor','none');
    leg_h(end+1) = plot(ax1,scales,mu,'Color',clr,'LineWidth',2.2*lw,...
        'DisplayName',sprintf('%s (n=%d)',rname,n_r)); %#ok<AGROW>

    xpos = xpos+1;
    v = mean_ent(idx); allv = [allv; v]; %#ok<AGROW>
    scatter(ax2,xpos+(rand(n_r,1)-0.5)*0.3,v,30*lw^2,clr,...
            'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
    plot(ax2,xpos+[-0.3 0.3],repmat(mean(v,'omitnan'),1,2),...
         '-','Color',clr*0.6,'LineWidth',2.5*lw);
    xt_pos(end+1)=xpos; xt_lbl{end+1}=rname; %#ok<AGROW> (n is in the legend)
end

xlim(ax1,[scales(1) scales(end)]);
xlabel(ax1,'Scale'); ylabel(ax1,entropyType,'Interpreter','none');
title(ax1,'Scale-resolved entropy by AAL lobe  (mean \pm SEM)');
set(ax1,'TickDir','out');

xlim(ax2,[0 xpos+1]);
allv = allv(isfinite(allv));
if numel(allv) > 1 && max(allv) > min(allv)   % pad so no point sits on the frame
    pad = 0.08*(max(allv)-min(allv)); ylim(ax2,[min(allv)-pad max(allv)+pad]);
end
set(ax2,'XTick',xt_pos,'XTickLabel',xt_lbl,'TickDir','out','XTickLabelRotation',30);
ylabel(ax2,['Mean ' entropyType],'Interpreter','none');
title(ax2,'Mean entropy by AAL lobe');
set([ax1 ax2],'FontSize',fs,'FontWeight','bold', ...
    'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1);
legend(ax1,leg_h,'Location','best','Box','off','FontSize',max(fs-1,8));
end


%% Purple-to-yellow colour map (anchors of matplotlib's plasma)
function c = plasma_map(n)
a = [0.050 0.030 0.528; 0.417 0.000 0.658; 0.692 0.165 0.565; 0.881 0.393 0.383; ...
     0.988 0.652 0.211; 0.940 0.975 0.131];
c = interp1(linspace(0,1,size(a,1)), a, linspace(0,1,n));
end


%% Dipole markers with a black edge (dots vanish on the grey MRI)
function outline_dipoles(fig, msz, mrk)
% Size and shape per dipole; dipplot tags each marker 'dipole<index>'
for h = findall(fig, 'Type', 'line', 'Marker', '.')'
    t = get(h,'Tag');
    if ~strncmp(t, 'dipole', 6), continue; end
    i = str2double(t(7:end));
    if isnan(i) || i > numel(msz), continue; end
    c = get(h, 'Color');
    set(h, 'Marker', mrk{i}, 'MarkerSize', msz(i) * ifelse(strcmp(mrk{i},'s'), 0.9, 1), ...
        'MarkerFaceColor', c, 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1);
end
end


%% dipplot '3d' summary: show only the slice each view faces
function hide_edge_on_slices(fig)
% Each of the three views shows one MRI slice face-on; hide the two seen
% edge-on and use exact axis views with an orthographic projection. (The bar
% above the coronal view was not one of these: see template_mri_without_slab.)
for ax = findall(fig,'Type','axes')'
    hs = findall(ax,'Type','surface','Tag','img');
    if isempty(hs), continue; end
    v = get(ax,'CameraPosition') - get(ax,'CameraTarget');
    [~, d] = max(abs(v));                        % coordinate the view looks along
    e = zeros(1,3); e(d) = sign(v(d)); view(ax, e);
    set(ax,'Projection','orthographic');
    w = zeros(numel(hs),1);
    for i = 1:numel(hs)
        c = get(hs(i), {'XData','YData','ZData'}); c = c{d};
        w(i) = max(c(:)) - min(c(:));             % slice extent along the view
    end
    [~, keep] = min(w);
    set(hs(setdiff(1:numel(hs), keep)), 'Visible','off');
end
end


%% EEGLAB template MRI without the stray tissue above the scalp
function mri = template_mri_without_slab(dipfit)
% The top slices of EEGLAB's template MRI (dipfit standard_mri.mat) hold stray
% tissue 5-8 mm above the scalp, which dipplot draws as a bar over the coronal
% view and blobs over the sagittal one. Blank every voxel above 60 mm that lies
% more than 4 mm above the scalp of the BEM head model. Any other case (other
% MRI or head model, rotated voxel grid, read error) returns the file name, so
% dipplot shows the MRI as before.
mri = dipfit.mrifile;
try
    M = load('-mat', dipfit.mrifile);  A = M.mri.anatomy;  T = M.mri.transform;
    if any(any(abs(T(1:3,1:3) - diag(diag(T(1:3,1:3)))) > 1e-9)), return; end
    H = load('-mat', dipfit.hdmfile);  P = [];
    for b = H.vol.bnd(:)'                      % scalp: the outermost surface
        if isfield(b,'pos'), q = b.pos; else, q = b.pnt; end
        if isempty(P) || max(q(:,3)) > max(P(:,3)), P = q; end
    end
    x = T(1,1)*(1:size(A,1)) + T(1,4);         % voxel centres in mm (1-based)
    y = T(2,2)*(1:size(A,2)) + T(2,4);
    z = T(3,3)*(1:size(A,3)) + T(3,4);
    Q = sortrows(P(P(:,3) > 0, :), -3);       % upper half; highest point first
    [~, iu] = unique(Q(:,1:2), 'rows', 'stable');  Q = Q(iu,:);   % one point per column
    [X, Y] = ndgrid(x, y);
    ztop = griddata(Q(:,1), Q(:,2), Q(:,3), X, Y, 'linear');   % scalp height per column
    ztop(isnan(ztop)) = -Inf;                  % outside the scalp footprint
    A(bsxfun(@gt, reshape(z,1,1,[]), max(ztop + 4, 60))) = 0;
    mri = M.mri;  mri.anatomy = A;
catch
end
end


%% Zero border of one cell around a matrix (padarray needs a toolbox)
function B = padarray0(A)
B = zeros(size(A) + 2);  B(2:end-1, 2:end-1) = A;
end


%% Heatmap ticks
function set_heat_ticks(ax, scales, nChan, isICA, chanlocs, fs)
% As many tick labels as fit the axes at font size fs (points): at most 20
% signal labels, and round scale values with the first scale added.
drawnow;
pp  = getpixelposition(ax);
ptx = fs*96/72;                                   % font height in pixels
if isICA, ylab = arrayfun(@(x) sprintf('IC%d',x), 1:nChan, 'UniformOutput', false);
else,     ylab = {chanlocs.labels};
end
nY = max(2, min([20, nChan, floor(pp(4)/(1.6*ptx))]));
yt = unique(round(linspace(1, nChan, nY)));
set(ax,'YTick',yt,'YTickLabel',ylab(yt));

nS = numel(scales);
if iscell(scales) || nS < 2 || any(diff(double(scales(:))) <= 0)
    xl = scales; if ~iscell(xl), xl = arrayfun(@num2str, xl, 'UniformOutput', false); end
    xt = 1:nS; if nS > 30, xt = round(linspace(1,nS,20)); end
    set(ax,'XTick',xt,'XTickLabel',xl(xt),'XTickLabelRotation',45);
    return
end
sv   = double(scales(:)');
maxX = max(3, floor(pp(3)/(2.5*ptx)));
for step = [1 2 5 10 20 25 50 100 200 500 1000]
    tv = step*ceil(sv(1)/step) : step : sv(end);
    if numel(tv) <= maxX, break; end
end
if isempty(tv) || tv(1) - sv(1) >= step/2, tv = [sv(1) tv]; end
set(ax,'XTick',interp1(sv, 1:nS, tv),'XTickLabel',arrayfun(@num2str, tv, 'UniformOutput', false), ...
    'XTickLabelRotation',0);
end


%% Topography for one cluster panel
function plot_topo(ax, vals, chanlocs, ems)
fin = vals(isfinite(vals));
if isempty(fin)
    axis(ax,'off');
    text(0.5,0.5,'No finite values','Parent',ax,...
         'HorizontalAlignment','center','VerticalAlignment','middle');
    return
end
try
    axes(ax); %#ok<LAXES>
    topoplot(vals, chanlocs, 'emarker', {'.','k',ems,1}, 'electrodes','on');
    if min(fin) < max(fin), clim(ax,[min(fin) max(fin)]); end
    colormap(ax,'parula');
catch
    bar(ax, vals); box(ax,'on'); ylabel(ax,'Entropy');
end
end


%% Local helpers
function out = ifelse(cond,a,b)
if cond, out=a; else, out=b; end
end
