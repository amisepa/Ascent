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
% Copyright (C) - ASCENT EEGLAB PLUGIN - Cedric Cannard, 2021-2025

isICA         = opt.isICA;
icawinv       = opt.icawinv;
clusterThresh = opt.clusterThresh;
dipfit        = opt.dipfit;
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
    ic_idx_1 = [];  pc_1 = 1;   % cluster 1, for the dipole figure

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
        if k == 1, ic_idx_1 = find(chanMask_k); pc_1 = pc_k; end

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
            xlim(ax_bar, [0 nChan+1]);
            xlabel(ax_bar,'IC'); ylabel(ax_bar,'Entropy'); box(ax_bar,'on');
            if usedCluster
                title(ax_bar, sprintf('%d ICs x %d scales, peak %s sc %s', nCh_k, nSc_k, ...
                      sigName, scLabel_k), 'Interpreter','none');
            else
                title(ax_bar, sprintf('Peak %s, sc %s', sigName, scLabel_k), 'Interpreter','none');
            end
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
        plot_regions(entropyData, xvals, dipfit, ic_idx_1, pc_1, entropyType, fs, clrCluster{1});
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


%% Region entropy comparison + dipole sources (ICA + dipfit)
function plot_regions(ent, scales, dipfit, ic_idx, pc_k, entropyType, fs, clu_color, rv_thresh)
if nargin < 9, rv_thresh = 0.15; end
lw = max(1, fs/9);

n_ics = size(ent, 1);
brain_ics = [];  ic_pos = [];

for ik = 1:n_ics
    if ik > numel(dipfit.model), continue; end
    m = dipfit.model(ik);
    if isempty(m.posxyz) || m.rv > rv_thresh, continue; end
    brain_ics(end+1) = ik;           %#ok<AGROW>
    ic_pos(end+1,:)  = m.posxyz(1,:);%#ok<AGROW>
end

if numel(brain_ics) < 3
    warning('ascent_plot: only %d ICs passed RV<%.0f%% — skipping region plot.', ...
            numel(brain_ics), rv_thresh*100);
    return;
end
fprintf('Region plot: %d ICs (RV<%.0f%%)\n', numel(brain_ics), rv_thresh*100);

% No blue or orange here: the dipole panel uses them for the other and the cluster ICs
region_names = {'Frontal','Temporal','Central','Parietal','Occipital'};
region_clrs  = {[0.13 0.63 0.30],[0.45 0.30 0.85],...
                [0.85 0.60 0.00],[0.85 0.20 0.55],[0.45 0.45 0.45]};

ic_regions = cell(numel(brain_ics),1);
for ii = 1:numel(brain_ics)
    x=ic_pos(ii,1); y=ic_pos(ii,2); z=ic_pos(ii,3);
    if     y >  20,            ic_regions{ii}='Frontal';
    elseif y < -65,            ic_regions{ii}='Occipital';
    elseif abs(x)>35 && z<25, ic_regions{ii}='Temporal';
    elseif abs(y)<20 && z>30, ic_regions{ii}='Central';
    else,                      ic_regions{ii}='Parietal';
    end
end

ent_sel  = ent(brain_ics,:);
mean_ent = mean(ent_sel, 2, 'omitnan');

% ---- Run dipplot into temp figure, capture as image --------------------
dip_img = [];
if ~isempty(dipfit) && isfield(dipfit,'mrifile') && ~isempty(ic_idx)
    try
        % All dipolar ICs of the decomposition (the set the region plots below use),
        % cluster ICs in orange, their centroid in red, the others in blue. Few
        % cluster ICs may be dipolar, so plotting only those can leave 3-4 dots.
        oth_color  = [0.30 0.55 1.00];
        cent_color = [0.85 0.10 0.10];
        sources = struct('posxyz',{},'momxyz',{},'rv',{});
        colors  = {};
        for ik = 1:numel(dipfit.model)
            m = dipfit.model(ik);
            if isempty(m.posxyz) || isempty(m.rv) || m.rv > rv_thresh, continue; end
            sources(end+1).posxyz = m.posxyz(1,:); %#ok<AGROW>
            sources(end).momxyz   = [0 0 0];
            sources(end).rv       = m.rv;
            if ik == pc_k, colors{end+1} = cent_color;          %#ok<AGROW>
            elseif ismember(ik, ic_idx), colors{end+1} = clu_color; %#ok<AGROW>
            else, colors{end+1} = oth_color;                    %#ok<AGROW>
            end
        end
        n_clu_dip = sum(cellfun(@(c) ~isequal(c, oth_color), colors));
        if ~isempty(sources)
            pre_figs = findall(0,'Type','figure');
            dipplot(sources,'mri',template_mri_without_slab(dipfit),'coordformat','MNI','color',colors,...
                'summary','3d','dipolesize',15,'dipolelength',0,'gui','off','mesh','off');
            drawnow; pause(0.3);
            new_figs = setdiff(findall(0,'Type','figure'),pre_figs);
            if ~isempty(new_figs)
                tmp_fig = new_figs(end);
                hide_edge_on_slices(tmp_fig);
                set(tmp_fig,'Color','k','InvertHardcopy','off');
                try
                    dip_img = print(tmp_fig,'-RGBImage','-r300');  % sharper than getframe
                catch
                    dip_img = getframe(tmp_fig).cdata;
                end
                close(tmp_fig);
            end
        end
    catch ME
        warning('ascent_plot: dipplot capture failed — %s', ME.message);
    end
end

% ---- Combined figure layout --------------------------------------------
has_dip = ~isempty(dip_img);
fig = figure('Color','w','InvertHardCopy','off','NumberTitle','off',...
             'Name','Dipole sources + entropy by region');
try icadefs; set(fig,'color',BACKCOLOR); catch; end %#ok<NODEF>
set(fig,'Units','pixels'); p0 = get(fig,'Position');
if has_dip
    sz = [1000 720]; rowH = 0.50;
else
    sz = [1000 420]; rowH = 1;
end
set(fig,'Position',[p0(1) max(40, p0(2)+p0(4)-sz(2)) sz]);
if has_dip
    ax_dip = axes(fig,'Position',[0.01 0.50 0.98 0.43]);
    image(ax_dip, dip_img); axis(ax_dip,'image','off');
    title(ax_dip, sprintf('Dipole sources, RV<%.0f%%: %d ICs (cluster ICs in orange: %d)', ...
          rv_thresh*100, numel(sources), n_clu_dip), 'FontSize',fs,'FontWeight','bold');
end
% OuterPosition: MATLAB fits the tick labels, axis labels and titles inside
ax1 = axes(fig,'OuterPosition',[0    0 0.64 rowH]);
ax2 = axes(fig,'OuterPosition',[0.64 0 0.36 rowH]);
hold(ax1,'on'); box(ax1,'on');
hold(ax2,'on'); box(ax2,'on');

xpos=0; xt_pos=[]; xt_lbl={}; leg_h=gobjects(0); allv=[];
min_n = 3;
for rr = 1:numel(region_names)
    rname = region_names{rr};
    idx   = strcmp(ic_regions, rname);
    n_r   = sum(idx);
    if n_r < min_n, continue; end
    clr = region_clrs{rr};

    mu  = mean(ent_sel(idx,:),1,'omitnan');
    sem = std( ent_sel(idx,:),0,1,'omitnan') / sqrt(n_r);
    fill(ax1,[scales fliplr(scales)],[mu+sem fliplr(mu-sem)],...
         clr,'FaceAlpha',0.15,'EdgeColor','none');
    leg_h(end+1) = plot(ax1,scales,mu,'Color',clr,'LineWidth',2.2*lw,...
        'DisplayName',sprintf('%s (n=%d)',rname,n_r)); %#ok<AGROW>

    xpos = xpos+1;
    vals = mean_ent(idx); allv = [allv; vals]; %#ok<AGROW>
    scatter(ax2,xpos+(rand(n_r,1)-0.5)*0.3,vals,30*lw^2,clr,...
            'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
    plot(ax2,xpos+[-0.3 0.3],repmat(mean(vals,'omitnan'),1,2),...
         '-','Color',clr*0.6,'LineWidth',2.5*lw);
    xt_pos(end+1)=xpos; xt_lbl{end+1}=rname; %#ok<AGROW> (n is in the legend)
end

xlim(ax1,[scales(1) scales(end)]);
xlabel(ax1,'Scale'); ylabel(ax1,entropyType,'Interpreter','none');
title(ax1,'Scale-resolved entropy by region  (mean \pm SEM)');
set(ax1,'TickDir','out');

xlim(ax2,[0 xpos+1]);
allv = allv(isfinite(allv));
if numel(allv) > 1 && max(allv) > min(allv)   % pad so no point sits on the frame
    pad = 0.08*(max(allv)-min(allv)); ylim(ax2,[min(allv)-pad max(allv)+pad]);
end
set(ax2,'XTick',xt_pos,'XTickLabel',xt_lbl,'TickDir','out','XTickLabelRotation',30);
ylabel(ax2,['Mean ' entropyType],'Interpreter','none');
title(ax2,'Mean entropy by region');
set([ax1 ax2],'FontSize',fs,'FontWeight','bold', ...
    'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1);
legend(ax1,leg_h,'Location','best','Box','off','FontSize',max(fs-1,8));
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
