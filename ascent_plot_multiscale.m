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
    % Multiscale + multichannel
    % Layout: 3 rows x 4 cols
    %   Cols 1-2 (all rows): heatmap
    %   Col 3-4, row k:      topo | curve for cluster k
    %   Up to 3 clusters shown; unused rows are blanked.
    % -------------------------------------------------------

    finiteMask=isfinite(entropyData);
    if ~any(finiteMask(:))
        warning('ascent_plot: all values are NaN; nothing to plot.'); return
    end

    % Cluster detection (up to 3)
    [clusterMasks,nClusters,centroids,usedCluster] = ...
        find_entropy_clusters(entropyData, clusterThresh, 3);
    nShow = min(nClusters,3);

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

    nScales=size(entropyData,2);

    % Heatmap: cols 1-2 across all 3 rows -> subplot indices [1 2 5 6 9 10]
    ax_heat=subplot(3,4,[1 2 5 6 9 10]);
    imagesc(ax_heat,1:nScales,1:nChan,entropyData); axis(ax_heat,'tight');
    set(ax_heat,'TickDir','out'); box(ax_heat,'on');

    if isICA
        Yticks=arrayfun(@(x){sprintf('IC%d',x)},1:nChan);
    else
        Yticks={chanlocs.labels};
    end
    newY=1:nChan; if nChan>30, newY=round(linspace(1,nChan,20)); end
    set(ax_heat,'YTick',newY,'YTickLabel',Yticks(newY),'FontWeight','normal');

    if iscell(scales)
        Xticks=scales; nX=numel(scales);
    else
        Xticks=arrayfun(@(x){num2str(x)},scales); nX=numel(scales);
    end
    newX=1:nX; if nX>30, newX=round(linspace(1,nX,20)); end
    set(ax_heat,'XTick',newX,'XTickLabel',Xticks(newX),'FontWeight','normal');
    ax_heat.TickLabelInterpreter='none';
    ax_heat.XTickLabelRotation=45;
    ax_heat.PositionConstraint='outerposition';
    drawnow;
    outer=ax_heat.OuterPosition; ti=ax_heat.TightInset;
    newBottom=max(outer(2),ti(2)+0.02);
    newHeight=max(outer(4)-(newBottom-outer(2))-(ti(4)+0.01),0.1);
    ax_heat.OuterPosition=[outer(1),newBottom,outer(3),newHeight];

    colormap(ax_heat,'parula');
    ch=colorbar(ax_heat); ylabel(ch,'Entropy','FontWeight','bold','FontSize',9);
    xlabel(ax_heat,'Scales');
    ylabel(ax_heat,ifelse(isICA,'ICs','EEG channels'));
    title(ax_heat,entropyType,'Interpreter','none');

    % Overlay cluster boundaries using the shared cluster colours.
    % The binary mask is Gaussian-smoothed before contouring to fill small
    % channel-membership gaps and produce a single clean boundary per cluster.
    if usedCluster
        sigma = 1.5;
        hw    = ceil(3 * sigma);
        [gx, gy] = meshgrid(-hw:hw, -hw:hw);
        kern  = exp(-(gx.^2 + gy.^2) / (2 * sigma^2));
        kern  = kern / sum(kern(:));
        hold(ax_heat,'on');
        for k = 1:nShow
            mk_smooth = conv2(double(clusterMasks(:,:,k)), kern, 'same');
            contour(ax_heat, 1:nScales, 1:nChan, mk_smooth, [0.35 0.35], ...
                    'Color', clrCluster{k}, 'LineWidth', 1.8);
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

    % Right-side subplot slot indices in a 3x4 grid
    % Row 1: 3,4 | Row 2: 7,8 | Row 3: 11,12
    topoSlots  = [3  7  11];
    curveSlots = [4  8  12];

    ic_idx_k = [];  pc_k = 1;   % initialised here; set inside each branch

    if isICA && nShow == 1
        % ---------------------------------------------------------------
        % Single ICA cluster: vertical 3-panel layout in cols 3-4
        %   Row 1  [3 4]  : IC bar chart
        %   Row 2  [7 8]  : Entropy curve
        %   Row 3 [11 12] : Back-projected scalp topo
        % ---------------------------------------------------------------
        k = 1;
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
        if usedCluster
            rowLabel = sprintf('Cluster 1  (%d ICs x %d sc, IC%d sc%s)', ...
                               nCh_k, nSc_k, pc_k, scLabel_k);
        else
            rowLabel = sprintf('Peak  IC%d, sc %s', pc_k, scLabel_k);
        end
        if usedCluster
            fprintf('Cluster 1: %d ICs x %d scales  (centroid IC%d, scale %s)\n',...
                    nCh_k, nSc_k, pc_k, scLabel_k);
        else
            fprintf('No cluster found - peak: IC%d, scale %s\n', pc_k, scLabel_k);
        end

        % -- Row 1: IC bar --
        ax_bar = subplot(3,4,[3 4]);
        bar(ax_bar, topo_vals_k, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none');
        hold(ax_bar,'on');
        if usedCluster
            ic_idx = find(chanMask_k);
            bar(ax_bar, ic_idx, topo_vals_k(ic_idx), ...
                'FaceColor', clrCluster{k}, 'EdgeColor', 'none');
        end
        bar(ax_bar, pc_k, topo_vals_k(pc_k), 'FaceColor', 'r', 'EdgeColor', 'none');
        xlabel(ax_bar,'IC'); ylabel(ax_bar,'Entropy');
        title(ax_bar, rowLabel, 'Interpreter','none','FontSize',9);
        box(ax_bar,'on');

        % -- Row 2: Entropy curve --
        ax_curve = subplot(3,4,[7 8]);
        hold(ax_curve,'on'); box(ax_curve,'on');
        if all(~isfinite(curve_k)), curve_k = nan(1,nScales); end
        xvals = ifelse(isnumeric(scales), scales, 1:nScales);
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
        plot(ax_curve, xvals, curve_k, 'LineWidth', 2, 'Color', clrCluster{k});
        xlim(ax_curve,[xvals(1) xvals(end)]);
        xlabel(ax_curve,'Scale');
        curveTitle = ifelse(usedCluster, ...
            sprintf('Cluster 1 mean  (%d ICs)', nCh_k), 'Peak IC curve');
        title(ax_curve, curveTitle, 'Interpreter','none','FontSize',9);

        % -- Row 3: Back-projected scalp topo --
        ax_topo_s = subplot(3,4,[11 12]);
        bp_vals = icawinv * topo_vals_k;
        finBP   = bp_vals(isfinite(bp_vals));
        if ~isempty(finBP)
            try
                axes(ax_topo_s); %#ok<LAXES>
                topoplot(bp_vals, chanlocs, 'emarker', {'.','k',7,1}, 'electrodes','on');
                if min(finBP) < max(finBP), clim(ax_topo_s,[min(finBP) max(finBP)]); end
                colormap(ax_topo_s,'parula');
            catch
                bar(ax_topo_s, bp_vals); box(ax_topo_s,'on');
                ylabel(ax_topo_s,'Weight');
            end
        end
        sc_idx_k = find(scaleMask_k);
        if isnumeric(scales)
            scRangeStr = sprintf('sc %g-%g', scales(sc_idx_k(1)), scales(sc_idx_k(end)));
        else
            scRangeStr = sprintf('sc %s-%s', scales{sc_idx_k(1)}, scales{sc_idx_k(end)});
        end
        title(ax_topo_s, sprintf('Back-projected topo  (cluster mean, %s)', scRangeStr), ...
              'Interpreter','none','FontSize',9);
        ic_idx_k = find(chanMask_k);

    else
        % ---------------------------------------------------------------
        % Standard per-row layout: one row per cluster.
        % For ICA: topo slot = back-projected scalp map, curve slot = curve.
        % For non-ICA: topo slot = channel entropy topo, curve slot = curve.
        % ---------------------------------------------------------------
        for k = 1:3
            ax_topo  = subplot(3,4, topoSlots(k));
            ax_curve = subplot(3,4, curveSlots(k));

            if k > nShow
                axis(ax_topo,'off'); axis(ax_curve,'off'); continue
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

            if usedCluster
                centStr  = ifelse(isICA, sprintf('IC%d sc%s',pc_k,scLabel_k), ...
                                         sprintf('%s sc%s',chanlocs(pc_k).labels,scLabel_k));
                rowLabel = sprintf('Cluster %d  (%d %s x %d sc, %s)', k, nCh_k, ...
                                   ifelse(isICA,'ICs','ch'), nSc_k, centStr);
            else
                rowLabel = ifelse(isICA, ...
                    sprintf('Peak  IC%d, sc %s',pc_k,scLabel_k), ...
                    sprintf('Peak  %s, sc %s',chanlocs(pc_k).labels,scLabel_k));
            end

            % Console summary
            if usedCluster
                if isICA
                    fprintf('Cluster %d: %d ICs x %d scales  (centroid IC%d, scale %s)\n',...
                        k,nCh_k,nSc_k,pc_k,scLabel_k);
                else
                    fprintf('Cluster %d: %d ch x %d scales  (centroid %s, scale %s)\n',...
                        k,nCh_k,nSc_k,chanlocs(pc_k).labels,scLabel_k);
                end
            elseif k==1
                if isICA
                    fprintf('No cluster found - peak: IC%d, scale %s\n',pc_k,scLabel_k);
                else
                    fprintf('No cluster found - peak: ch %s, scale %s\n',...
                        chanlocs(pc_k).labels,scLabel_k);
                end
            end

            % Topo slot
            finV = topo_vals_k(isfinite(topo_vals_k));
            if isICA
                % Back-projected scalp topo
                bp_vals = icawinv * topo_vals_k;
                finBP   = bp_vals(isfinite(bp_vals));
                if ~isempty(finBP)
                    try
                        axes(ax_topo); %#ok<LAXES>
                        topoplot(bp_vals, chanlocs, 'emarker', {'.','k',7,1}, 'electrodes','on');
                        if min(finBP) < max(finBP), clim(ax_topo,[min(finBP) max(finBP)]); end
                        colormap(ax_topo,'parula');
                    catch
                        bar(ax_topo, bp_vals); box(ax_topo,'on');
                        ylabel(ax_topo,'Weight');
                    end
                else
                    axis(ax_topo,'off');
                end
                title(ax_topo, rowLabel, 'Interpreter','none','FontSize',9);
            else
                if isempty(finV)
                    axis(ax_topo,'off');
                    text(0.5,0.5,'No finite values','Parent',ax_topo,...
                         'HorizontalAlignment','center','VerticalAlignment','middle');
                else
                    try
                        axes(ax_topo); %#ok<LAXES>
                        topoplot(topo_vals_k,chanlocs,'emarker',{'.','k',7,1},'electrodes','on');
                        lo=min(finV); hi=max(finV);
                        if isfinite(lo)&&isfinite(hi)&&lo<hi, clim(ax_topo,[lo hi]); end
                        colormap(ax_topo,'parula');
                        title(ax_topo,rowLabel,'Interpreter','none','FontSize',9);
                    catch
                        bar(ax_topo,topo_vals_k); box(ax_topo,'on');
                        title(ax_topo,rowLabel,'Interpreter','none','FontSize',9);
                        ylabel(ax_topo,'Entropy');
                    end
                end
            end

            % Curve slot
            axes(ax_curve); %#ok<LAXES>
            hold(ax_curve,'on'); box(ax_curve,'on');
            if all(~isfinite(curve_k)), curve_k=nan(1,nScales); end
            xvals = ifelse(isnumeric(scales), scales, 1:nScales);

            if doShading && isnumeric(scales) && any(scaleMask_k)
                x_lo = xvals(find(scaleMask_k,1,'first'));
                x_hi = xvals(find(scaleMask_k,1,'last'));
                finC = curve_k(isfinite(curve_k));
                if ~isempty(finC)
                    yl_pre=[min(finC)*0.95 max(finC)*1.05];
                    if yl_pre(1)==yl_pre(2), yl_pre=yl_pre+[-0.1 0.1]; end
                    patch(ax_curve,[x_lo x_hi x_hi x_lo],...
                          [yl_pre(1) yl_pre(1) yl_pre(2) yl_pre(2)],...
                          clrCluster{k},'FaceAlpha',0.15,'EdgeColor','none');
                end
            end

            plot(ax_curve,xvals,curve_k,'LineWidth',2,'Color',clrCluster{k});
            xlim(ax_curve,[xvals(1) xvals(end)]);
            xlabel(ax_curve,'Scale'); ylabel(ax_curve,'Entropy');
            curveTitle = ifelse(usedCluster,...
                sprintf('Cluster %d mean  (%d %s)',k,nCh_k,ifelse(isICA,'ICs','ch')),...
                'Peak channel curve');
            title(ax_curve,curveTitle,'Interpreter','none','FontSize',9);
        end
    end

    colormap(hFig2,parula);
    set(findall(hFig2,'type','axes'),'FontSize',9,'FontWeight','bold');

    % Region entropy + dipole sources figure (ICA + dipfit only)
    if isICA && ~isempty(dipfit)
        xvals_num = ifelse(isnumeric(scales), scales, 1:size(entropyData,2));
        plot_regions(entropyData, xvals_num, dipfit, ic_idx_k, pc_k);
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
        c.Label.String='Entropy'; c.Label.FontSize=11; c.Label.FontWeight='bold';
    end
    set(findall(hFig3,'type','axes'),'FontSize',10,'FontWeight','bold');

elseif multiscale && ~multiChan
    % Single channel/IC curve
    hFig4=figure('Color','w','InvertHardCopy','off',...
                 'Name',entropyType,'Toolbar','none','Menu','none','NumberTitle','Off');
    try icadefs; set(hFig4,'color',BACKCOLOR); catch; end %#ok<NODEF>
    ax_s=axes(hFig4);
    xvals = ifelse(isnumeric(scales), scales, 1:numel(scales));
    plot(ax_s,xvals,entropyData,'LineWidth',2);
    xlabel(ax_s,'Scale'); ylabel(ax_s,'Entropy');
    title(ax_s,entropyType,'Interpreter','none');
    axis(ax_s,'tight');
    set(findall(hFig4,'type','axes'),'FontSize',10,'FontWeight','bold');

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
    set(findall(hFig5,'type','axes'),'FontSize',10,'FontWeight','bold');

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
function plot_regions(rcmfe, scales, dipfit, ic_idx, pc_k, rv_thresh)
if nargin < 6, rv_thresh = 0.15; end

n_ics = size(rcmfe, 1);
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

region_names = {'Frontal','Temporal','Central','Parietal','Occipital'};
region_clrs  = {[0.20 0.55 0.90],[0.10 0.72 0.35],...
                [0.90 0.40 0.10],[0.75 0.20 0.75],[0.50 0.50 0.50]};

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

rcmfe_sel = rcmfe(brain_ics,:);
mean_ent  = mean(rcmfe_sel, 2, 'omitnan');

% ---- Run dipplot into temp figure, capture as image --------------------
dip_img = [];
if ~isempty(dipfit) && isfield(dipfit,'mrifile') && nargin >= 4 && ~isempty(ic_idx)
    try
        dot_color  = [0.30 0.55 1.00];
        cent_color = [0.85 0.10 0.10];
        sources = struct('posxyz',{},'momxyz',{},'rv',{},'color',{});
        for ii = 1:numel(ic_idx)
            ik = ic_idx(ii);
            if ik > numel(dipfit.model), continue; end
            m = dipfit.model(ik);
            if isempty(m.posxyz) || m.rv > rv_thresh, continue; end
            sources(end+1).posxyz = m.posxyz(1,:); %#ok<AGROW>
            sources(end).momxyz   = [0 0 0];
            sources(end).rv       = m.rv;
            sources(end).color    = ifelse(ik==pc_k, cent_color, dot_color);
        end
        if ~isempty(sources)
            pre_figs = findall(0,'Type','figure');
            dipplot(sources,'mri',dipfit.mrifile,'coordformat','MNI',...
                'summary','3d','dipolesize',12,'dipolelength',0,'gui','off','mesh','off');
            drawnow; pause(0.3);
            new_figs = setdiff(findall(0,'Type','figure'),pre_figs);
            if ~isempty(new_figs)
                tmp_fig  = new_figs(end);
                dip_img  = getframe(tmp_fig).cdata;
                close(tmp_fig);
            end
        end
    catch ME
        warning('ascent_plot: dipplot capture failed — %s', ME.message);
    end
end

% ---- Combined figure layout --------------------------------------------
has_dip = ~isempty(dip_img);
fig = figure('Color','w','NumberTitle','off',...
             'Name','Dipole sources + RCMFE by region');
if has_dip
    fig.Position = [60 60 1100 700];
    ax_dip = axes(fig,'Position',[0.03 0.50 0.94 0.46]);
    image(ax_dip, dip_img); axis(ax_dip,'off');   % no 'image' — fills panel without clipping
    title(ax_dip, sprintf('Dipole sources  (cluster ICs, RV<%.0f%%)', rv_thresh*100),...
          'FontSize',10,'FontWeight','bold');
    bot = 0.13; ht = 0.33;
else
    fig.Position = [60 60 1100 420];
    bot = 0.18; ht = 0.72;
end
ax1 = axes(fig,'Position',[0.07 bot 0.57 ht]);
ax2 = axes(fig,'Position',[0.71 bot 0.26 ht]);
hold(ax1,'on'); box(ax1,'on');
hold(ax2,'on'); box(ax2,'on');

xpos=0; xt_pos=[]; xt_lbl={}; leg_h=gobjects(0);
min_n = 3;
for rr = 1:numel(region_names)
    rname = region_names{rr};
    idx   = strcmp(ic_regions, rname);
    n_r   = sum(idx);
    if n_r < min_n, continue; end
    clr = region_clrs{rr};

    mu  = mean(rcmfe_sel(idx,:),1,'omitnan');
    sem = std( rcmfe_sel(idx,:),0,1,'omitnan') / sqrt(n_r);
    fill(ax1,[scales fliplr(scales)],[mu+sem fliplr(mu-sem)],...
         clr,'FaceAlpha',0.15,'EdgeColor','none');
    leg_h(end+1) = plot(ax1,scales,mu,'Color',clr,'LineWidth',2.2,...
        'DisplayName',sprintf('%s (n=%d)',rname,n_r)); %#ok<AGROW>

    xpos = xpos+1;
    vals = mean_ent(idx);
    scatter(ax2,xpos+(rand(n_r,1)-0.5)*0.3,vals,30,clr,...
            'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
    plot(ax2,xpos+[-0.3 0.3],repmat(mean(vals,'omitnan'),1,2),...
         '-','Color',clr*0.6,'LineWidth',2.5);
    xt_pos(end+1)=xpos; xt_lbl{end+1}=sprintf('%s (n=%d)',rname,n_r); %#ok<AGROW>
end

xlabel(ax1,'Scale'); ylabel(ax1,'RCMFE');
title(ax1,'Scale-resolved entropy by region  (mean \pm SEM)');
legend(ax1,leg_h,'Location','best','Box','off','FontSize',9);
set(ax1,'TickDir','out');

xlim(ax2,[0 xpos+1]);
set(ax2,'XTick',xt_pos,'XTickLabel',xt_lbl,'TickDir','out','XTickLabelRotation',30);
ylabel(ax2,'Mean RCMFE');
title(ax2,'Mean entropy by region');
set(findall(fig,'type','axes'),'FontSize',10,'FontWeight','bold');
end


%% Local helpers
function out = ifelse(cond,a,b)
if cond, out=a; else, out=b; end
end
