function ascent_plot(entropyData, chanlocs, entropyType, scales, varargin)
% ascent_plot  Visualize uni/multi/time-resolved entropy with robust NaN handling.
%
% entropyData : [chan x scale] or [chan x scale x time]
% chanlocs    : EEGLAB channel locations
% entropyType : label, e.g., 'RCMFEsigma' or 'MSE (std)'
% scales      : cellstr of scale labels or numeric 1:S
% varargin    : context-dependent optional inputs:
%   ICA name-value pairs (can appear anywhere in varargin):
%       'ICA',          true/false
%       'icawinv',      EEG.icawinv  [nChans x nICs]
%   Cluster:
%       'ClusterThresh', scalar in (0,1)  percentile threshold (default 0.75)
%   Aperiodic (entropyType contains 'aperiodic'):
%       varargin{1} = offset        [chan x 1] or [chan x nTimes]
%       varargin{2} = freqs         [1 x freq]
%       varargin{3} = psd           [chan x freq] or [chan x freq x nTimes]
%       varargin{4} = psd_corrected [chan x freq] or [chan x freq x nTimes] or []
%       varargin{5} = times         [1 x nTimes] (time-varying only)
%   Time-resolved:
%       varargin{1} = time_sec      [1 x time]

% Strip named pairs before positional parsing
isICA         = false;
icawinv       = [];
clusterThresh = 0.75;
dipfit        = [];
rmIdx         = [];
for vi = 1:(numel(varargin)-1)
    if ischar(varargin{vi}) || isstring(varargin{vi})
        switch lower(char(varargin{vi}))
            case 'ica'
                isICA = logical(varargin{vi+1});
                rmIdx = [rmIdx vi vi+1];
            case 'icawinv'
                icawinv = varargin{vi+1};
                rmIdx   = [rmIdx vi vi+1];
            case 'clusterthresh'
                clusterThresh = varargin{vi+1};
                rmIdx         = [rmIdx vi vi+1];
            case 'dipfit'
                dipfit = varargin{vi+1};
                rmIdx  = [rmIdx vi vi+1];
        end
    end
end
varargin(rmIdx) = [];

if isICA
    % icawinv is needed for back-projection (non-aperiodic and static aperiodic).
    % Time-varying aperiodic only uses isICA to relabel axes, no back-projection.
    isTimeVaryingAperiodic = contains(lower(entropyType),'aperiodic') && size(entropyData,2) > 1;
    if isempty(icawinv) && ~isTimeVaryingAperiodic
        error('ascent_plot: ICA=true requires ''icawinv'' (EEG.icawinv).');
    end
    if ~isempty(icawinv) && size(icawinv,2) ~= size(entropyData,1)
        error('ascent_plot: icawinv has %d ICs but entropyData has %d rows.', ...
              size(icawinv,2), size(entropyData,1));
    end
end

% Positional varargin
time_sec=[];  offset_in=[];  freqs_in=[];
psd_in=[];    psd_corrected_in=[];  times_in=[];

isAperiodic = contains(lower(entropyType),'aperiodic');
if isAperiodic
    if numel(varargin)>=1, offset_in        = varargin{1}; end
    if numel(varargin)>=2, freqs_in         = varargin{2}; end
    if numel(varargin)>=3, psd_in           = varargin{3}; end
    if numel(varargin)>=4, psd_corrected_in = varargin{4}; end
    if numel(varargin)>=5, times_in         = varargin{5}; end
else
    if ~isempty(varargin), time_sec = varargin{1}; end
end

% =========================================================
% APERIODIC BRANCH
% =========================================================
if isAperiodic
    exponent_mat  = entropyData;
    isTimeVarying = size(exponent_mat,2) > 1;
    if isTimeVarying
        plot_aperiodic_timecourse(exponent_mat,offset_in,freqs_in,psd_in,psd_corrected_in,times_in,isICA);
        return
    end
    exponent_in = exponent_mat(:);

    clrBlue=[0.18 0.45 0.87]; clrOrang=[0.90 0.38 0.15]; clrReg=[0.65 0.08 0.08];
    fAlpha=0.18; lw=2; fsAx=12; fsTtl=12; bgCol=[0.97 0.97 0.97];

    hFig=figure('Color','w','Position',[100 50 820 980],...
                'Name','Aperiodic visualization','NumberTitle','off',...
                'Toolbar','none','Menu','none');
    try icadefs; set(hFig,'color',BACKCOLOR); catch; end

    ax1=subplot(4,2,1);
    finE=exponent_in(isfinite(exponent_in));
    if ~isempty(finE)
        if isICA
            bp_exp=icawinv*exponent_in(:);
            axes(ax1); topoplot(bp_exp,chanlocs,'emarker',{'.','k',7,1},'electrodes','on');
            finBP=bp_exp(isfinite(bp_exp));
            if min(finBP)<max(finBP), clim(ax1,[min(finBP) max(finBP)]); end
        else
            axes(ax1); topoplot(exponent_in,chanlocs,'emarker',{'.','k',7,1},'electrodes','on');
            if min(finE)<max(finE), clim(ax1,[min(finE) max(finE)]); end
        end
        c1=colorbar(ax1,'Location','eastoutside');
        c1.Label.String='Exponent'; c1.Label.FontSize=fsAx-1; c1.TickDirection='out';
    end
    title(ax1,ifelse(isICA,'Exponent (back-projected)','Aperiodic Exponent'),...
          'FontSize',fsTtl,'FontWeight','bold');

    ax2=subplot(4,2,2);
    if ~isempty(offset_in)
        finO=offset_in(isfinite(offset_in));
        if ~isempty(finO)
            if isICA
                bp_off=icawinv*offset_in(:);
                axes(ax2); topoplot(bp_off,chanlocs,'emarker',{'.','k',7,1},'electrodes','on');
                finBO=bp_off(isfinite(bp_off));
                if min(finBO)<max(finBO), clim(ax2,[min(finBO) max(finBO)]); end
            else
                axes(ax2); topoplot(offset_in,chanlocs,'emarker',{'.','k',7,1},'electrodes','on');
                if min(finO)<max(finO), clim(ax2,[min(finO) max(finO)]); end
            end
            c2=colorbar(ax2,'Location','eastoutside');
            c2.Label.String='Offset'; c2.Label.FontSize=fsAx-1; c2.TickDirection='out';
        end
    end
    title(ax2,ifelse(isICA,'Offset (back-projected)','Aperiodic Offset'),...
          'FontSize',fsTtl,'FontWeight','bold');

    ax4=subplot(4,2,[3 4]); hold(ax4,'on');
    if ~isempty(offset_in)
        mask=isfinite(exponent_in)&isfinite(offset_in);
        scatter(ax4,exponent_in(mask),offset_in(mask),50,clrBlue,'filled',...
                'MarkerFaceAlpha',0.65,'MarkerEdgeColor','none');
        if isICA
            [~,ic_max]=max(exponent_in);
            scatter(ax4,exponent_in(ic_max),offset_in(ic_max),100,'r',...
                    'filled','MarkerEdgeColor','k','LineWidth',1.2);
            text(ax4,exponent_in(ic_max),offset_in(ic_max),...
                 sprintf('  IC%d',ic_max),'FontSize',fsAx-1,'Color','r');
        end
        if sum(mask)>2
            cc=corrcoef(exponent_in(mask),offset_in(mask)); r=cc(1,2);
            lm=polyfit(exponent_in(mask),offset_in(mask),1);
            xr=linspace(min(exponent_in(mask)),max(exponent_in(mask)),60);
            plot(ax4,xr,polyval(lm,xr),'-','Color',clrReg,'LineWidth',1.8);
            title(ax4,sprintf('Exponent vs Offset  (r = %.2f)',r),'FontSize',fsTtl,'FontWeight','bold');
        else
            title(ax4,'Exponent vs Offset','FontSize',fsTtl,'FontWeight','bold');
        end
        xlabel(ax4,'Exponent','FontSize',fsAx); ylabel(ax4,'Offset','FontSize',fsAx);
        box(ax4,'on'); ax4.TickDir='out'; set(ax4,'Color',bgCol,'LineWidth',0.8);
    end

    if ~isempty(psd_in) && ~isempty(freqs_in)
        log_psd=log10(psd_in); mu_raw=mean(log_psd,1,'omitnan'); log_f=log10(freqs_in);
        ax3=subplot(4,2,[5 6]); hold(ax3,'on');
        fill(ax3,[freqs_in fliplr(freqs_in)],...
             [mu_raw+std(log_psd,0,1,'omitnan') fliplr(mu_raw-std(log_psd,0,1,'omitnan'))],...
             clrBlue,'FaceAlpha',fAlpha,'EdgeColor','none');
        hRaw  = plot(ax3,freqs_in,mu_raw,'-','Color',clrBlue,'LineWidth',lw);
        ap_mean=mean(offset_in,'omitnan')-mean(exponent_in,'omitnan').*log_f;
        hSlope=plot(ax3,freqs_in,ap_mean,'--','Color',clrBlue,'LineWidth',lw-0.4);
        if ~isempty(psd_corrected_in)
            log_pc=log10(psd_corrected_in); mu_c=mean(log_pc,1,'omitnan'); sd_c=std(log_pc,0,1,'omitnan');
            fill(ax3,[freqs_in fliplr(freqs_in)],[mu_c+sd_c fliplr(mu_c-sd_c)],...
                 clrOrang,'FaceAlpha',fAlpha,'EdgeColor','none');
            hCor  = plot(ax3,freqs_in,mu_c,'-','Color',clrOrang,'LineWidth',lw);
            hSlopeC= plot(ax3,freqs_in,zeros(size(freqs_in)),'--','Color',clrOrang,'LineWidth',lw-0.4);
            legend(ax3, [hRaw hSlope hCor hSlopeC], ...
                   {'Raw PSD', '1/f slope (raw)', 'Corrected PSD', '1/f slope (corrected)'}, ...
                   'Location', 'eastoutside', 'Box', 'off', 'FontSize', fsAx-1);
            title(ax3,'PSD: raw vs corrected  (mean +/- SD)','FontSize',fsTtl,'FontWeight','bold');
        else
            legend(ax3,[hRaw hSlope],{'Raw PSD','1/f slope'},...
                   'Location','northeast','Box','off','FontSize',fsAx-1);
            title(ax3,'Raw PSD + aperiodic fit  (mean +/- SD)','FontSize',fsTtl,'FontWeight','bold');
        end
        xlabel(ax3,'Frequency (Hz)','FontSize',fsAx);
        ylabel(ax3,'Power (log_{10} \muV^2/Hz)','FontSize',fsAx);
        xlim(ax3,[freqs_in(1) freqs_in(end)]); box(ax3,'on'); ax3.TickDir='out';
        set(ax3,'Color',bgCol,'LineWidth',0.8);
    end

    ax_pr=subplot(4,2,7); ax_pc=subplot(4,2,8);
    if ~isempty(psd_in) && ~isempty(freqs_in)
        aMask=freqs_in>=8 & freqs_in<=13;
        if any(aMask)
            alpha_raw=mean(log10(psd_in(:,aMask)),2,'omitnan');
            if isICA
                bp_alpha=icawinv*alpha_raw(:); finP=bp_alpha(isfinite(bp_alpha));
                if ~isempty(finP)
                    axes(ax_pr); topoplot(bp_alpha,chanlocs,'emarker',{'.','k',7,1},'electrodes','on');
                    if min(finP)<max(finP), clim(ax_pr,[min(finP) max(finP)]); end
                    c3=colorbar(ax_pr,'Location','eastoutside');
                    c3.Label.String='log_{10} \muV^2/Hz'; c3.Label.FontSize=fsAx-1; c3.TickDirection='out';
                end
                title(ax_pr,'Alpha power (back-projected)','FontSize',fsTtl,'FontWeight','bold');
                [~,ic_max]=max(exponent_in);
                axes(ax_pc); topoplot(icawinv(:,ic_max),chanlocs,'emarker',{'.','k',7,1},'electrodes','on');
                finIC=icawinv(:,ic_max); finIC=finIC(isfinite(finIC));
                if min(finIC)<max(finIC), clim(ax_pc,[min(finIC) max(finIC)]); end
                c4=colorbar(ax_pc,'Location','eastoutside');
                c4.Label.String='Weight'; c4.Label.FontSize=fsAx-1; c4.TickDirection='out';
                title(ax_pc,sprintf('IC%d scalp map (max exponent)',ic_max),'FontSize',fsTtl,'FontWeight','bold');
            else
                finP=alpha_raw(isfinite(alpha_raw));
                if ~isempty(finP)
                    axes(ax_pr); topoplot(alpha_raw,chanlocs,'emarker',{'.','k',7,1},'electrodes','on');
                    if min(finP)<max(finP), clim(ax_pr,[min(finP) max(finP)]); end
                    c3=colorbar(ax_pr,'Location','eastoutside');
                    c3.Label.String='log_{10} \muV^2/Hz'; c3.Label.FontSize=fsAx-1; c3.TickDirection='out';
                end
                title(ax_pr,'Alpha power (8-13 Hz) - raw','FontSize',fsTtl,'FontWeight','bold');
                if ~isempty(psd_corrected_in)
                    alpha_cor=mean(log10(psd_corrected_in(:,aMask)),2,'omitnan');
                    finPC=alpha_cor(isfinite(alpha_cor));
                    if ~isempty(finPC)
                        axes(ax_pc); topoplot(alpha_cor,chanlocs,'emarker',{'.','k',7,1},'electrodes','on');
                        if min(finPC)<max(finPC), clim(ax_pc,[min(finPC) max(finPC)]); end
                        c4=colorbar(ax_pc,'Location','eastoutside');
                        c4.Label.String='log_{10} \muV^2/Hz'; c4.Label.FontSize=fsAx-1; c4.TickDirection='out';
                    end
                    title(ax_pc,'Alpha power (8-13 Hz) - corrected','FontSize',fsTtl,'FontWeight','bold');
                else
                    axes(ax_pc); axis off;
                    text(0.5,0.5,'Correction not applied','HorizontalAlignment','center',...
                         'VerticalAlignment','middle','FontSize',fsAx,'Color',[0.55 0.55 0.55]);
                end
            end
        else
            warning('ascent_plot: alpha band not covered by freqs_in.');
            axes(ax_pr); axis off; axes(ax_pc); axis off;
        end
    else
        axes(ax_pr); axis off; axes(ax_pc); axis off;
    end

    colormap(hFig,parula);
    allAx=findall(hFig,'type','axes');
    for k=1:numel(allAx)
        pos=allAx(k).Position;
        allAx(k).Position=[pos(1),pos(2)+0.02,pos(3),pos(4)-0.04];
    end
    colormap(hFig,parula);
    set(findall(hFig,'type','axes'),'FontSize',fsAx,'FontWeight','bold');
    return
end


% =========================================================
% ALL OTHER MEASURES
% =========================================================

isTimeResolved = ndims(entropyData)==3;
if isTimeResolved
    entropyData3D = entropyData;
    entropyData   = mean(entropyData3D,3,'omitnan');
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
    try icadefs; set(hFig2,'color',BACKCOLOR); catch; end

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
    try icadefs; set(hFig3,'color',BACKCOLOR); catch; end
    vals=entropyData(:); finiteVals=vals(isfinite(vals));
    ax_uni=axes(hFig3);
    if isempty(finiteVals)
        axis(ax_uni,'off');
        text(0.5,0.5,'No finite values to plot','Parent',ax_uni,...
             'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold');
    else
        if isICA
            bp_vals=icawinv*vals(:); finBP=bp_vals(isfinite(bp_vals));
            axes(ax_uni); topoplot(bp_vals,chanlocs,'emarker',{'.','k',15,1},'electrodes','labels');
            if min(finBP)<max(finBP), clim(ax_uni,[min(finBP)*0.95 max(finBP)*1.05]); end
            title(ax_uni,[entropyType ' (back-projected)'],'Interpreter','none');
        else
            axes(ax_uni); topoplot(vals,chanlocs,'emarker',{'.','k',15,1},'electrodes','labels');
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
    try icadefs; set(hFig4,'color',BACKCOLOR); catch; end
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
    try icadefs; set(hFig5,'color',BACKCOLOR); catch; end
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



%% Aperiodic timecourse plot (SPRiNT-style)
function plot_aperiodic_timecourse(exponent_t, offset_t, freqs, psd_t, psd_corr_t, times, isICA)
% Layout (single figure):
%   Row 1       : Raw PSD spectrogram heatmap     (freq x time, mean across channels/ICs)
%   Row 2       : Corrected PSD spectrogram heatmap (if available)
%   Row 3       : Exponent(t) +/- SD [left axis] | Offset(t) +/- SD [right axis]
%   Row 4       : Band-averaged PSD time courses - theta / alpha / beta overlaid
if nargin < 7, isICA = false; end
if isICA
    unitLbl = 'IC(s)';
    unitPlu = 'ICs';
else
    unitLbl = 'channel(s)';
    unitPlu = 'channels';
end

[nChan, ~, nTimes] = size(psd_t);
hasCorrected = ~isempty(psd_corr_t) && size(psd_corr_t, 3) == nTimes;

if isempty(times) || numel(times) ~= nTimes
    times  = 1:nTimes;
    xLabel = 'Window';
else
    times  = times(:)';
    xLabel = 'Time (s)';
end

% Mean log PSD across channels for heatmap
log_psd_mean  = squeeze(mean(log10(psd_t + eps), 1, 'omitnan'));   % [nFreqs x nTimes]
if hasCorrected
    log_corr_mean = squeeze(mean(log10(psd_corr_t + eps), 1, 'omitnan'));
end

% Aperiodic parameter statistics across channels
exp_mean = mean(exponent_t, 1, 'omitnan');
exp_sd   = std(exponent_t,  0, 1, 'omitnan');
off_mean = mean(offset_t,   1, 'omitnan');
off_sd   = std(offset_t,    0, 1, 'omitnan');

clrExp   = [0.18 0.45 0.87];
clrOff   = [0.40 0.40 0.40];
fAlpha   = 0.15;

nRows = 3 + hasCorrected;   % 3 always + 1 optional corrected heatmap
figH  = 160 + nRows * 170;

hFig = figure('Color', 'w', 'Position', [80 50 920 figH], ...
    'Name', 'Aperiodic time course (SPRiNT-style)', ...
    'NumberTitle', 'off', 'Toolbar', 'none', 'Menu', 'none');
try icadefs; set(hFig, 'color', BACKCOLOR); catch; end %#ok<NODEF>

row = 0;

% Row 1: Raw PSD spectrogram
row  = row + 1;
ax1  = subplot(nRows, 1, row);
imagesc(ax1, times, freqs, log_psd_mean);
axis(ax1, 'xy');                            % frequency increases upward
set(ax1, 'TickDir', 'out');
ylabel(ax1, 'Frequency (Hz)');
title(ax1, sprintf('Raw PSD  (mean across %d %s)', nChan, unitLbl), 'FontWeight', 'bold');
cb1 = colorbar(ax1, 'Location', 'eastoutside');
ylabel(cb1, 'log_{10} \muV^2/Hz', 'FontSize', 8);
colormap(ax1, 'parula');
xlim(ax1, [times(1) times(end)]);

% Row 2: Corrected PSD spectrogram (optional)
if hasCorrected
    row  = row + 1;
    ax2  = subplot(nRows, 1, row);
    imagesc(ax2, times, freqs, log_corr_mean);
    axis(ax2, 'xy');
    set(ax2, 'TickDir', 'out');
    ylabel(ax2, 'Frequency (Hz)');
    title(ax2, 'Corrected PSD  (aperiodic component removed)', 'FontWeight', 'bold');
    cb2 = colorbar(ax2, 'Location', 'eastoutside');
    ylabel(cb2, 'log_{10} \muV^2/Hz', 'FontSize', 8);
    colormap(ax2, 'parula');
    xlim(ax2, [times(1) times(end)]);
end

% Row 3: Exponent + Offset time series
row  = row + 1;
ax3  = subplot(nRows, 1, row);
hold(ax3, 'on');

yyaxis(ax3, 'left');
fill(ax3, [times fliplr(times)], ...
     [exp_mean + exp_sd  fliplr(exp_mean - exp_sd)], ...
     clrExp, 'FaceAlpha', fAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
hExp = plot(ax3, times, exp_mean, '-', 'Color', clrExp, 'LineWidth', 2, 'DisplayName', 'Exponent');
ylabel(ax3, 'Exponent');
ax3.YColor = clrExp;

yyaxis(ax3, 'right');
fill(ax3, [times fliplr(times)], ...
     [off_mean + off_sd  fliplr(off_mean - off_sd)], ...
     clrOff, 'FaceAlpha', fAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
hOff = plot(ax3, times, off_mean, '--', 'Color', clrOff, 'LineWidth', 2, 'DisplayName', 'Offset');
ylabel(ax3, 'Offset');
ax3.YColor = clrOff;

legend(ax3, [hExp hOff], 'Location', 'eastoutside', 'Box', 'off', 'FontSize', 10);
title(ax3, sprintf('Aperiodic parameters  (mean \\pm SD across %s)', unitPlu), 'FontWeight', 'bold');
xlim(ax3, [times(1) times(end)]);
set(ax3, 'TickDir', 'out'); box(ax3, 'on');

% Row 4: Simplified band time courses (theta / alpha / beta)
row  = row + 1;
ax4  = subplot(nRows, 1, row);
hold(ax4, 'on');

bandNames  = {'Theta', 'Alpha', 'Beta'};
bandEdges  = [4 8; 8 13; 13 30];
clrBands   = {[0.47 0.67 0.19], [0.85 0.33 0.10], [0.49 0.18 0.56]};
src        = psd_t;
srcLabel   = 'Raw';
if hasCorrected
    src      = psd_corr_t;
    srcLabel = 'Corrected';
end

hLinesRaw  = gobjects(1, 3);
hLines     = gobjects(1, 3);
for ib = 1:3
    fMask = freqs >= bandEdges(ib, 1) & freqs <= bandEdges(ib, 2);
    if ~any(fMask), continue; end

    % --- Raw (plotted first, underneath) ---
    band_raw_avg = squeeze(mean(log10(psd_t(:, fMask, :) + eps), 2, 'omitnan'));
    if nChan == 1, band_raw_avg = band_raw_avg(:)'; end
    raw_mean = mean(band_raw_avg, 1, 'omitnan');
    raw_sd   = std(band_raw_avg,  0, 1, 'omitnan');
    fill(ax4, [times fliplr(times)], ...
         [raw_mean + raw_sd  fliplr(raw_mean - raw_sd)], ...
         clrBands{ib} * 0.55 + 0.45, 'FaceAlpha', fAlpha, ...
         'EdgeColor', 'none', 'HandleVisibility', 'off');
    hLinesRaw(ib) = plot(ax4, times, raw_mean, '--', ...
        'Color', clrBands{ib} * 0.55 + 0.45, ...
        'LineWidth', 1.2, 'DisplayName', bandNames{ib});

    % --- Corrected (plotted on top) ---
    if hasCorrected
        band_freq_avg = squeeze(mean(log10(src(:, fMask, :) + eps), 2, 'omitnan'));
        if nChan == 1, band_freq_avg = band_freq_avg(:)'; end
        band_mean = mean(band_freq_avg, 1, 'omitnan');
        band_sd   = std(band_freq_avg,  0, 1, 'omitnan');
        fill(ax4, [times fliplr(times)], ...
             [band_mean + band_sd  fliplr(band_mean - band_sd)], ...
             clrBands{ib}, 'FaceAlpha', fAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        hLines(ib) = plot(ax4, times, band_mean, '-', 'Color', clrBands{ib}, ...
            'LineWidth', 1.8, 'DisplayName', [bandNames{ib} ' (corrected)']);
    end
end

% Legend: raw first, then corrected - placed outside to the right
validRaw  = isgraphics(hLinesRaw);
validCorr = isgraphics(hLines);
allH = [hLinesRaw(validRaw) hLines(validCorr)];
if ~isempty(allH)
    leg = legend(ax4, allH, 'Location', 'eastoutside', 'Box', 'off', 'FontSize', 10);
end

xlabel(ax4, xLabel);
ylabel(ax4, 'log_{10} \muV^2/Hz');
title(ax4, sprintf('%s PSD band averages  (mean \\pm SD)', srcLabel), 'FontWeight', 'bold');
xlim(ax4, [times(1) times(end)]);
set(ax4, 'TickDir', 'out'); box(ax4, 'on');

% Final styling
set(findall(hFig, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');
colormap(hFig, 'parula');
try
    sgtitle(hFig, 'Aperiodic time course', 'FontSize', 14, 'FontWeight', 'bold');
catch; end

% Align all plot axes to the same left edge and width so every row lines up.
% Strategy: take the tightest common bounds across all axes (max left edge,
% min right edge), then apply uniformly. Colorbars are repositioned manually.
drawnow;
if hasCorrected
    plotAxes = [ax1 ax2 ax3 ax4];
    cbAxes   = [ax1 ax2];
    cbList   = [cb1 cb2];
else
    plotAxes = [ax1 ax3 ax4];
    cbAxes   = ax1;
    cbList   = cb1;
end
leftEdges  = arrayfun(@(a) a.Position(1),              plotAxes);
rightEdges = arrayfun(@(a) a.Position(1)+a.Position(3), plotAxes);
refLeft  = max(leftEdges);
refWidth = min(rightEdges) - refLeft;
for ax = plotAxes
    p = ax.Position;
    ax.Position = [refLeft, p(2), refWidth, p(4)];
end
cbW   = 0.025;
cbGap = 0.008;
for ii = 1:numel(cbAxes)
    p = cbAxes(ii).Position;
    cbList(ii).Position = [refLeft+refWidth+cbGap, p(2)+p(4)*0.1, cbW, p(4)*0.8];
end
end

%% Dipole info text panel (subplot [12] in single-cluster ICA layout)
function show_dipole_info(ax, dipfit, ic_idx, centroid_ic)
% show_dipole_info  Show MNI coordinates and RV for the centroid IC.
% The full 3D dipole visualization is rendered by dipplot in a separate figure.

axis(ax,'off');

if isempty(dipfit) || ~isfield(dipfit,'model') || isempty(dipfit.model)
    text(0.5, 0.5, sprintf('Pass ''dipfit'',EEG.dipfit\nto show dipole sources'), ...
         'Parent',ax,'HorizontalAlignment','center','VerticalAlignment','middle', ...
         'FontSize',8,'Color',[.5 .5 .5],'Interpreter','none');
    return;
end

if centroid_ic > numel(dipfit.model) || isempty(dipfit.model(centroid_ic).posxyz)
    text(0.5, 0.5, sprintf('IC%d: no dipole fitted', centroid_ic), ...
         'Parent',ax,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8);
    return;
end

m        = dipfit.model(centroid_ic);
pos      = m.posxyz(1,:);
rv_pct   = m.rv * 100;
rvThresh = 0.15;

n_valid = sum(arrayfun(@(i) i<=numel(dipfit.model) && ...
    ~isempty(dipfit.model(i).posxyz) && dipfit.model(i).rv < rvThresh, ic_idx));

lines = { sprintf('Centroid IC%d', centroid_ic), ...
          sprintf('MNI  [%.0f  %.0f  %.0f] mm', pos(1), pos(2), pos(3)), ...
          sprintf('RV = %.1f%%', rv_pct), ...
          '', ...
          sprintf('%d/%d cluster ICs', n_valid, numel(ic_idx)), ...
          sprintf('with RV < %.0f%%', rvThresh*100), ...
          '', ...
          '(See dipplot figure)' };

ypos = linspace(0.88, 0.12, numel(lines));
for ii = 1:numel(lines)
    if ii == 1,     fw = 'bold';   fs = 9;
    elseif ii == 8, fw = 'normal'; fs = 7.5;
    else,           fw = 'normal'; fs = 8.5;
    end
    text(0.5, ypos(ii), lines{ii}, 'Parent', ax, ...
         'HorizontalAlignment','center','VerticalAlignment','middle', ...
         'FontSize',fs,'FontWeight',fw,'Interpreter','none');
end
title(ax, 'Dipole source', 'FontSize', 9, 'FontWeight','bold');
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
        'DisplayName',sprintf('%s (n=%d)',rname,n_r));

    xpos = xpos+1;
    vals = mean_ent(idx);
    scatter(ax2,xpos+(rand(n_r,1)-0.5)*0.3,vals,30,clr,...
            'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
    plot(ax2,xpos+[-0.3 0.3],repmat(mean(vals,'omitnan'),1,2),...
         '-','Color',clr*0.6,'LineWidth',2.5);
    xt_pos(end+1)=xpos; xt_lbl{end+1}=sprintf('%s (n=%d)',rname,n_r);
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