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

%% --- Strip named pairs before positional parsing ---
isICA         = false;
icawinv       = [];
clusterThresh = 0.50;
rmIdx         = [];
for vi = 1:2:(numel(varargin)-1)
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
        end
    end
end
varargin(rmIdx) = [];

if isICA
    if isempty(icawinv)
        error('ascent_plot: ICA=true requires ''icawinv'' (EEG.icawinv).');
    end
    if size(icawinv,2) ~= size(entropyData,1)
        error('ascent_plot: icawinv has %d ICs but entropyData has %d rows.', ...
              size(icawinv,2), size(entropyData,1));
    end
end

%% --- Positional varargin ---
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

%% =========================================================
%% APERIODIC BRANCH
%% =========================================================
if isAperiodic
    exponent_mat  = entropyData;
    isTimeVarying = size(exponent_mat,2) > 1;
    if isTimeVarying
        plot_aperiodic_timecourse(exponent_mat,offset_in,freqs_in,psd_in,psd_corrected_in,times_in);
        return
    end
    exponent_in = exponent_mat(:);

    clrBlue=[0.18 0.45 0.87]; clrOrang=[0.90 0.38 0.15]; clrReg=[0.65 0.08 0.08];
    fAlpha=0.18; lw=2; fsAx=11; fsTtl=12; bgCol=[0.97 0.97 0.97];

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
        hRaw  =plot(ax3,freqs_in,mu_raw,'-','Color',clrBlue,'LineWidth',lw);
        ap_mean=mean(offset_in,'omitnan')-mean(exponent_in,'omitnan').*log_f;
        hSlope=plot(ax3,freqs_in,ap_mean,'--','Color',clrBlue,'LineWidth',lw-0.4);
        if ~isempty(psd_corrected_in)
            log_pc=log10(psd_corrected_in); mu_c=mean(log_pc,1,'omitnan'); sd_c=std(log_pc,0,1,'omitnan');
            fill(ax3,[freqs_in fliplr(freqs_in)],[mu_c+sd_c fliplr(mu_c-sd_c)],...
                 clrOrang,'FaceAlpha',fAlpha,'EdgeColor','none');
            hCor  =plot(ax3,freqs_in,mu_c,'-','Color',clrOrang,'LineWidth',lw);
            hSlopeC=plot(ax3,freqs_in,zeros(size(freqs_in)),'--','Color',clrOrang,'LineWidth',lw-0.4);
            legend(ax3,[hRaw hSlope hCor hSlopeC],...
                   {'Raw PSD','1/f slope (raw)','Corrected PSD','1/f slope (corrected)'},...
                   'Location','northeast','Box','off','FontSize',fsAx-1);
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


%% =========================================================
%% ALL OTHER MEASURES
%% =========================================================

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
    %% -------------------------------------------------------
    %% Multiscale + multichannel
    %% Layout: 3 rows x 4 cols
    %%   Cols 1-2 (all rows): heatmap
    %%   Col 3-4, row k:      topo | curve for cluster k
    %%   Up to 3 clusters shown; unused rows are blanked.
    %% -------------------------------------------------------

    finiteMask=isfinite(entropyData);
    if ~any(finiteMask(:))
        warning('ascent_plot: all values are NaN; nothing to plot.'); return
    end

    % Cluster detection (up to 3)
    [clusterMasks,nClusters,centroids,usedCluster] = ...
        find_entropy_clusters(entropyData, clusterThresh, 3);
    nShow = min(nClusters,3);

    % Contour colours: white / yellow / cyan (for heatmap overlay on dark bg)
    clrBoundary = {[1 1 1], [1 0.90 0.10], [0 0.90 0.90]};
    % Line colours for curve plots on white background (must stay visible)
    clrLine     = {[0.18 0.45 0.87], [0.80 0.45 0], [0 0.55 0.55]};

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

    % Overlay cluster bounding boxes (colour-coded).
    % A bounding box (min-to-max channel row, min-to-max scale) is used
    % instead of contour() to avoid multiple separate loops appearing when
    % channel membership inside a scale-range segment is non-contiguous.
    if usedCluster
        hold(ax_heat,'on');
        for k = 1:nShow
            rows_in = find(any(clusterMasks(:,:,k), 2));   % channels present
            cols_in = find(any(clusterMasks(:,:,k), 1));   % scales present
            if isempty(rows_in) || isempty(cols_in), continue; end
            rectangle(ax_heat, ...
                'Position', [cols_in(1)-0.5,  rows_in(1)-0.5, ...
                             cols_in(end)-cols_in(1)+1, rows_in(end)-rows_in(1)+1], ...
                'EdgeColor', clrBoundary{k}, 'LineWidth', 2, 'FaceColor', 'none', ...
                'Curvature', [0.12 0.12]);
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

    for k=1:3
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
        if iscell(scales)
            scLabel_k = scales{ps_k};
        else
            scLabel_k = num2str(scales(ps_k));
        end

        if usedCluster
            centStr = ifelse(isICA, sprintf('IC%d sc%s',pc_k,scLabel_k), ...
                                    sprintf('%s sc%s',chanlocs(pc_k).labels,scLabel_k));
            rowLabel = sprintf('Cluster %d  (%d %s x %d sc, %s)', k, nCh_k, ...
                               ifelse(isICA,'ICs','ch'), nSc_k, centStr);
        else
            rowLabel = ifelse(isICA, ...
                sprintf('Peak  IC%d, sc %s',pc_k,scLabel_k), ...
                sprintf('Peak  %s, sc %s',chanlocs(pc_k).labels,scLabel_k));
        end

        % Console summary (once per cluster)
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

        %% -- Topo / bar --
        finV = topo_vals_k(isfinite(topo_vals_k));
        if isICA
            axes(ax_topo); %#ok<LAXES>
            bar(ax_topo, topo_vals_k,'FaceColor',[0.18 0.45 0.87],'EdgeColor','none');
            hold(ax_topo,'on');
            if usedCluster
                ic_idx=find(chanMask_k);
                bar(ax_topo,ic_idx,topo_vals_k(ic_idx),...
                    'FaceColor',[0.05 0.25 0.65],'EdgeColor','none');
            end
            bar(ax_topo,pc_k,topo_vals_k(pc_k),'FaceColor','r','EdgeColor','none');
            xlabel(ax_topo,'IC'); ylabel(ax_topo,'Entropy');
            title(ax_topo,rowLabel,'Interpreter','none','FontSize',9);
            box(ax_topo,'on');
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

        %% -- Curve --
        axes(ax_curve); %#ok<LAXES>
        hold(ax_curve,'on'); box(ax_curve,'on');
        if all(~isfinite(curve_k)), curve_k=nan(1,nScales); end
        if isnumeric(scales)
            xvals=scales;
        else
            xvals=1:nScales;
        end

        % Shade cluster scale extent only when clusters differ in scale coverage
        if doShading && isnumeric(scales) && any(scaleMask_k)
            x_lo = xvals(find(scaleMask_k,1,'first'));
            x_hi = xvals(find(scaleMask_k,1,'last'));
            finC = curve_k(isfinite(curve_k));
            if ~isempty(finC)
                yl_pre=[min(finC)*0.95 max(finC)*1.05];
                if yl_pre(1)==yl_pre(2), yl_pre=yl_pre+[-0.1 0.1]; end
                patch(ax_curve,[x_lo x_hi x_hi x_lo],...
                      [yl_pre(1) yl_pre(1) yl_pre(2) yl_pre(2)],...
                      clrLine{k},'FaceAlpha',0.15,'EdgeColor','none');
            end
        end

        plot(ax_curve,xvals,curve_k,'LineWidth',2,'Color',clrLine{k});
        xlim(ax_curve,[xvals(1) xvals(end)]);
        xlabel(ax_curve,'Scale'); ylabel(ax_curve,'Entropy');
        curveTitle = ifelse(usedCluster,...
            sprintf('C%d mean  (%d %s)',k,nCh_k,ifelse(isICA,'ICs','ch')),...
            'Peak channel curve');
        title(ax_curve,curveTitle,'Interpreter','none','FontSize',9);
    end

    colormap(hFig2,parula);
    set(findall(hFig2,'type','axes'),'FontSize',9,'FontWeight','bold');

elseif ~multiscale && multiChan
    %% --- Uniscale topo ---
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
    %% --- Single channel/IC curve ---
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
    %% --- Scalar ---
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


%% =========================================================
%% CLUSTER DETECTION  (scale-first, channel-order agnostic)
%% =========================================================
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
    [rr, cc_] = find(mk);
    w  = data(mk); wS = sum(w, 'omitnan');
    if wS > 0 && isfinite(wS)
        pc = round(sum(rr  .* w, 'omitnan') / wS);
        ps = round(sum(cc_ .* w, 'omitnan') / wS);
    else
        [~, mi] = max(data(mk)); idxs = find(mk);
        [pc, ps] = ind2sub([nChan nScale], idxs(mi));
    end
    centroids(k,:) = [max(1,min(nChan,pc)), max(1,min(nScale,ps))];
end
nClusters = nFound;
end


%% =========================================================
%% APERIODIC TIME-COURSE SUBFUNCTION
%% =========================================================
function plot_aperiodic_timecourse(exponent_t,offset_t,freqs,psd_t,psd_corr_t,times)
bands={'Delta','Theta','Alpha','Beta','Gamma'};
bandEdges=[1 4;4 8;8 13;13 30;30 40];
fMin=freqs(1); fMax=freqs(end);
keep=bandEdges(:,2)>fMin & bandEdges(:,1)<fMax;
bands=bands(keep); bandEdges=bandEdges(keep,:); nBands=numel(bands);
if nBands==0
    warning('ascent_plot: no standard bands within freqs range [%.1f %.1f Hz].',fMin,fMax); return
end
[nChan,~,nTimes]=size(psd_t);
hasCorrected=~isempty(psd_corr_t)&&size(psd_corr_t,3)==nTimes;
if isempty(times)||numel(times)~=nTimes
    times=1:nTimes; xLabel='Window';
else
    times=times(:)'; xLabel='Time (s)';
end
exp_mean=mean(exponent_t,1,'omitnan'); off_mean=mean(offset_t,1,'omitnan');
clrRaw=[0.18 0.45 0.87]; clrCorr=[0.90 0.38 0.15];
clrExp=[0.15 0.15 0.15]; clrOff=[0.50 0.50 0.50]; fAlpha=0.15;
figH=max(420,80+nBands*130);
hFig=figure('Color','w','Position',[80 50 920 figH],...
    'Name','Aperiodic - Time course','NumberTitle','off','Toolbar','none','Menu','none');
try icadefs; set(hFig,'color',BACKCOLOR); catch; end
for iBand=1:nBands
    ax=subplot(nBands,1,iBand);
    fMask=freqs>=bandEdges(iBand,1)&freqs<=bandEdges(iBand,2);
    if ~any(fMask), title(ax,sprintf('%s (no data)',bands{iBand})); axis(ax,'off'); continue; end
    psd_band=log10(psd_t(:,fMask,:)+eps);
    bp_raw=squeeze(mean(psd_band,2,'omitnan')); if nChan==1, bp_raw=bp_raw(:)'; end
    mu_raw=mean(bp_raw,1,'omitnan'); sd_raw=std(bp_raw,0,1,'omitnan');
    yyaxis(ax,'left'); hold(ax,'on');
    fill(ax,[times fliplr(times)],[mu_raw+sd_raw fliplr(mu_raw-sd_raw)],...
         clrRaw,'FaceAlpha',fAlpha,'EdgeColor','none','HandleVisibility','off');
    hRaw=plot(ax,times,mu_raw,'-','Color',clrRaw,'LineWidth',1.8);
    ax.YColor=clrRaw; ylabel(ax,'log_{10} \muV^2/Hz','FontSize',8.5);
    if hasCorrected
        psd_cband=log10(psd_corr_t(:,fMask,:)+eps);
        bp_corr=squeeze(mean(psd_cband,2,'omitnan')); if nChan==1, bp_corr=bp_corr(:)'; end
        mu_corr=mean(bp_corr,1,'omitnan'); sd_corr=std(bp_corr,0,1,'omitnan');
        fill(ax,[times fliplr(times)],[mu_corr+sd_corr fliplr(mu_corr-sd_corr)],...
             clrCorr,'FaceAlpha',fAlpha,'EdgeColor','none','HandleVisibility','off');
        hCorr=plot(ax,times,mu_corr,'-','Color',clrCorr,'LineWidth',1.8);
    end
    yyaxis(ax,'right'); hold(ax,'on');
    hExp=plot(ax,times,exp_mean,'--','Color',clrExp,'LineWidth',0.8);
    hOff=plot(ax,times,off_mean,':','Color',clrOff,'LineWidth',0.8);
    ax.YColor=[0.3 0.3 0.3]; ylabel(ax,'Exp. / Offset','FontSize',8.5);
    xlim(ax,[times(1) times(end)]); box(ax,'on'); ax.TickDir='out';
    title(ax,sprintf('%s  (%g-%g Hz)',bands{iBand},bandEdges(iBand,1),bandEdges(iBand,2)),...
          'FontSize',10,'FontWeight','bold');
    if iBand==1
        if hasCorrected
            legend(ax,[hRaw hCorr hExp hOff],{'Raw PSD','Corrected PSD','Exponent','Offset'},...
                'Location','northeast','Box','off','FontSize',8);
        else
            legend(ax,[hRaw hExp hOff],{'Raw PSD','Exponent','Offset'},...
                'Location','northeast','Box','off','FontSize',8);
        end
    end
    if iBand==nBands, xlabel(ax,xLabel,'FontSize',10); end
end
set(findall(hFig,'type','axes'),'FontSize',9,'FontWeight','bold');
try
    sgtitle(hFig,sprintf('Aperiodic time course  (n = %d channel(s),  mean +/- SD)',nChan),...
        'FontSize',11,'FontWeight','bold');
catch; end
end


%% --- Local helper ---
function out = ifelse(cond,a,b)
if cond, out=a; else, out=b; end
end