function ascent_plot(entropyData, chanlocs, entropyType, scales, varargin)
% ascent_plot  Visualize uni/multi/time-resolved entropy with robust NaN handling.
%
% entropyData : [chan x scale] or [chan x scale x time]
% chanlocs    : EEGLAB channel locations
% entropyType : label, e.g., 'RCMFEσ' or 'MSE (std)'
% scales      : cellstr of scale labels or numeric 1:S
% varargin    : context-dependent optional inputs:
%   ICA name-value pairs (can appear anywhere in varargin):
%       'ICA',     true/false        — enable ICA back-projection mode
%       'icawinv', EEG.icawinv       — [nChans x nICs] inverse weight matrix
%   Aperiodic mode (entropyType contains 'aperiodic'):
%       varargin{1} = offset        [chan x 1]
%       varargin{2} = freqs         [1 x freq]
%       varargin{3} = psd           [chan x freq]  linear
%       varargin{4} = psd_corrected [chan x freq]  or [] if not computed
%   Time-resolved mode:
%       varargin{1} = time_sec      [1 x time]

%% --- Strip ICA name-value pairs before positional varargin parsing ---
isICA   = false;
icawinv = [];
rmIdx   = [];
for vi = 1:2:(numel(varargin)-1)
    if ischar(varargin{vi}) || isstring(varargin{vi})
        switch lower(char(varargin{vi}))
            case 'ica'
                isICA = logical(varargin{vi+1});
                rmIdx = [rmIdx vi vi+1];
            case 'icawinv'
                icawinv = varargin{vi+1};
                rmIdx   = [rmIdx vi vi+1];
        end
    end
end
varargin(rmIdx) = [];

% Validate ICA inputs
if isICA
    if isempty(icawinv)
        error('ascent_plot: ICA=true requires ''icawinv'' (EEG.icawinv, [nChans x nICs]).');
    end
    nICs_winv = size(icawinv, 2);
    nICs_data = size(entropyData, 1);
    if nICs_winv ~= nICs_data
        error('ascent_plot: icawinv has %d ICs but entropyData has %d rows. They must match.', ...
              nICs_winv, nICs_data);
    end
end

%% --- Parse positional varargin ---
time_sec         = [];
offset_in        = [];
freqs_in         = [];
psd_in           = [];
psd_corrected_in = [];

isAperiodic = contains(lower(entropyType), 'aperiodic');

if isAperiodic
    if numel(varargin) >= 1, offset_in        = varargin{1}; end
    if numel(varargin) >= 2, freqs_in         = varargin{2}; end
    if numel(varargin) >= 3, psd_in           = varargin{3}; end
    if numel(varargin) >= 4, psd_corrected_in = varargin{4}; end
else
    if ~isempty(varargin), time_sec = varargin{1}; end
end

%% =========================================================
%% APERIODIC BRANCH
%% =========================================================
if isAperiodic

    exponent_in = entropyData(:);   % [nIC/nChan x 1]

    % Style constants
    clrBlue  = [0.18 0.45 0.87];
    clrOrang = [0.90 0.38 0.15];
    clrReg   = [0.65 0.08 0.08];
    fAlpha   = 0.18;
    lw       = 2;
    fsAx     = 11;
    fsTtl    = 12;
    bgCol    = [0.97 0.97 0.97];

    hFig = figure('Color','w','Position',[100 50 820 980], ...
                  'Name','Aperiodic visualization','NumberTitle','off', ...
                  'Toolbar','none','Menu','none');
    try icadefs; set(hFig,'color',BACKCOLOR); catch; end

    % --- Row 1: exponent topo (col 1) | offset topo (col 2) ---
    ax1 = subplot(4,2,1);
    finE = exponent_in(isfinite(exponent_in));
    if ~isempty(finE)
        if isICA
            bp_exp = icawinv * exponent_in(:);
            axes(ax1); topoplot(bp_exp, chanlocs, 'emarker',{'.','k',7,1},'electrodes','on');
            finBP = bp_exp(isfinite(bp_exp));
            if min(finBP) < max(finBP), clim(ax1,[min(finBP) max(finBP)]); end
        else
            axes(ax1); topoplot(exponent_in, chanlocs, 'emarker',{'.','k',7,1},'electrodes','on');
            if min(finE) < max(finE), clim(ax1,[min(finE) max(finE)]); end
        end
        c1 = colorbar(ax1,'Location','eastoutside');
        c1.Label.String = 'Exponent'; c1.Label.FontSize = fsAx-1; c1.TickDirection = 'out';
    end
    title(ax1, ifelse(isICA,'Exponent (back-projected)','Aperiodic Exponent'), ...
          'FontSize',fsTtl,'FontWeight','bold');

    ax2 = subplot(4,2,2);
    if ~isempty(offset_in)
        finO = offset_in(isfinite(offset_in));
        if ~isempty(finO)
            if isICA
                bp_off = icawinv * offset_in(:);
                axes(ax2); topoplot(bp_off, chanlocs, 'emarker',{'.','k',7,1},'electrodes','on');
                finBO = bp_off(isfinite(bp_off));
                if min(finBO) < max(finBO), clim(ax2,[min(finBO) max(finBO)]); end
            else
                axes(ax2); topoplot(offset_in, chanlocs, 'emarker',{'.','k',7,1},'electrodes','on');
                if min(finO) < max(finO), clim(ax2,[min(finO) max(finO)]); end
            end
            c2 = colorbar(ax2,'Location','eastoutside');
            c2.Label.String = 'Offset'; c2.Label.FontSize = fsAx-1; c2.TickDirection = 'out';
        end
    end
    title(ax2, ifelse(isICA,'Offset (back-projected)','Aperiodic Offset'), ...
          'FontSize',fsTtl,'FontWeight','bold');

    % --- Row 2: scatter exponent vs offset ---
    ax4 = subplot(4,2,[3 4]); hold(ax4,'on');
    if ~isempty(offset_in)
        mask = isfinite(exponent_in) & isfinite(offset_in);
        scatter(ax4, exponent_in(mask), offset_in(mask), 50, clrBlue, 'filled', ...
                'MarkerFaceAlpha',0.65,'MarkerEdgeColor','none');
        if isICA
            [~, ic_max] = max(exponent_in);
            scatter(ax4, exponent_in(ic_max), offset_in(ic_max), 100, 'r', ...
                    'filled','MarkerEdgeColor','k','LineWidth',1.2);
            text(ax4, exponent_in(ic_max), offset_in(ic_max), ...
                 sprintf('  IC%d',ic_max),'FontSize',fsAx-1,'Color','r');
        end
        if sum(mask) > 2
            cc = corrcoef(exponent_in(mask), offset_in(mask));
            r  = cc(1,2);
            lm = polyfit(exponent_in(mask), offset_in(mask), 1);
            xr = linspace(min(exponent_in(mask)), max(exponent_in(mask)), 60);
            plot(ax4, xr, polyval(lm,xr), '-','Color',clrReg,'LineWidth',1.8);
            title(ax4, sprintf('Exponent vs Offset  (r = %.2f)',r),'FontSize',fsTtl,'FontWeight','bold');
        else
            title(ax4,'Exponent vs Offset','FontSize',fsTtl,'FontWeight','bold');
        end
        xlabel(ax4,'Exponent','FontSize',fsAx);
        ylabel(ax4,'Offset','FontSize',fsAx);
        box(ax4,'on'); ax4.TickDir = 'out';
        set(ax4,'Color',bgCol,'LineWidth',0.8);
    end

    % --- Row 3: PSD line plot ---
    if ~isempty(psd_in) && ~isempty(freqs_in)
        log_psd = log10(psd_in);
        mu_raw  = mean(log_psd, 1, 'omitnan');
        sd_raw  = std(log_psd,  0, 1, 'omitnan'); %#ok<NASGU>
        log_f   = log10(freqs_in);

        ax3 = subplot(4,2,[5 6]); hold(ax3,'on');
        fill(ax3, [freqs_in fliplr(freqs_in)], ...
             [mu_raw+std(log_psd,0,1,'omitnan'), fliplr(mu_raw-std(log_psd,0,1,'omitnan'))], ...
             clrBlue,'FaceAlpha',fAlpha,'EdgeColor','none');
        hRaw   = plot(ax3, freqs_in, mu_raw, '-', 'Color',clrBlue,'LineWidth',lw);
        ap_mean = mean(offset_in,'omitnan') - mean(exponent_in,'omitnan') .* log_f;
        hSlope  = plot(ax3, freqs_in, ap_mean, '--','Color',clrBlue,'LineWidth',lw-0.4);

        if ~isempty(psd_corrected_in)
            log_pc = log10(psd_corrected_in);
            mu_c   = mean(log_pc, 1, 'omitnan');
            sd_c   = std(log_pc,  0, 1, 'omitnan');
            fill(ax3, [freqs_in fliplr(freqs_in)], ...
                 [mu_c+sd_c, fliplr(mu_c-sd_c)], ...
                 clrOrang,'FaceAlpha',fAlpha,'EdgeColor','none');
            hCor    = plot(ax3, freqs_in, mu_c,  '-', 'Color',clrOrang,'LineWidth',lw);
            hSlopeC = plot(ax3, freqs_in, zeros(size(freqs_in)), '--','Color',clrOrang,'LineWidth',lw-0.4);
            legend(ax3, [hRaw hSlope hCor hSlopeC], ...
                   {'Raw PSD','1/f slope (raw)','Corrected PSD','1/f slope (corrected)'}, ...
                   'Location','northeast','Box','off','FontSize',fsAx-1);
            title(ax3,'PSD: raw vs corrected  (mean \pm SD)','FontSize',fsTtl,'FontWeight','bold');
        else
            legend(ax3, [hRaw hSlope], {'Raw PSD','1/f slope'}, ...
                   'Location','northeast','Box','off','FontSize',fsAx-1);
            title(ax3,'Raw PSD + aperiodic fit  (mean \pm SD)','FontSize',fsTtl,'FontWeight','bold');
        end
        xlabel(ax3,'Frequency (Hz)','FontSize',fsAx);
        ylabel(ax3,'Power (log_{10} \muV^2/Hz)','FontSize',fsAx);
        xlim(ax3,[freqs_in(1) freqs_in(end)]);
        box(ax3,'on'); ax3.TickDir = 'out';
        set(ax3,'Color',bgCol,'LineWidth',0.8);
    end

    % --- Row 4: alpha power topos ---
    ax_pr = subplot(4,2,7);
    ax_pc = subplot(4,2,8);
    if ~isempty(psd_in) && ~isempty(freqs_in)
        aMask = freqs_in >= 8 & freqs_in <= 13;
        if any(aMask)
            alpha_raw = mean(log10(psd_in(:,aMask)), 2, 'omitnan');
            if isICA
                bp_alpha = icawinv * alpha_raw(:);
                finP = bp_alpha(isfinite(bp_alpha));
                if ~isempty(finP)
                    axes(ax_pr); topoplot(bp_alpha, chanlocs, 'emarker',{'.','k',7,1},'electrodes','on');
                    if min(finP) < max(finP), clim(ax_pr,[min(finP) max(finP)]); end
                    c3 = colorbar(ax_pr,'Location','eastoutside');
                    c3.Label.String = 'log_{10} \muV^2/Hz'; c3.Label.FontSize = fsAx-1; c3.TickDirection = 'out';
                end
                title(ax_pr,'Alpha power (back-projected)','FontSize',fsTtl,'FontWeight','bold');

                [~, ic_max] = max(exponent_in);
                axes(ax_pc); topoplot(icawinv(:,ic_max), chanlocs, 'emarker',{'.','k',7,1},'electrodes','on');
                finIC = icawinv(:,ic_max); finIC = finIC(isfinite(finIC));
                if min(finIC) < max(finIC), clim(ax_pc,[min(finIC) max(finIC)]); end
                c4 = colorbar(ax_pc,'Location','eastoutside');
                c4.Label.String = 'Weight'; c4.Label.FontSize = fsAx-1; c4.TickDirection = 'out';
                title(ax_pc, sprintf('IC%d scalp map (max exponent)',ic_max), ...
                      'FontSize',fsTtl,'FontWeight','bold');
            else
                finP = alpha_raw(isfinite(alpha_raw));
                if ~isempty(finP)
                    axes(ax_pr); topoplot(alpha_raw, chanlocs, 'emarker',{'.','k',7,1},'electrodes','on');
                    if min(finP) < max(finP), clim(ax_pr,[min(finP) max(finP)]); end
                    c3 = colorbar(ax_pr,'Location','eastoutside');
                    c3.Label.String = 'log_{10} \muV^2/Hz'; c3.Label.FontSize = fsAx-1; c3.TickDirection = 'out';
                end
                title(ax_pr,'Alpha power (8–13 Hz) — raw','FontSize',fsTtl,'FontWeight','bold');

                if ~isempty(psd_corrected_in)
                    alpha_cor = mean(log10(psd_corrected_in(:,aMask)), 2, 'omitnan');
                    finPC = alpha_cor(isfinite(alpha_cor));
                    if ~isempty(finPC)
                        axes(ax_pc); topoplot(alpha_cor, chanlocs, 'emarker',{'.','k',7,1},'electrodes','on');
                        if min(finPC) < max(finPC), clim(ax_pc,[min(finPC) max(finPC)]); end
                        c4 = colorbar(ax_pc,'Location','eastoutside');
                        c4.Label.String = 'log_{10} \muV^2/Hz'; c4.Label.FontSize = fsAx-1; c4.TickDirection = 'out';
                    end
                    title(ax_pc,'Alpha power (8–13 Hz) — corrected','FontSize',fsTtl,'FontWeight','bold');
                else
                    axes(ax_pc); axis off;
                    text(0.5,0.5,'Correction not applied', ...
                         'HorizontalAlignment','center','VerticalAlignment','middle', ...
                         'FontSize',fsAx,'Color',[0.55 0.55 0.55]);
                end
            end
        else
            warning('ascent_plot: alpha band (8-13 Hz) not covered by freqs_in — skipping alpha topos.');
            axes(ax_pr); axis off;
            axes(ax_pc); axis off;
        end
    else
        axes(ax_pr); axis off;
        axes(ax_pc); axis off;
    end

    % Set colormaps after all topoplot calls to avoid overwrite
    % colormap(ax1,  parula);
    % colormap(ax2,  parula);
    % colormap(ax_pr,parula);
    % colormap(ax_pc,parula);
    colormap(hFig, parula);

    % Space out subplots to avoid title/xlabel overlap
    allAx = findall(hFig, 'type', 'axes');
    for k = 1:numel(allAx)
        pos = allAx(k).Position;
        % shrink height by 4%, shift up by 2%
        allAx(k).Position = [pos(1), pos(2)+0.02, pos(3), pos(4)-0.04];
    end

    colormap(hFig, parula);
    set(findall(hFig,'type','axes'),'FontSize',fsAx,'FontWeight','bold');

    return
end


%% =========================================================
%% ALL OTHER MEASURES
%% =========================================================

isTimeResolved = ndims(entropyData) == 3;
if isTimeResolved
    entropyData3D = entropyData;
    entropyData   = mean(entropyData3D, 3, 'omitnan');
end
multiscale = size(entropyData,2) > 1;

if ~multiscale
    entropyData(entropyData==0) = NaN;
end

nChan     = size(entropyData,1);
multiChan = nChan > 1;

if multiscale && multiChan
    %% --- Multiscale heatmap ---
    hFig2 = figure('Color','w','InvertHardCopy','off', ...
                   'Name','Multiscale entropy visualization', ...
                   'Toolbar','none','Menu','none','NumberTitle','Off');
    try icadefs; set(hFig2,'color',BACKCOLOR); catch; end

    nScales = size(entropyData,2);

    ax_heat = subplot(3,3,[1 2 4 5 7 8]);
    imagesc(ax_heat, 1:nScales, 1:nChan, entropyData); axis(ax_heat,'tight');
    set(ax_heat,'TickDir','out'); box(ax_heat,'on');

    if isICA
        Yticks = arrayfun(@(x){sprintf('IC%d',x)}, 1:nChan);
    else
        Yticks = {chanlocs.labels};
    end
    newY = 1:nChan;
    if nChan > 30, newY = round(linspace(1,nChan,20)); end
    set(ax_heat,'YTick',newY,'YTickLabel',Yticks(newY),'FontWeight','normal');

    if iscell(scales)
        Xticks = scales; nX = numel(scales);
    else
        Xticks = arrayfun(@(x){num2str(x)},1:nScales); nX = nScales;
    end
    newX = 1:nX;
    if nX > 30, newX = round(linspace(1,nX,20)); end
    set(ax_heat,'XTick',newX,'XTickLabel',Xticks(newX),'FontWeight','normal');
    ax_heat.TickLabelInterpreter = 'none';
    ax_heat.XTickLabelRotation   = 45;
    ax_heat.PositionConstraint   = 'outerposition';
    drawnow;
    outer     = ax_heat.OuterPosition;
    ti        = ax_heat.TightInset;
    newBottom = max(outer(2), ti(2) + 0.02);
    newHeight = max(outer(4) - (newBottom - outer(2)) - (ti(4) + 0.01), 0.1);
    ax_heat.OuterPosition = [outer(1), newBottom, outer(3), newHeight];

    colormap(ax_heat,'parula');
    c = colorbar(ax_heat); ylabel(c,'Entropy','FontWeight','bold','FontSize',9);
    xlabel(ax_heat,'Scales');
    ylabel(ax_heat, ifelse(isICA,'ICs','EEG channels'));
    title(ax_heat, entropyType,'Interpreter','none');

    % Find peak
    finiteMask = isfinite(entropyData);
    if ~any(finiteMask(:))
        warning('ascent_plot: all values are NaN; nothing to plot.');
        return
    end
    tmp = entropyData; tmp(~finiteMask) = -Inf;
    [peak_value, linear_idx] = max(tmp(:));
    [peak_channel, peak_scale] = ind2sub(size(tmp), linear_idx);
    if ~any(isfinite(entropyData(:,peak_scale)))
        [~, peak_scale] = max(sum(isfinite(entropyData),1));
    end
    if iscell(scales), sclabel = scales{peak_scale}; else, sclabel = num2str(peak_scale); end

    if isICA
        fprintf('Peak entropy: %.3f at Scale %s (index %d), IC%d.\n', ...
            peak_value, sclabel, peak_scale, peak_channel);
    else
        fprintf('Peak entropy: %.3f at Scale %s (index %d), Channel %s.\n', ...
            peak_value, sclabel, peak_scale, chanlocs(peak_channel).labels);
    end

    % Subplot 6: per-channel/IC entropy curve at peak scale
    ax6 = subplot(3,3,6); hold(ax6,'on'); box(ax6,'on');
    row = entropyData(peak_channel,:);
    if all(~isfinite(row)), row = nan(1,nScales); end
    plot(ax6, 1:nScales, row, 'LineWidth',2);
    xlim(ax6,[1 nScales]);
    xlabel(ax6,'Scale'); ylabel(ax6,'Entropy');
    title(ax6, ifelse(isICA, sprintf('IC%d',peak_channel), ...
          sprintf('Channel %s',chanlocs(peak_channel).labels)), 'Interpreter','none');

    % Subplot 3: bar chart (ICA) | entropy topo at peak scale (chan)
    vals = entropyData(:, peak_scale);
    ax3  = subplot(3,3,3);
    if isICA
        bar(ax3, vals, 'FaceColor',[0.18 0.45 0.87],'EdgeColor','none');
        hold(ax3,'on');
        bar(ax3, peak_channel, vals(peak_channel), 'FaceColor','r','EdgeColor','none');
        xlabel(ax3,'IC'); ylabel(ax3,'Entropy');
        title(ax3, sprintf('All ICs @ scale %s',sclabel),'Interpreter','none');
        box(ax3,'on');
    else
        finiteVals = vals(isfinite(vals));
        if isempty(finiteVals)
            axis(ax3,'off');
            text(0.5,0.5,sprintf('No finite values @ scale %s',sclabel), ...
                'Parent',ax3,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontWeight','bold');
        else
            try
                axes(ax3); topoplot(vals, chanlocs, 'emarker',{'.','k',8,1},'electrodes','on');
                lo = min(finiteVals); hi = max(finiteVals);
                if isfinite(lo) && isfinite(hi) && lo < hi, clim(ax3,[lo hi]); end
                colormap(ax3,'parula');
                title(ax3, sprintf('%s @ scale %s',entropyType,sclabel),'Interpreter','none');
            catch
                bar(ax3, vals); xlim(ax3,[0 numel(vals)+1]); box(ax3,'on');
                title(ax3, sprintf('Entropy @ scale %s',sclabel),'Interpreter','none');
                ylabel(ax3,'Entropy');
            end
        end
    end

    % Subplot 9: time course | back-projected topo (ICA) | off (chan)
    ax9 = subplot(3,3,9);
    if isTimeResolved
        hold(ax9,'on'); box(ax9,'on');
        tr = squeeze(entropyData3D(peak_channel, peak_scale, :));
        if isempty(tr) || all(~isfinite(tr))
            axis(ax9,'off');
            text(0.5,0.5,'No finite time-resolved values', ...
                'Parent',ax9,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontWeight','bold');
        else
            if isempty(time_sec) || numel(time_sec) ~= numel(tr)
                t = 1:numel(tr); xlabel(ax9,'Window #');
            else
                t = time_sec(:); xlabel(ax9,'Time (s)');
            end
            plot(ax9, t, tr, 'LineWidth',1.5);
            ylabel(ax9,'Entropy');
            title(ax9, ifelse(isICA, ...
                  sprintf('IC%d, scale %s — time course',peak_channel,sclabel), ...
                  sprintf('Time course @ %s, scale %s',chanlocs(peak_channel).labels,sclabel)), ...
                  'Interpreter','none');
        end
    elseif isICA
        bp_vals = icawinv * vals(:);
        finBP   = bp_vals(isfinite(bp_vals));
        if ~isempty(finBP)
            axes(ax9); topoplot(bp_vals, chanlocs, 'emarker',{'.','k',8,1},'electrodes','on');
            lo = min(finBP); hi = max(finBP);
            if isfinite(lo) && isfinite(hi) && lo < hi, clim(ax9,[lo hi]); end
            colormap(ax9,'parula');
            title(ax9, sprintf('Back-projected entropy @ scale %s',sclabel),'Interpreter','none');
        else
            axis(ax9,'off');
        end
    else
        axis(ax9,'off');
    end

    colormap(hFig2, parula);
    set(findall(hFig2,'type','axes'),'FontSize',10,'FontWeight','bold');

elseif ~multiscale && multiChan
    %% --- Uniscale topo ---
    hFig3 = figure('Color','w','InvertHardCopy','off', ...
                   'Name','Uniscale entropy visualization', ...
                   'Toolbar','none','Menu','none','NumberTitle','Off');
    try icadefs; set(hFig3,'color',BACKCOLOR); catch; end

    vals       = entropyData(:);
    finiteVals = vals(isfinite(vals));
    ax_uni     = axes(hFig3);
    if isempty(finiteVals)
        axis(ax_uni,'off');
        text(0.5,0.5,'No finite values to plot','Parent',ax_uni, ...
             'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold');
    else
        if isICA
            bp_vals = icawinv * vals(:);
            finBP   = bp_vals(isfinite(bp_vals));
            axes(ax_uni); topoplot(bp_vals, chanlocs, 'emarker',{'.','k',15,1},'electrodes','labels');
            if min(finBP) < max(finBP), clim(ax_uni,[min(finBP)*0.95 max(finBP)*1.05]); end
            title(ax_uni,[entropyType ' (back-projected)'],'Interpreter','none');
        else
            axes(ax_uni); topoplot(vals, chanlocs, 'emarker',{'.','k',15,1},'electrodes','labels');
            if min(finiteVals) < max(finiteVals)
                clim(ax_uni,[min(finiteVals)*0.95 max(finiteVals)*1.05]);
            end
            title(ax_uni, entropyType,'Interpreter','none');
        end
        % colormap(ax_uni,'parula');
        colormap(hFig3, parula);
        c = colorbar(ax_uni);
        c.Label.String = 'Entropy'; c.Label.FontSize = 11; c.Label.FontWeight = 'bold';
    end
    set(findall(hFig3,'type','axes'),'FontSize',10,'FontWeight','bold');

elseif multiscale && ~multiChan
    %% --- Single channel/IC curve ---
    hFig4 = figure('Color','w','InvertHardCopy','off');
    try icadefs; set(hFig4,'color',BACKCOLOR); catch; end
    ax_s = axes(hFig4);
    plot(ax_s, scales, entropyData, 'LineWidth',2);
    ylabel(ax_s,'Entropy'); xlabel(ax_s,'Scales');
    title(ax_s, entropyType,'Interpreter','none');
    axis(ax_s,'tight');
    set(findall(hFig4,'type','axes'),'FontSize',10,'FontWeight','bold');

else
    error('ascent_plot: data format not recognized.')
end
end


%% --- Local helper ---
function out = ifelse(cond, a, b)
if cond, out = a; else, out = b; end
end