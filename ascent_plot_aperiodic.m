function ascent_plot_aperiodic(exponent_mat, chanlocs, opt)
% ascent_plot_aperiodic  Aperiodic (1/f) figures for Ascent.
%
% Called by ascent_plot; not intended for direct use. See ascent_plot for the
% public interface and for what opt carries.
%
% Static layout (4x2):
%   Row 1 : exponent topo | offset topo
%   Row 2 : exponent vs offset scatter
%   Row 3 : PSD raw vs corrected (mean +/- SD)
%   Row 4 : alpha power raw | alpha power corrected
%
% Time-varying layout: see plot_aperiodic_timecourse below.
%
% This function performs no analysis. Alpha power comes from opt.alphaPower
% (computed and saved by ascent_compute); when a direct caller omits it, it is
% derived with compute_AperiodicBandPower -- the same function ascent_compute
% uses -- so the numbers are identical either way.
%
% Copyright (C) - ASCENT EEGLAB PLUGIN - Cedric Cannard, 2021-2025

isICA            = opt.isICA;
icawinv          = opt.icawinv;
offset_in        = opt.offset;
freqs_in         = opt.freqs;
psd_in           = opt.psd;
psd_corrected_in = opt.psd_corrected;
times_in         = opt.times;

isTimeVarying = size(exponent_mat,2) > 1;
if isTimeVarying
    plot_aperiodic_timecourse(exponent_mat,offset_in,freqs_in,psd_in, ...
                              psd_corrected_in,times_in,isICA);
    return
end
exponent_in = exponent_mat(:);

% Top IC = steepest aperiodic exponent. Resolved once here so the scatter
% marker and the alpha topos below always refer to the same component.
ic_top = [];
if isICA, [~,ic_top] = max(exponent_in); end

clrBlue=[0.18 0.45 0.87]; clrOrang=[0.90 0.38 0.15]; clrReg=[0.65 0.08 0.08];
fAlpha=0.18; lw=2; fsAx=12; fsTtl=12; bgCol=[0.97 0.97 0.97];

hFig=figure('Color','w','Position',[100 50 820 980],...
            'Name','Aperiodic visualization','NumberTitle','off',...
            'Toolbar','none','Menu','none');
try icadefs; set(hFig,'color',BACKCOLOR); catch; end %#ok<NODEF>

ax1=subplot(4,2,1);
finE=exponent_in(isfinite(exponent_in));
if ~isempty(finE)
    if isICA
        bp_exp=icawinv*exponent_in(:);
        axes(ax1); topoplot(bp_exp,chanlocs,'emarker',{'.','k',7,1},'electrodes','on'); %#ok<LAXES>
        finBP=bp_exp(isfinite(bp_exp));
        if min(finBP)<max(finBP), clim(ax1,[min(finBP) max(finBP)]); end
    else
        axes(ax1); topoplot(exponent_in,chanlocs,'emarker',{'.','k',7,1},'electrodes','on'); %#ok<LAXES>
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
            axes(ax2); topoplot(bp_off,chanlocs,'emarker',{'.','k',7,1},'electrodes','on'); %#ok<LAXES>
            finBO=bp_off(isfinite(bp_off));
            if min(finBO)<max(finBO), clim(ax2,[min(finBO) max(finBO)]); end
        else
            axes(ax2); topoplot(offset_in,chanlocs,'emarker',{'.','k',7,1},'electrodes','on'); %#ok<LAXES>
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
        scatter(ax4,exponent_in(ic_top),offset_in(ic_top),100,'r',...
                'filled','MarkerEdgeColor','k','LineWidth',1.2);
        text(ax4,exponent_in(ic_top),offset_in(ic_top),...
             sprintf('  IC%d',ic_top),'FontSize',fsAx-1,'Color','r');
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
        lgd = legend(ax3, [hRaw hSlope hCor hSlopeC], ...
               {'Raw PSD', '1/f fit (raw)', 'Corrected PSD', '1/f fit (corrected)'}, ...
               'Location', 'eastoutside', 'Box', 'off', 'FontSize', fsAx-1);
        title(ax3,'PSD: raw vs corrected  (mean +/- SD)','FontSize',fsTtl,'FontWeight','bold');
    else
        lgd = legend(ax3,[hRaw hSlope],{'Raw PSD','1/f fit'},...
               'Location','northeast','Box','off','FontSize',fsAx-1);
        title(ax3,'Raw PSD + aperiodic fit  (mean +/- SD)','FontSize',fsTtl,'FontWeight','bold');
    end
    xlabel(ax3,'Frequency (Hz)','FontSize',fsAx);
    ylabel(ax3,'Power (log_{10} \muV^2/Hz)','FontSize',fsAx);
    xlim(ax3,[freqs_in(1) freqs_in(end)]); box(ax3,'on'); ax3.TickDir='out';
    set(ax3,'Color',bgCol,'LineWidth',0.8);
end

% ---------------------------------------------------------------------
% Bottom row: alpha power before vs after aperiodic correction.
% Both panels are absolute power in uV^2/Hz on a SHARED colour scale, so the
% intensity drop from left to right is the aperiodic contribution at alpha.
% (psd_corrected is deliberately not used here: it is a ratio to the fit, and
% cannot share a scale with an absolute power.)
% ---------------------------------------------------------------------
ax_pr=subplot(4,2,7); ax_pc=subplot(4,2,8);
ap = resolve_alpha_power(opt, exponent_in, offset_in);

if isempty(ap)
    axes(ax_pr); axis off; axes(ax_pc); axis off; %#ok<LAXES>
else
    bandStr = sprintf('%g-%g Hz', ap.band(1), ap.band(2));
    if isICA
        % A single IC's power reaches the scalp through the SQUARED mixing
        % weights (icawinv(:,ic).^2 * P_ic): power adds, amplitude does not.
        w2       = icawinv(:,ic_top).^2;
        proj_raw = safelog10(w2*ap.raw(ic_top));
        proj_osc = safelog10(w2*ap.osc(ic_top));
        unitLbl  = 'log_{10} \muV^2/Hz';
        ttlRaw   = sprintf('IC%d alpha (%s) - raw',       ic_top, bandStr);
        ttlCor   = sprintf('IC%d alpha (%s) - corrected', ic_top, bandStr);
        if ~isempty(ap.osc) && isfinite(ap.osc(ic_top)) && ap.osc(ic_top) > 0
            fprintf(['Top IC (max exponent): IC%d  alpha %.3g -> %.3g uV^2/Hz ' ...
                     '(%.0f%% of alpha-band power was aperiodic)\n'], ...
                    ic_top, ap.raw(ic_top), ap.osc(ic_top), ...
                    100*(1-ap.osc(ic_top)/ap.raw(ic_top)));
        else
            fprintf('Top IC (max exponent): IC%d  alpha %.3g uV^2/Hz (raw)\n', ...
                    ic_top, ap.raw(ic_top));
        end
    else
        proj_raw = safelog10(ap.raw);
        proj_osc = safelog10(ap.osc);
        unitLbl  = 'log_{10} \muV^2/Hz';
        ttlRaw   = sprintf('Alpha power (%s) - raw',       bandStr);
        ttlCor   = sprintf('Alpha power (%s) - corrected', bandStr);
    end

    % One colour scale across both panels
    shared = [proj_raw; proj_osc]; shared = shared(isfinite(shared));
    cl = [];
    if ~isempty(shared) && min(shared) < max(shared), cl = [min(shared) max(shared)]; end

    if any(isfinite(proj_raw))
        axes(ax_pr); topoplot(proj_raw,chanlocs,'emarker',{'.','k',7,1},'electrodes','on'); %#ok<LAXES>
        if ~isempty(cl), clim(ax_pr,cl); end
        c3=colorbar(ax_pr,'Location','eastoutside');
        c3.Label.String=unitLbl; c3.Label.FontSize=fsAx-1; c3.TickDirection='out';
    else
        axes(ax_pr); axis off; %#ok<LAXES>
    end
    title(ax_pr,ttlRaw,'FontSize',fsTtl,'FontWeight','bold');

    if any(isfinite(proj_osc))
        axes(ax_pc); topoplot(proj_osc,chanlocs,'emarker',{'.','k',7,1},'electrodes','on'); %#ok<LAXES>
        if ~isempty(cl), clim(ax_pc,cl); end
        c4=colorbar(ax_pc,'Location','eastoutside');
        c4.Label.String=unitLbl; c4.Label.FontSize=fsAx-1; c4.TickDirection='out';
        title(ax_pc,ttlCor,'FontSize',fsTtl,'FontWeight','bold');
    else
        axes(ax_pc); axis off; %#ok<LAXES>
        if isempty(ap.osc) || all(isnan(ap.osc))
            msg = 'Correction not applied';
        else
            msg = sprintf('1/f fit exceeds\nalpha-band power');
        end
        text(0.5,0.5,msg,'HorizontalAlignment','center',...
             'VerticalAlignment','middle','FontSize',fsAx,'Color',[0.55 0.55 0.55]);
    end
end

colormap(hFig,parula);
allAx=findall(hFig,'type','axes');
for k=1:numel(allAx)
    pos=allAx(k).Position;
    allAx(k).Position=[pos(1),pos(2)+0.02,pos(3),pos(4)-0.04];
end
colormap(hFig,parula);
set(findall(hFig,'type','axes'),'FontSize',fsAx,'FontWeight','bold');

% Scatter and PSD rows share one footprint, with a right margin sized to the
% legend. 'eastoutside' alone shrinks only the PSD axes (leaving the two rows
% different widths) and still clips the labels at this figure width.
drawnow;
rowAx = gobjects(0);
if exist('ax4','var') && isgraphics(ax4), rowAx(end+1) = ax4; end
if exist('ax3','var') && isgraphics(ax3), rowAx(end+1) = ax3; end
if ~isempty(rowAx)
    rowLeft = max(arrayfun(@(a) a.Position(1), rowAx));
    haveLgd = exist('lgd','var') && isgraphics(lgd);
    if haveLgd
        set(lgd,'FontSize',fsAx-1);   % the axes FontSize sweep above can bump it
        drawnow;
        rowWidth = 1 - rowLeft - lgd.Position(3) - 0.045;
    else
        rowWidth = 1 - rowLeft - 0.045;
    end
    rowWidth = max(rowWidth, 0.30);
    for a = rowAx
        pp = a.Position;
        a.Position = [rowLeft, pp(2), rowWidth, pp(4)];
    end
    if haveLgd
        pp = ax3.Position; lp = lgd.Position;
        lgd.Position = [rowLeft+rowWidth+0.015, pp(2)+pp(4)/2-lp(4)/2, lp(3), lp(4)];
    end
end
end


%% Alpha power: prefer what ascent_compute saved, else derive it the same way
function ap = resolve_alpha_power(opt, exponent_in, offset_in)
ap = [];
if ~isempty(opt.alphaPower) && isstruct(opt.alphaPower)
    ap = opt.alphaPower;
    if ~isfield(ap,'band') || isempty(ap.band), ap.band = [8 13]; end
    ap.raw = ap.raw(:);
    ap.osc = ap.osc(:);
    return
end
% Direct caller did not pass it: derive with the same function ascent_compute
% uses, so the figure cannot disagree with a saved result.
if isempty(opt.psd) || isempty(opt.freqs) || isempty(offset_in), return; end
[raw, osc] = compute_AperiodicBandPower(opt.freqs, opt.psd, exponent_in, offset_in);
ap = struct('raw',raw(:), 'osc',osc(:), 'band',[8 13]);
end


%% log10 that maps non-positive/non-finite power to NaN instead of complex/-Inf
function y = safelog10(x)
y = x(:);
y(~isfinite(y) | y <= 0) = NaN;
y = log10(y);
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
    ylabel(cb2, 'log_{10} (power / aperiodic fit)', 'FontSize', 8);
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
    legend(ax4, allH, 'Location', 'eastoutside', 'Box', 'off', 'FontSize', 10);
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
    pp = ax.Position;
    ax.Position = [refLeft, pp(2), refWidth, pp(4)];
end
cbW   = 0.025;
cbGap = 0.008;
for ii = 1:numel(cbAxes)
    pp = cbAxes(ii).Position;
    cbList(ii).Position = [refLeft+refWidth+cbGap, pp(2)+pp(4)*0.1, cbW, pp(4)*0.8];
end
end


%% Local helpers
function out = ifelse(cond,a,b)
if cond, out=a; else, out=b; end
end
