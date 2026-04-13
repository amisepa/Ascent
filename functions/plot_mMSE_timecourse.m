function plot_mMSE_timecourse(mse_time, time_sec, scales, varargin)
% plot_mMSE_timecourse  Plot mMSE time course: mean ± SE across channels.
%
%   plot_mMSE_timecourse(mse_time, time_sec, scales)
%   plot_mMSE_timecourse(mse_time, time_sec, scales, 'Title', 'my title')
%   plot_mMSE_timecourse(mse_time, time_sec, scales, 'ScaleIdx', [2 4 6])
%   plot_mMSE_timecourse(mse_time, time_sec, scales, 'Alpha', 0.25)
%
% Inputs:
%   mse_time  : [nChan x nScales x nTime]  — from info.mse_time
%   time_sec  : [1 x nTime]                — from info.time_sec
%   scales    : cell array of scale labels — from compute_mMSE output
%
% Optional name-value:
%   'ScaleIdx' : indices of scales to plot (default: all)
%   'Title'    : figure title string (default: 'mMSE time course')
%   'Alpha'    : shaded region transparency 0-1 (default: 0.20)
%   'LineWidth': line width (default: 1.8)
%   'YLabel'   : y-axis label (default: 'SampEn')
%   'XLabel'   : x-axis label (default: 'Time (s)')
%
% Output: figure handle (returned implicitly)

p2 = inputParser;
p2.addParameter('ScaleIdx', [], @(x) isnumeric(x) && isvector(x));
p2.addParameter('Title',    'mMSE time course', @ischar);
p2.addParameter('Alpha',    0.20, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p2.addParameter('LineWidth',1.8,  @(x) isnumeric(x) && isscalar(x) && x>0);
p2.addParameter('YLabel',   'SampEn', @ischar);
p2.addParameter('XLabel',   'Time (s)', @ischar);
p2.parse(varargin{:});
opt = p2.Results;

% mse_time: [nChan x nScales x nTime]
[nChan, nScales, nTime] = size(mse_time);

if isempty(opt.ScaleIdx)
    scIdx = 1:nScales;
else
    scIdx = opt.ScaleIdx(opt.ScaleIdx >= 1 & opt.ScaleIdx <= nScales);
end
nPlot = numel(scIdx);

if nPlot == 0
    warning('plot_mMSE_timecourse: no valid scale indices to plot.'); return
end

% Mean and SE across channels (dim 1)
% Extract subset first to avoid repeated indexing: [nChan x nPlot x nTime]
sub = mse_time(:, scIdx, :);
mu  = squeeze(mean(sub, 1, 'omitnan'));       % [nPlot x nTime]
sd  = squeeze(std( sub, 0, 1, 'omitnan'));    % [nPlot x nTime]
n   = squeeze(sum(isfinite(sub), 1));         % [nPlot x nTime]
se  = sd ./ sqrt(max(1, n));                  % [nPlot x nTime]

% Handle edge cases: single scale or single time point
if nPlot == 1, mu = mu(:)'; se = se(:)'; end
if nTime == 1, mu = mu(:);  se = se(:);  end

% Colormap: one color per scale
cmap = parula(nPlot);

figure('Color','w','Name', opt.Title, 'NumberTitle','off');
ax = axes; hold(ax,'on'); box(ax,'on');

for ii = 1:nPlot
    s_idx = scIdx(ii);
    c = cmap(ii,:);
    t = time_sec(:)';
    m = mu(ii,:);
    e = se(ii,:);

    % Shaded SE region — excluded from legend
    tFill = [t, fliplr(t)];
    yFill = [m + e, fliplr(m - e)];
    if any(isfinite(yFill))
        yFill(~isfinite(yFill)) = 0;
        fill(ax, tFill, yFill, c, 'FaceAlpha', opt.Alpha, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
    end

    % Mean line — appears in legend
    if s_idx <= numel(scales)
        lbl = scales{s_idx};
    else
        lbl = sprintf('Scale %d', s_idx);
    end
    plot(ax, t, m, '-', 'Color', c, 'LineWidth', opt.LineWidth, 'DisplayName', lbl);
end
xlabel(ax, opt.XLabel, 'FontSize', 12);
ylabel(ax, opt.YLabel, 'FontSize', 12);
title(ax, opt.Title, 'Interpreter', 'none', 'FontSize', 13, 'FontWeight', 'bold');
legend(ax, 'Location', 'best', 'Box', 'off', 'FontSize', 9, 'Interpreter', 'none');
set(ax, 'FontSize', 11, 'FontWeight', 'bold', 'TickDir', 'out', 'LineWidth', 1.0);
xlim(ax, [time_sec(1) time_sec(end)]);

end % plot_mMSE_timecourse