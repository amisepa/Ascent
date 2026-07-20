function [HD, HA, HDA, info] = compute_ExSEnt(data, varargin)
% compute_ExSEnt  Extrema-Segmented Entropy (ExSEnt) for multichannel data.
%
%   [HD, HA, HDA, info] = compute_ExSEnt(data, 'm', 2, ...
%                                        'r', 0.15, 'lambda', 0.01, ...
%                                        'Parallel', true, 'Progress', true, ...
%                                        'Plot', false)
%
% Inputs:
%   data        : matrix [n_channels x n_samples]
%   'm'         : embedding dimension (default = 2)
%   'r'         : r scaling (r = r * std), default = 0.15
%   'lambda'    : threshold factor for extrema detection (θ = λ * IQR(diff(x))), default = 0.01
%   'Parallel'  : logical true/false to enable parfor over channels (default = true)
%   'Progress'  : logical true/false to show progress (default = true)
%                 • If Parallel==true  and Progress==true → text only (parfor-safe)
%                 • If Parallel==false and Progress==true → text + waitbar (fallback to text if headless)
%   'Plot'      : plot mode (default = false)
%                 false / 'none'      → no plot
%                 true / 'validation' → signal + D/A sequences for channel 1
%                                       (use during lambda tuning on new data)
%                 'bifurcation'       → bifurcation-style scatter across all
%                                       channels, coloured by H_DA
%
% Outputs:
%   HD, HA, HDA : [n_channels x 1] entropies (duration, amplitude, joint)
%   info        : struct with per-channel diagnostics:
%                 .M, .rD, .rA, .rDA, .rangeD, .rangeA, .segment_idx
%
% Notes:
%   • Durations (D) are in samples; amplitudes (A) are net differences.
%   • Joint SampEn uses interleaved z-scored [D A] (embedding = 2*m).
%   • SampEn estimation is delegated to the shared compute_SampEn engine with
%     'ZScore', false, so rD/rA/r are applied as absolute Chebyshev tolerances
%     on the sequences ExSEnt has already prepared (single source of truth).
%
% Example:
%   [HD,HA,HDA,info] = compute_ExSEnt(EEG.data, 'm', 2, 'r', 0.2, ...
%                                     'lambda', 0.01, 'Parallel', true, ...
%                                     'Progress', true, 'Plot', 'bifurcation', ...
%                                     'ChanLabels', {EEG.chanlocs.labels});
%
% Reference:
%   Kamali, S., Baroni, F., & Varona, P. (2025). ExSEnt: Extrema-Segmented
%   Entropy Analysis of Time Series. arXiv:2509.07751
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Ascent plugin — Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

% ---------------- Parse inputs ----------------
p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ndims(x) == 2);
p.addParameter('m', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('r', 0.15, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('lambda', 0.01, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('Parallel', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Plot', false);                % logical | 'none' | 'validation' | 'bifurcation'
p.addParameter('ChanLabels', {}, @iscell);   % optional channel label cell array for bifurcation plot
p.parse(data, varargin{:});

m            = p.Results.m;
r            = p.Results.r;
lambda       = p.Results.lambda;
parallelMode = p.Results.Parallel;
showProgress = p.Results.Progress;
plotArg      = p.Results.Plot;
chanLabels   = p.Results.ChanLabels;

% Normalise plotArg → plotMode string
if islogical(plotArg) || isnumeric(plotArg)
    if plotArg, plotMode = 'validation'; else, plotMode = 'none'; end
else
    plotMode = lower(char(plotArg));   % 'none' | 'validation' | 'bifurcation'
end

% ---------------- Shape ----------------
if size(data,1) > size(data,2)
    data = data.'; % [n_channels x n_samples]
end
[nchan, ~] = size(data);

HD  = nan(nchan,1);
HA  = nan(nchan,1);
HDA = nan(nchan,1);
info = struct('M',nan(nchan,1),'rD',nan(nchan,1),'rA',nan(nchan,1), ...
              'rDA',nan(nchan,1),'rangeD',nan(nchan,1),'rangeA',nan(nchan,1), ...
              'segment_idx',[]);

% ---------------- Progress header ----------------
if showProgress
    if parallelMode
        fprintf('ExSEnt: %d channel(s) | m=%g, r=%.3f, lambda=%.4f | parallel=on (text only)\n', ...
                nchan, m, r, lambda);
        fprintf('Progress:\n');
    else
        fprintf('ExSEnt: %d channel(s) | m=%g, r=%.3f, lambda=%.4f | parallel=off (text + waitbar)\n', ...
                nchan, m, r, lambda);
    end
end

% ---------------- Compute per channel ----------------
if parallelMode && ~isempty(ver('parallel'))
    parfor ch = 1:nchan
        [HD(ch),HA(ch),HDA(ch),info(ch)] = exsent_channel(data(ch,:), m, r, lambda);
        if showProgress
            fprintf('  ch %3d/%3d: M=%4d | HD=%.6f HA=%.6f HDA=%.6f\n', ...
                    ch, nchan, info(ch).M, HD(ch), HA(ch), HDA(ch));
        end
    end
else
    useWB = showProgress && usejava('desktop');
    hWB = [];
    if useWB
        try hWB = waitbar(0,'Computing ExSEnt...','Name','compute_ExSEnt'); catch, hWB = []; end
    end
    for ch = 1:nchan
        [HD(ch),HA(ch),HDA(ch),info(ch)] = exsent_channel(data(ch,:), m, r, lambda);
        if showProgress
            fprintf('  ch %3d/%3d: M=%4d | HD=%.6f HA=%.6f HDA=%.6f\n', ...
                    ch, nchan, info(ch).M, HD(ch), HA(ch), HDA(ch));
            if ~isempty(hWB) && isvalid(hWB)
                try, waitbar(ch/nchan, hWB, sprintf('Computing ExSEnt... (%d/%d)', ch, nchan)); end
            end
        end
    end
    if ~isempty(hWB) && isvalid(hWB), try, close(hWB); end, end
end

% ---------------- Plot ----------------
switch plotMode

    % ---- (1) Single-channel diagnostic (lambda tuning) ------------------
    case 'validation'
        x = data(1,:);
        [D_vals, A_vals, segment_idx] = extract_DA(x, lambda);

        figure('Name','ExSEnt Validation','Color','w');
        subplot(2,1,1);
        plot(x,'k'); hold on;
        plot(segment_idx, x(segment_idx),'ro','MarkerFaceColor','r');
        xlabel('Samples'); ylabel('Amplitude');
        title('Signal with detected extrema');
        grid on;

        subplot(2,1,2);
        yyaxis left;  plot(D_vals,'b-o','LineWidth',1.2); ylabel('Durations (samples)');
        yyaxis right; plot(A_vals,'r-s','LineWidth',1.2); ylabel('Amplitudes (Δ)');
        xlabel('Segment index');
        title('Duration and amplitude sequences (D and A)');
        grid on;

    % ---- (2) Bifurcation-style summary across all channels --------------
    case 'bifurcation'

        if isempty(chanLabels)
            chanLabels = arrayfun(@(i) sprintf('Ch%d',i), 1:nchan, 'UniformOutput', false);
        end

        % Collect (channel, signal value, HDA) — one row per sample
        temp_bif = cell(1, nchan);
        for ch = 1:nchan
            sig = data(ch,:);
            temp_bif{ch} = [repmat(ch, numel(sig), 1), sig(:), repmat(HDA(ch), numel(sig), 1)];
        end
        bif_data = vertcat(temp_bif{:});
        bif_ch   = bif_data(:,1);
        bif_x    = bif_data(:,2);
        bif_HDA  = bif_data(:,3);

        % Normalise H_DA for colourmap
        Hnorm_scatter = (bif_HDA - min(HDA)) ./ (max(HDA) - min(HDA) + eps);

        figure('Color','w','Name','ExSEnt Bifurcation');

        % Top: scatter coloured by H_DA
        ax1 = axes('Position', [0.10, 0.35, 0.85, 0.60]);
        scatter(ax1, bif_ch, bif_x, 4, Hnorm_scatter, 'filled');
        colormap(ax1, jet);
        cl = colorbar(ax1);
        cl.Label.String = 'H_{DA} (normalised)';
        set(ax1, 'XTick', 1:nchan, 'XTickLabel', [], ...
                 'XLim', [0.5, nchan+0.5], 'FontSize', 22, 'LineWidth', 1.5);
        ylabel(ax1, 'Amplitude');
        title(ax1, sprintf('ExSEnt  |  m=%d, r=%.2f, \\lambda=%.3f', m, r, lambda), ...
            'FontWeight', 'bold', 'FontSize', 13);

        % Bottom: HD, HA, HDA lines per channel
        ax2 = axes('Position', [0.10, 0.10, 0.85, 0.20]);
        hold(ax2, 'on');
        plot(ax2, 1:nchan, HD,  'Color', [0.10, 0.35, 0.65], 'LineWidth', 1.7);
        plot(ax2, 1:nchan, HA,  'Color', [0.85, 0.20, 0.25], 'LineWidth', 1.7);
        plot(ax2, 1:nchan, HDA, 'Color', [0.65, 0.20, 0.75], 'LineWidth', 1.7);
        hold(ax2, 'off');
        legend(ax2, {'H_D','H_A','H_{DA}'}, 'Location', 'best');
        ylabel(ax2, 'ExSEnt');
        xlabel(ax2, 'Channel');
        set(ax2, 'XTick', 1:nchan, 'XTickLabel', chanLabels, ...
                 'XTickLabelRotation', 90, 'TickLabelInterpreter', 'none', ...
                 'XLim', [0.5, nchan+0.5], 'FontSize', 16, 'LineWidth', 1.5);

        linkaxes([ax1, ax2], 'x');

    case 'none'
        % no plot

end
end

% ========================= Channel-level helper ========================= %
function [HD,HA,HDA,info] = exsent_channel(x, m, r, lambda)
% Segment signal into durations/amplitudes
[D_vals, A_vals, segment_idx] = extract_DA(x, lambda);
M = numel(D_vals);

% Clean non-finite / zero-variance edge cases
D_vals = D_vals(isfinite(D_vals));
A_vals = A_vals(isfinite(A_vals));

if M <= m+1 || numel(D_vals) <= m+1 || numel(A_vals) <= m+1 ...
   || std(D_vals) == 0 || std(A_vals) == 0
    [HD,HA,HDA] = deal(NaN);
    info = struct('M',M,'rD',NaN,'rA',NaN,'rDA',NaN, ...
                  'rangeD',range_safe(D_vals),'rangeA',range_safe(A_vals), ...
                  'segment_idx',segment_idx);
    return;
end

% Tolerances (Azami-style on raw sequences)
rD = r * std(D_vals);
rA = r * std(A_vals);

% Delegate to the shared SampEn engine. ZScore=false so rD/rA are applied as
% absolute Chebyshev tolerances on the raw duration/amplitude sequences.
HD = compute_SampEn(D_vals, 'm', m, 'r', rD, 'tau', 1, ...
    'ZScore', false, 'Parallel', false, 'Progress', false);
HA = compute_SampEn(A_vals, 'm', m, 'r', rA, 'tau', 1, ...
    'ZScore', false, 'Parallel', false, 'Progress', false);

% Normalize then interleave [D_norm, A_norm] → 1D series, embed with 2*m
D_norm = normalize(D_vals);
A_norm = normalize(A_vals);
joint_data = reshape([D_norm(:) A_norm(:)]', [], 1);

% Joint r after z-scoring (std≈1): use r directly (ZScore=false, absolute r)
HDA = compute_SampEn(joint_data, 'm', 2*m, 'r', r, 'tau', 1, ...
    'ZScore', false, 'Parallel', false, 'Progress', false);

% Diagnostics
info = struct('M',M,'rD',rD,'rA',rA,'rDA',r, ...
              'rangeD',range_safe(D_vals),'rangeA',range_safe(A_vals), ...
              'segment_idx',segment_idx);
end

% ========================= Segmentation helper ========================== %
function [D_vals, A_vals, segment_idx] = extract_DA(signal, lambda)
% Extract segment durations and net amplitudes via robust sign-change rule.
if nargin < 2, lambda = 0.001; end
signal = signal(:)';                 % row

dx = diff(signal);
threshold = lambda * iqr(dx);

D_vals = [];
A_vals = [];
segment_idx = [];

start_idx = 1;
for i = 2:length(dx)
    if sign(dx(i)) ~= sign(dx(i-1)) && abs(dx(i) - dx(i-1)) > threshold
        duration  = i - start_idx;
        amplitude = signal(i) - signal(start_idx);

        D_vals(end+1)      = duration;   %#ok<AGROW>
        A_vals(end+1)      = amplitude;  %#ok<AGROW>
        segment_idx(end+1) = i;          %#ok<AGROW>

        start_idx = i;
    end
end
end

% ========================= misc ========================= %
function r = range_safe(x)
x = x(isfinite(x));
if isempty(x), r = NaN; else, r = max(x) - min(x); end
end