function ascent_gui_data(fig, action)
% ascent_gui_data - Header and Data row of the first ASCENT window (ascent_compute_gui).
%
% Draws the header (logo, name, version), loads a dataset or the sample dataset
% shipped with the plugin into the window, and shows which file is loaded. The
% dataset is kept in the figure's 'userdata', which inputgui returns on OK.
%
%   ascent_gui_data(fig, 'init')    % wider window, header band, dataset description
%   ascent_gui_data(fig, 'load')    % pop_loadset file dialog
%   ascent_gui_data(fig, 'sample')  % ascent_sample_data.set
%   ascent_gui_data(fig, 'show')    % refresh the description only
%
% Copyright (C) - ASCENT EEGLAB PLUGIN - Cedric Cannard, 2026

root = fileparts(which('eegplugin_ascent.m'));
E = [];
switch action
    case 'init'
        band = make_room(fig, 125);
        draw_header(fig, root, band);
    case 'load'
        try
            E = pop_loadset();
        catch ME
            errordlg(ME.message, 'ASCENT'); return
        end
    case 'sample'
        E = pop_loadset('filename', 'ascent_sample_data.set', 'filepath', root);
        % No file name: ascent_compute saves its outputs over the file a dataset
        % was loaded from, and the sample file must stay as shipped
        E.filename = ''; E.filepath = '';
        E.setname = 'ASCENT sample dataset';
        E.etc.ascent_sample = true;
end
if has_data(E)
    set(fig, 'userdata', E);
    fprintf('ASCENT: loaded %s\n', source_name(E));
end

% Two lines: which file, then what it holds
EEG = get(fig, 'userdata');
if has_data(EEG)
    if isempty(EEG.icaweights), ica = 'no ICA'; else, ica = sprintf('%d ICs', size(EEG.icaweights,1)); end
    if EEG.trials > 1, ep = sprintf(', %d epochs', EEG.trials); else, ep = ''; end
    str = {['Loaded: ' source_name(EEG)], ...
           sprintf('%d channels, %g Hz, %.1f min%s, %s', EEG.nbchan, EEG.srate, ...
                   EEG.pnts*EEG.trials/EEG.srate/60, ep, ica)};
    col = [0 0.45 0];  wt = 'bold';
    tip = fullfile(EEG.filepath, EEG.filename);   % full path when the name is cut off
else
    str = {'No dataset loaded.', 'Load one, or the sample dataset.'};
    col = [0.35 0.35 0.35];  wt = 'normal';  tip = '';
end
set(findobj(fig, 'tag', 'ascent_data'), 'string', str, 'foregroundcolor', col, ...
    'fontweight', wt, 'tooltipstring', tip);
end

function s = source_name(EEG)
% File name, or what the dataset is when it has none
if isfield(EEG, 'etc') && isfield(EEG.etc, 'ascent_sample')
    s = 'ascent_sample_data.set (sample, not saved)';
elseif ~isempty(EEG.filename)
    s = EEG.filename;
elseif ~isempty(EEG.setname)
    s = [strtrim(EEG.setname) ' (not saved)'];
else
    s = 'unsaved dataset';
end
end

function band = make_room(fig, H)
% Widen the window (supergui controls are in normalized units, so they stretch)
% and add a band of H pixels at the top for the header: the controls are held
% in pixels while the window grows upward, then go back to normalized units.
% Returns the band [x y w h] in pixels.
kids = allchild(fig);
kids = kids(isprop(kids, 'Units'));
set(fig, 'Units', 'pixels');
p = get(fig, 'Position');
W = max(p(3), 700);
x = max(10, p(1) - (W - p(3))/2);
set(fig, 'Position', [x p(2) W p(4)]);
set(kids, 'Units', 'pixels');
set(fig, 'Position', [x max(40, p(2) - H) W p(4) + H]);
set(kids, 'Units', 'normalized');
band = [0 p(4) W H];
end

function draw_header(fig, root, band)
% Logo, name and version in the band at the top of the window
if ~isempty(findobj(fig, 'tag', 'ascent_logo')), return; end
bg = get(fig, 'color');
fg = [0 0 0.4];
try icadefs; fg = GUITEXTCOLOR; catch, end %#ok<NODEF>
sz = band(4) - 20;  y = band(2) + 8;  objs = [];
try
    % circular logo with transparency, blended with the window colour
    [img, ~, alpha] = imread(fullfile(root, 'ascent_logo.png'));
    a = double(alpha) / 255;
    img = double(img) / 255 .* a + reshape(bg, 1, 1, 3) .* (1 - a);
    ax = axes('Parent', fig, 'Units', 'pixels', 'Position', [20 y sz sz], 'Tag', 'ascent_logo');
    hi = image(ax, img); axis(ax, 'image', 'off');
    try set(hi, 'Interpolation', 'bilinear'); catch, end   % R2023b+
    objs = ax;
catch
    sz = 0;
end
x = 20 + sz + 18;  w = band(3) - x - 20;
objs(end+1) = uicontrol(fig, 'Style', 'text', 'Units', 'pixels', 'Position', [x y+0.52*sz w 0.40*sz], ...
    'String', 'ASCENT', 'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', ...
    'BackgroundColor', bg, 'ForegroundColor', fg);
objs(end+1) = uicontrol(fig, 'Style', 'text', 'Units', 'pixels', 'Position', [x y+0.04*sz w 0.46*sz], ...
    'String', {'Aperiodic Signal Complexity Estimation', ...
               sprintf('for Neurophysiological Time series  (v%s)', eegplugin_ascent)}, ...
    'FontSize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', bg, 'ForegroundColor', fg);
set(objs, 'Units', 'normalized');   % follow the window if it is resized
end

function ok = has_data(E)
ok = ~isempty(E) && isstruct(E) && isfield(E, 'data') && ~isempty(E.data);
end
