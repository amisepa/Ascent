function [measType, chanlist, tau, m, coarseType, nScales, filtData, n, vis, parallelComp, progressTrack] = ascent_compute_gui(EEG)
% ascent_compute_gui
% Two GUIs returning compute parameters.
% - GUI #1 includes visualize, track progress, parallel computing (defaults ON)
% - GUI #2 includes coarse/nScales and fuzzy power, greyed out when irrelevant
% - X close behaves like Cancel

% Defaults
measType = []; chanlist = []; tau = []; m = [];
coarseType = []; nScales = []; filtData = []; n = [];
vis = true; parallelComp = true; progressTrack = true;

if nargin < 1 || isempty(EEG)
    error('EEG input is required.');
end

% ---------------------------
% GUI #1
% ---------------------------
meas = {'SampEn' 'FuzzEn' 'ExSEnt' 'FracDim' 'MSE' 'mMSE' 'MFE' 'RCMFE (default)' 'RCmvMFE'};

setappdata(0,'ascent_gui1_cancel',0);
oldCloseFcn = get(0,'DefaultFigureCloseRequestFcn');
set(0,'DefaultFigureCloseRequestFcn', "setappdata(0,'ascent_gui1_cancel',1); delete(gcbf);");

% Boxes: 2 + 3 + 2 + 2 + 3 = 12
uigeom = { ...
    [0.5 0.5] ...
    [0.5 0.35 0.15] ...
    [0.5 0.5] ...
    [0.5 0.5] ...
    [0.34 0.33 0.33] ...
    };

uilist = {
    {'style' 'text' 'string' 'Measure to compute:' 'fontweight' 'bold'}               % 1
    {'style' 'popupmenu' 'string' meas 'value' 8}                                      % 2

    {'style' 'text' 'string' 'M/EEG channels selection:' 'fontweight' 'bold'}          % 3
    {'style' 'edit' 'string' ''}                                                       % 4
    {'style' 'pushbutton' 'string' '...' 'enable' 'on' ...
     'callback' "tmpEEG=get(gcbf,'userdata'); tmpchanlocs=tmpEEG.chanlocs; " + ...
                "[tmp,tmpval]=pop_chansel({tmpchanlocs.labels},'withindex','on'); " + ...
                "set(findobj(gcbf,'style','edit'),'string',tmpval); " + ...
                "clear tmp tmpEEG tmpchanlocs tmpval" }                                % 5

    {'style' 'text' 'string' 'Time lag (tau):' 'fontweight' 'bold'}                    % 6
    {'style' 'edit' 'string' '1'}                                                      % 7

    {'style' 'text' 'string' 'Embedding dimension (m):' 'fontweight' 'bold'}           % 8
    {'style' 'edit' 'string' '2'}                                                      % 9

    {'style' 'checkbox' 'string' 'Visualize outputs'  'value' 1 'fontweight' 'bold'}   % 10
    {'style' 'checkbox' 'string' 'Track progress'     'value' 1 'fontweight' 'bold'}   % 11
    {'style' 'checkbox' 'string' 'Parallel computing' 'value' 1 'fontweight' 'bold'}   % 12
    };

param = inputgui(uigeom, uilist, 'pophelp(''ascent_compute'')', 'Ascent EEGLAB plugin', EEG);

set(0,'DefaultFigureCloseRequestFcn', oldCloseFcn);

if getappdata(0,'ascent_gui1_cancel') == 1 || isempty(param)
    measType = [];
    return;
end

% Parse GUI #1 using fixed indices
measType = meas{param{1}};
if contains(measType,'default')
    measType = extractBefore(measType,' (default)');
end

if ~isempty(param{2})
    chanlist = split(string(param{2}));
    chanlist = cellstr(chanlist(:));
else
    chanlist = {EEG.chanlocs.labels}';
end

tau = str2double(param{3});
m   = str2double(param{4});

vis          = logical(param{5});
progressTrack= logical(param{6});
parallelComp = logical(param{7});

% ---------------------------
% GUI #2
% ---------------------------
isMS = contains(lower(measType), {'mse','mmse','mfe','rcmfe','rcmvmfe'});
isFZ = contains(lower(measType), {'fuzzen','mfe','rcmfe','rcmvmfe'});

if isMS, enMS = 'on'; else, enMS = 'off'; end
if isFZ, enFZ = 'on'; else, enFZ = 'off'; end

cTypes = {'Mean' 'Median' 'Std' 'Variance (default)'};

setappdata(0,'ascent_gui2_cancel',0);
oldCloseFcn2 = get(0,'DefaultFigureCloseRequestFcn');
set(0,'DefaultFigureCloseRequestFcn', "setappdata(0,'ascent_gui2_cancel',1); delete(gcbf);");

% Boxes: 2 + 1 + 2 + 1 + 2 = 8
uigeom2 = { ...
    [0.5 0.5] ...
    [1] ...
    [0.5 0.5] ...
    [1] ...
    [0.5 0.5] ...
    };

uilist2 = {
    {'style' 'text' 'string' 'Coarse graining method:' 'enable' enMS}
    {'style' 'popupmenu' 'string' cTypes 'value' 4 'enable' enMS}
    {}
    {'style' 'text' 'string' 'Number of scale factors:' 'enable' enMS}
    {'style' 'edit' 'string' '30' 'enable' enMS}
    {}
    {'style' 'text' 'string' 'Fuzzy power:' 'enable' enFZ}
    {'style' 'edit' 'string' '2' 'enable' enFZ}
    };

param2 = inputgui(uigeom2, uilist2, 'pophelp(''ascent_compute'')', 'Ascent EEGLAB plugin', EEG);

set(0,'DefaultFigureCloseRequestFcn', oldCloseFcn2);

if getappdata(0,'ascent_gui2_cancel') == 1 || isempty(param2)
    measType = [];
    return;
end

% Parse GUI #2
p = 1;
if isMS
    coarseType = cTypes{param2{p}}; p = p + 1;
    if contains(lower(coarseType),'var')
        coarseType = 'var';
    end
    nScales = str2double(param2{p}); p = p + 1;
else
    coarseType = [];
    nScales = [];
end

if isFZ
    n = str2double(param2{p});
else
    n = [];
end

% Filter removed throughout
filtData = [];

end
