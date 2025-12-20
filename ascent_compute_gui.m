function [measType, chanlist, tau, m, coarseType, nScales, filtData, n, vis, parallelComp] = ascent_compute_gui(EEG)
% Minimal wrapper that produces 3 GUIs and returns values.

measType = []; chanlist = []; tau = []; m = [];
coarseType = []; nScales = []; filtData = []; n = [];
vis = true; parallelComp = false;

% --- GUI #1: type + channels + tau + m + plot
meas = {'SampEn' 'FuzzEn' 'ExSEnt' 'FracDim' 'MSE' 'mMSE' 'MFE' 'RCMFE (default)' 'RCmvMFE'};
uigeom = { [.5 .9] .5 [.5 .4 .2] .5 [.5 .1] .5 [.5 .1] .5 .5};
uilist = {
    {'style' 'text' 'string' 'Measure to compute:' 'fontweight' 'bold'}
    {'style' 'popupmenu' 'string' meas 'tag' 'etype' 'value' 8}
    {}
    {'style' 'text' 'string' 'M/EEG channels selection:' 'fontweight' 'bold'}
    {'style' 'edit' 'tag' 'chanlist'}
    {'style' 'pushbutton' 'string'  '...', 'enable' 'on' ...
     'callback' "tmpEEG = get(gcbf, 'userdata'); tmpchanlocs = tmpEEG.chanlocs; [tmp,tmpval] = pop_chansel({tmpchanlocs.labels},'withindex','on'); set(findobj(gcbf,'tag','chanlist'),'string',tmpval); clear tmp tmpEEG tmpchanlocs tmpval" }
    {}
    {'style' 'text' 'string' 'Time lag (tau):' 'fontweight' 'bold'}
    {'style' 'edit' 'string' '1' 'tag' 'tau'}
    {}
    {'style' 'text' 'string' 'Embedding dimension (m):' 'fontweight' 'bold'}
    {'style' 'edit' 'string' '2' 'tag' 'm'}
    {}
    {'style' 'checkbox' 'string' 'Visualize outputs' 'tag' 'vis' 'value' 1 'fontweight' 'bold'}
    };
param = inputgui(uigeom, uilist, 'pophelp(''ascent_compute'')','Ascent EEGLAB plugin', EEG);
if isempty(param), return; end
measType = meas{param{1}};
if ~isempty(param{2})
    chanlist = split(param{2});
else
    chanlist = {EEG.chanlocs.labels}';
end
tau  = str2double(param{3});
m    = str2double(param{4});
vis  = logical(param{5});

% --- GUI #2: coarse + nScales + fuzzy power (greyed out when not relevant)
isMS = contains(lower(measType), {'mse','mmse','mfe','rcmfe','rcmvmfe'});
isFZ = contains(lower(measType), {'fuzzen','mfe','rcmfe','rcmvmfe'});
if isMS, enMS = 'on'; else, enMS = 'off'; end
if isFZ, enFZ = 'on'; else, enFZ = 'off'; end

cTypes = {'Mean' 'Median' 'Std' 'Variance (default)'};

% Geometry:
% rows = 5
% boxes = 2 + 1 + 2 + 1 + 2 = 8
uigeom = { ...
    [0.5 0.5] ...  % coarse label + popup
    [1] ...        % spacer
    [0.5 0.5] ...  % nScales label + edit
    [1] ...        % spacer
    [0.5 0.5] ...  % fuzzy label + edit
    };

uilist = {
    {'style' 'text' 'string' 'Coarse graining method:' 'tag' 'coarseLbl' 'enable' enMS}
    {'style' 'popupmenu' 'string' cTypes 'tag' 'stype' 'value' 4 'enable' enMS}

    {}

    {'style' 'text' 'string' 'Number of scale factors:' 'tag' 'nScalesLbl' 'enable' enMS}
    {'style' 'edit' 'string' '30' 'tag' 'nScales' 'enable' enMS}

    {}

    {'style' 'text' 'string' 'Fuzzy power:' 'tag' 'fuzzyLbl' 'enable' enFZ}
    {'style' 'edit' 'string' '2' 'tag' 'fuzzyN' 'enable' enFZ}
    };

param = inputgui( ...
    uigeom, uilist, ...
    'pophelp(''ascent_compute'')', ...
    'Ascent EEGLAB plugin', EEG);

if isempty(param)
    return;
end

% --- Parse outputs (order matters)
p = 1;
if isMS
    coarseType = cTypes{param{p}}; p = p + 1;
    if contains(lower(coarseType),'var')
        coarseType = 'var';
    end
    nScales = str2double(param{p}); p = p + 1;
else
    coarseType = [];
    nScales    = [];
end

if isFZ
    n = str2double(param{p});
else
    n = [];
end

end
