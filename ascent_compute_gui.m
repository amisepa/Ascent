function [measType, chanlist, tau, m, coarseType, nScales, extraParams, n, vis, parallelComp, progressTrack] = ascent_compute_gui(EEG)
% ascent_compute_gui - EEGLAB GUIs for ascent_compute parameter entry.
%
% Called internally by ascent_get_params. Not intended for direct use.
%
% Spawns 1 or 2 GUIs depending on the selected measure:
%   GUI #1 (always)     : measure, channels, vis/progress/parallel
%   GUI #2a (entropy)   : tau, m  (SampEn, ExSEnt, FracDim, HigFracDim)
%   GUI #2b (MS/fuzzy)  : tau, m, coarse-graining, scale factors, fuzzy power
%   GUI #2c (Aperiodic) : PSD settings, aperiodic fitting, correction toggle
%
% X / Cancel on any GUI aborts and returns measType = [].
%
% OUTPUTS:
%   measType      - selected measure string
%   chanlist      - cell array of channel labels (or all channels)
%   tau           - time lag
%   m             - embedding dimension
%   coarseType    - coarse-graining method string (multiscale only, else [])
%   nScales       - number of scale factors      (multiscale only, else [])
%   extraParams   - struct of aperiodic settings (Aperiodic only, else [])
%   n             - fuzzy power                  (fuzzy measures only, else [])
%   vis           - logical, visualize outputs
%   parallelComp  - logical, use parallel computing
%   progressTrack - logical, track progress in console

% Default outputs
chanlist = []; tau = []; m = [];
coarseType = []; nScales = []; extraParams = []; n = [];
vis = true; parallelComp = true; progressTrack = true;

if nargin < 1 || isempty(EEG)
    error('EEG input is required.');
end

% ---------------------------
% GUI #1 — all measures
% ---------------------------
meas = {'SampEn' 'FuzzEn' 'ExSEnt' 'FracDim' 'HigFracDim' 'Aperiodic' ...
        'MSE' 'mMSE' 'MFE' 'RCMFE' 'RCmvMFE'};

setappdata(0, 'ascent_gui1_cancel', 0);
oldCloseFcn = get(0, 'DefaultFigureCloseRequestFcn');
set(0, 'DefaultFigureCloseRequestFcn', "setappdata(0,'ascent_gui1_cancel',1); delete(gcbf);");

uigeom = { [0.5 0.5] [0.5 0.35 0.15] [0.34 0.33 0.33] };

uilist = {
    {'style' 'text'      'string' 'Measure to compute:'        'fontweight' 'bold'}   % 1
    {'style' 'popupmenu' 'string' meas 'value' 2}                                      % 2  default: FuzzEn

    {'style' 'text'      'string' 'M/EEG channels selection:'  'fontweight' 'bold'}   % 3
    {'style' 'edit'      'string' ''}                                                  % 4
    {'style' 'pushbutton' 'string' '...' 'enable' 'on' ...
     'callback' "tmpEEG=get(gcbf,'userdata'); tmpchanlocs=tmpEEG.chanlocs; " + ...
                "[tmp,tmpval]=pop_chansel({tmpchanlocs.labels},'withindex','on'); " + ...
                "set(findobj(gcbf,'style','edit'),'string',tmpval); " + ...
                "clear tmp tmpEEG tmpchanlocs tmpval" }                                % 5

    {'style' 'checkbox'  'string' 'Visualize outputs'  'value' 1 'fontweight' 'bold'} % 6
    {'style' 'checkbox'  'string' 'Track progress'     'value' 1 'fontweight' 'bold'} % 7
    {'style' 'checkbox'  'string' 'Parallel computing' 'value' 1 'fontweight' 'bold'} % 8
    };

param = inputgui(uigeom, uilist, 'pophelp(''ascent_compute'')', 'Ascent EEGLAB plugin', EEG);

set(0, 'DefaultFigureCloseRequestFcn', oldCloseFcn);

if isempty(param)
    measType = [];
    return;
end

% Parse GUI #1
measType = meas{param{1}};

if ~isempty(param{2})
    chanlist = cellstr(split(string(param{2})));
else
    chanlist = {EEG.chanlocs.labels}';
end

vis           = logical(param{3});
progressTrack = logical(param{4});
parallelComp  = logical(param{5});

% ---------------------------
% GUI #2b — multiscale / fuzzy measures (includes tau and m)
% ---------------------------
isMS = contains(lower(measType), {'mse','mmse','mfe','cmfe','rcmfe','rcmvmfe'});
isFZ = contains(lower(measType), {'fuzzen','mfe','rcmfe','rcmvmfe'});

if isMS || isFZ

    if isMS, enMS = 'on'; else, enMS = 'off'; end
    if isFZ, enFZ = 'on'; else, enFZ = 'off'; end

    cTypes = {'Mean' 'Median' 'Std' 'Variance'};

    uigeom2 = { [0.5 0.5] [1] [0.5 0.5] [1] [0.5 0.5] [1] [0.5 0.5] [1] [0.5 0.5] };

    uilist2 = {
        {'style' 'text' 'string' 'Time lag (tau):'          'fontweight' 'bold'}  % 1
        {'style' 'edit' 'string' '1'}                                              % 2
        {}                                                                         % 3
        {'style' 'text' 'string' 'Embedding dimension (m):' 'fontweight' 'bold'}  % 4
        {'style' 'edit' 'string' '2'}                                              % 5
        {}                                                                         % 6
        {'style' 'text'      'string' 'Coarse graining method:' 'fontweight' 'bold' 'enable' enMS}   % 7
        {'style' 'popupmenu' 'string' cTypes 'value' 1                             'enable' enMS}    % 8
        {}                                                                         % 9
        {'style' 'text'      'string' 'Number of scale factors:' 'fontweight' 'bold' 'enable' enMS}  % 10
        {'style' 'edit'      'string' '30'                                         'enable' enMS}    % 11
        {}                                                                         % 12
        {'style' 'text'      'string' 'Fuzzy power:'             'fontweight' 'bold' 'enable' enFZ}  % 13
        {'style' 'edit'      'string' '2'                                          'enable' enFZ}    % 14
        };

    param2 = inputgui(uigeom2, uilist2, 'pophelp(''ascent_compute'')', 'Ascent EEGLAB plugin', EEG);

    if isempty(param2)
        measType = [];
        return;
    end

    % Fixed-position parse (inputgui returns all fields regardless of enable state)
    tau = str2double(param2{1});
    m   = str2double(param2{2});

    if isMS
        coarseType = cTypes{param2{3}};
        if contains(lower(coarseType), 'var'), coarseType = 'var'; end
        nScales = str2double(param2{4});
    else
        coarseType = []; nScales = [];
    end

    if isFZ
        n = str2double(param2{5});
    else
        n = [];
    end

    return;
end

% ---------------------------
% GUI #2c — Aperiodic parameters (no tau/m needed)
% ---------------------------
if strcmpi(measType, 'Aperiodic')

    apModes  = {'Fixed (default)' 'Knee'};

    uigeomAP = { [1] ...
                 [0.55 0.45] [0.55 0.45] ...
                 [1] [1] ...
                 [0.55 0.45] [0.55 0.45] [0.55 0.45] [0.55 0.45] [0.55 0.45] [0.55 0.45] ...
                 [1] [1] [1] };

    uilistAP = {
        % --- PSD settings ---
        {'style' 'text' 'string' '— PSD settings —' 'fontweight' 'bold' 'fontsize' 10}

        {'style' 'text' 'string' 'Frequency range (Hz):'}
        {'style' 'edit' 'string' '1 40'}

        {'style' 'text' 'string' 'Window length (s):'}
        {'style' 'edit' 'string' '4'}

        % --- Aperiodic fitting ---
        {}
        {'style' 'text' 'string' '— Aperiodic fitting —' 'fontweight' 'bold' 'fontsize' 10}

        {'style' 'text' 'string' 'Aperiodic mode:'}
        {'style' 'popupmenu' 'string' apModes 'value' 1}

        {'style' 'text' 'string' 'Fitting freq range (Hz):'}
        {'style' 'edit' 'string' '1 40'}

        {'style' 'text' 'string' 'Max peaks:'}
        {'style' 'edit' 'string' '6'}

        {'style' 'text' 'string' 'Min peak height:'}
        {'style' 'edit' 'string' '0.05'}

        {'style' 'text' 'string' 'Peak threshold:'}
        {'style' 'edit' 'string' '2.0'}

        {'style' 'text' 'string' 'Peak width limits (Hz):'}
        {'style' 'edit' 'string' '0.5 12'}

        % --- Correction ---
        {}
        {'style' 'text'     'string' '— Aperiodic correction —' 'fontweight' 'bold' 'fontsize' 10}
        {'style' 'checkbox' 'string' 'Subtract aperiodic component from PSD' ...
                            'value' 1 'fontweight' 'bold'}
        };

    paramAP = inputgui(uigeomAP, uilistAP, 'pophelp(''ascent_compute'')', ...
                       'Ascent – Aperiodic settings', EEG);

    if isempty(paramAP)
        measType = [];
        return;
    end

    % Parse — positions match uilistAP control order above
    psdFreqRange = str2num(paramAP{1});                         %#ok<ST2NM>
    winSec       = str2double(paramAP{2});
    apModeStr    = strtrim(strsplit(lower(apModes{paramAP{3}}), ' (')); apModeStr = apModeStr{1};
    fitFreqRange = str2num(paramAP{4});                         %#ok<ST2NM>
    maxPeaks     = str2double(paramAP{5});
    minPeakHt    = str2double(paramAP{6});
    peakThresh   = str2double(paramAP{7});
    peakWidths   = str2num(paramAP{8});                         %#ok<ST2NM>
    correctAP    = logical(paramAP{9});

    % Pack into struct for ascent_get_params
    extraParams.freqRange        = psdFreqRange;
    extraParams.winSec           = winSec;
    extraParams.overlap          = 0.5;        % fixed default; tune via command line
    extraParams.window           = 'hann';     % fixed default; tune via command line
    extraParams.aperiodicMode    = apModeStr;
    extraParams.fitFreqRange     = fitFreqRange;
    extraParams.maxPeaks         = maxPeaks;
    extraParams.minPeakHeight    = minPeakHt;
    extraParams.peakThreshold    = peakThresh;
    extraParams.peakWidthLimits  = peakWidths;
    extraParams.correctAperiodic = correctAP;

    coarseType = []; nScales = []; n = [];
    return;
end

% ---------------------------
% GUI #2a — all other entropy measures (SampEn, ExSEnt, FracDim, HigFracDim)
% just tau and m
% ---------------------------
uigeom_tm = { [0.5 0.5] [1] [0.5 0.5] };

uilist_tm = {
    {'style' 'text' 'string' 'Time lag (tau):'          'fontweight' 'bold'}
    {'style' 'edit' 'string' '1'}
    {}
    {'style' 'text' 'string' 'Embedding dimension (m):' 'fontweight' 'bold'}
    {'style' 'edit' 'string' '2'}
    };

param_tm = inputgui(uigeom_tm, uilist_tm, 'pophelp(''ascent_compute'')', 'Ascent EEGLAB plugin', EEG);

if isempty(param_tm)
    measType = [];
    return;
end

tau = str2double(param_tm{1});
m   = str2double(param_tm{2});

coarseType = []; nScales = []; n = []; extraParams = [];

end