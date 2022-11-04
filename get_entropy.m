%% EEGLAB plugin to compute entropy-based measures on MEEG data.
% Also works on other types of biosignals (e.g., ECG, HRV).
% Founded on code developed by Costa et al. (2002) and Azami and Escudero
% (2017). 
% 
% INPUTS:
%   EEG - EEG structure in EEGLAB format
%   entropyType - 'Approximate entropy', 'Sample entropy', 'Fuzzy entropy',
%                   'Multiscale entropy', 'Multiscale fuzzy entropy',
%                   'Refined composite multiscale fuzzy entropy (default)'
%   chanlist - EEG channels of interest (empty will select all channels)
%   tau - time lag (default = 1)
%   m - embedding dimension (default = 2)
%   coarseType - coarse graining method for multiscale entropies: 'Mean', 
%                 'Standard deviation' (default), or 'Variance'
%   nScales - number of scale factors (default = 30)
%   filtData - apply band pass filters to each time scale to control for
%       broadband spectral bias (see Kosciessa et al. 2020 for more detail). 
%   vis - visualize entropy outputs (1, default) or not (0)
%
% Copyright - Cedric Cannard, 2022

function [entropy,scales] = get_entropy(EEG, entropyType, chanlist, tau, m, coarseType, nScales, filtData, n, vis)

entropy = [];
scales = [];

% add path to subfolders
mainpath = fileparts(which('get_entropy.m'));
addpath(fullfile(mainpath, 'functions'));

% Basic checks and warnings
if nargin < 1
    help pop_entropy; return; 
end
if isempty(EEG.data)
    error('Empty dataset.'); 
end
if isempty(EEG.chanlocs(1).labels)
    error('No channel labels.'); 
end
if vis
    if ~isfield(EEG.chanlocs, 'X') || isempty(EEG.chanlocs(1).X) 
        error("Electrode locations are required. " + ...
            "Go to 'Edit > Channel locations' and import the appropriate coordinates for your montage"); 
    end
end
if isempty(EEG.ref)
    warning(['EEG data not referenced! Referencing is highly recommended ' ...
        '(e.g., average- reference)!']); 
end

% Continuous/epoched data
if length(size(EEG.data)) == 2
    continuous = true;
else
    continuous = false; %%%%%%%%%%%%% ADD OPTION TO RUN ON EPOCHED DATA %%%%%%%%%
end

%% 1st GUI to select channels and type of entropy

if nargin == 1
    eTypes = {'Approximate entropy' 'Sample entropy' 'Fuzzy entropy' 'Multiscale entropy' 'Multiscale fuzzy entropy' 'Refined composite multiscale fuzzy entropy (default)'};
    uigeom = { [.5 .9] .5 [.5 .4 .2] .5 [.5 .1] .5 [.5 .1] .5 .5};
    uilist = {
        {'style' 'text' 'string' 'Entropy type:' 'fontweight' 'bold'} ...
        {'style' 'popupmenu' 'string' eTypes 'tag' 'etype' 'value' 6} ...
        {} ...
        {'style' 'text' 'string' 'Channel selection:' 'fontweight' 'bold'} ...
        {'style' 'edit' 'tag' 'chanlist'} ...
        {'style' 'pushbutton' 'string'  '...', 'enable' 'on' ...
        'callback' "tmpEEG = get(gcbf, 'userdata'); tmpchanlocs = tmpEEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels},'withindex','on'); set(findobj(gcbf,'tag','chanlist'),'string',tmpval); clear tmp tmpEEG tmpchanlocs tmpval" } ...
        {} ...
        {'style' 'text' 'string' 'Tau (time lag):' 'fontweight' 'bold'} ...
        {'style' 'edit' 'string' '1' 'tag' 'tau'} ...
        {} ...
        {'style' 'text' 'string' 'Embedding dimension:' 'fontweight' 'bold'} ...
        {'style' 'edit' 'string' '2' 'tag' 'm'}  ...
        {} ...
        {'style' 'checkbox' 'string' 'Plot entropy outputs' 'tag' 'vis' 'value' 1 'fontweight' 'bold'}  ...
            };
    param = inputgui(uigeom,uilist,'pophelp(''pop_entropy'')','entropy EEGLAB plugin',EEG);
    entropyType = eTypes{param{1}};
    if ~isempty(param{2})
        chanlist = split(param{2});
    else
        chanlist = {EEG.chanlocs.labels}';
    end
    tau  = str2double(param{3});
    m = str2double(param{4});
    vis = logical(param{5});
end


%% 2nd GUI to select additional parameters

if nargin == 1 && contains(entropyType, 'Multiscale')
    cTypes = {'Mean' 'Standard deviation (default)' 'Variance'};
    uigeom = { [.5 .6] .5 [.9 .3] .5 .5 };
    uilist = {
        {'style' 'text' 'string' 'Coarse graining method:'} ...
        {'style' 'popupmenu' 'string' cTypes 'tag' 'stype' 'value' 2} ...
        {} ...
        {'style' 'text' 'string' 'Number of scale factors:' } ...
        {'style' 'edit' 'string' '30' 'tag' 'n'}  ...
        {} ...
        {'style' 'checkbox' 'string' 'Bandpass filter each time scale (recommended to control for spectral bias)','tag' 'filtData','value',0}  ...
            };
    param = inputgui(uigeom,uilist,'pophelp(''pop_entropy'')','get_entropy() EEGLAB plugin',EEG);
    coarseType = cTypes{param{1}};
    nScales = str2double(param{2});
    filtData = logical(param{3});
end

%% 3rd GUI for fuzzy power

if nargin == 1 && contains(entropyType, 'Fuzzy')
    uigeom = { [.9 .3] };
    uilist = { {'style' 'text' 'string' 'Fuzzy power:' } ...
        {'style' 'edit' 'string' '2' 'tag' 'n'}  };
    param = inputgui(uigeom,uilist,'pophelp(''pop_entropy'')','entropy EEGLAB plugin',EEG);
    n = str2double(param{1});
end

%% Defaults if something was missed in command line

if ~exist('chanlist','var') || isempty(chanlist)
    disp('No channels were selected: selecting all channels (default)')
    chanlist = {EEG.chanlocs.labels}';
end
if ~exist('entropyType','var') || isempty(entropyType)
    disp('No entropy type selected: selecting Refined composite multiscale fuzzy entropy (default)')
    entropyType = 'Refined composite multiscale fuzzy entropy (default)';
end
if ~exist('tau','var') || isempty(tau)
    disp('No time lag selected: selecting tau = 1 (default).')
    tau = 1;
end
if ~exist('m','var') || isempty(m)
    disp('No embedding dimension selected: selecting m = 2 (default).')
    m = 2;
end
if ~exist('vis','var') || isempty(vis)
    disp('Plotting option not selected: turning plotting ON (default).')
    vis = true;
end
if contains(lower(entropyType), 'multiscale')
    if ~exist('coarseType','var') || isempty(coarseType)
        disp('No coarse graining method selected: selecting standard deviation (default).')
        coarseType = 'Standard deviation';
    end
    if ~exist('nScales','var') || isempty(nScales)
        disp('Number of scale factors not selected: selecting nScales = 15 (default).')
        nScales = 30;
    end
    if ~exist('filtData','var') || isempty(filtData)
%         disp('Selecting bandpass filtering at each time scale to control for the spectral bias (default).')
        filtData = false;
    end
end
if contains(lower(entropyType), 'fuzzy')
    if ~exist('n','var') || isempty(n)
        disp('No fuzzy power selected: selecting n = 2 (default).')
        n = 2;
    end
end

% Simplify entropy names 
if contains(lower(entropyType), 'approximate')
    entropyType = 'AE';
elseif contains(lower(entropyType), 'sample')
    entropyType = 'SE';
elseif strcmpi(entropyType, 'fuzzy entropy')
    entropyType = 'FE';
elseif strcmpi(entropyType, 'multiscale entropy')
    entropyType = 'MSE';
elseif strcmpi(entropyType, 'multiscale fuzzy entropy')
    entropyType = 'MFE';
elseif contains(lower(entropyType), 'refined')
    entropyType = 'RCMFE';
end

%% Compute entropy depending on choices

% Hardcode r to .15 because data are scaled to have SD of 1
r = .15;
nchan = length(chanlist);

switch entropyType

    case 'AE'
        disp('Computing approximate entropy (AE)...')
        entropy = nan(nchan,1);
        parfor ichan = 1:nchan
            entropy(ichan,:) = compute_ae(EEG.data(ichan,:), m, r);
            fprintf('   %s: %6.3f \n', EEG.chanlocs(ichan).labels, entropy(ichan,:))
        end
        if vis, plot_entropy(entropy, EEG.chanlocs); end

    case 'SE'
        disp('Computing sample entropy (SE)...')
        entropy = nan(nchan,1);
        if continuous && EEG.pnts <= 34000  % CHECK THRESHOLD
            disp('Computing sample entropy on continuous data (standard method)...')
            disp('If this takes too long, try the fast method (see main_script code)')
            parfor ichan = 1:nchan
                entropy(ichan,:) = compute_se(EEG.data(ichan,:),m,r,tau);  % standard method
                fprintf('   %s: %6.3f \n', EEG.chanlocs(ichan).labels, entropy(ichan,:))
            end

        else
            disp('Large continuous data detected, computing sample entropy using the fast method...')
            parfor ichan = 1:nchan
                entropy(ichan,:) = compute_se_fast(EEG.data(ichan,:),m,r); % fast method
                fprintf('   %s: %6.3f \n', EEG.chanlocs(ichan).labels, entropy(ichan,:))
            end
        end
        if vis, plot_entropy(entropy, EEG.chanlocs); end

    case 'FE'
        disp('Computing fuzzy entropy (FE)...')
        entropy = nan(nchan,1);
        parfor ichan = 1:nchan
            [entropy(ichan,:), p(ichan,:)] = compute_fe(EEG.data(ichan,:),m,r,n,tau);
            fprintf('   %s: %6.3f \n', EEG.chanlocs(ichan).labels, entropy(ichan,:))
        end        
        if vis, plot_entropy(entropy, EEG.chanlocs); end

    case 'MSE'
        disp('Computing multiscale entropy (MSE)...')
        parfor ichan = 1:nchan
            fprintf('Channel %d: \n', ichan)
            [entropy(ichan,:), scales] = compute_mse(EEG.data(ichan,:),m,r,tau,coarseType,nScales,filtData,EEG.srate);
        end
        if vis, plot_entropy(entropy, EEG.chanlocs); end

    case 'MFE'
        disp('Computing multiscale fuzzy entropy (MFE)...')
        parfor ichan = 1:nchan
            fprintf('Channel %d: \n', ichan)
            [mfe(ichan,:), scales] = compute_mfe(EEG.data(ichan,:),m,r,tau,coarseType,nScales,filtData,EEG.srate,n);
        end
        if vis, plot_entropy(entropy, EEG.chanlocs); end

    case 'RFCMFE'

%         disp('Computing refined composite multiscale fuzzy entropy (RCMFE)...')
%         for ichan = 1:nchan
%             fprintf('Channel %d: \n', ichan)
%             [rcmfe(ichan,:), scales] = compute_rcmfe(EEG.data(ichan,:),m,r,tau,coarseType,nScales,filtData,EEG.srate);

end
