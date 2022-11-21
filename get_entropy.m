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
% USAGE:
%   1) Import EEG datra into EEGLAB
%   2) Preprocess as necessary (e.g., reference, clean data with ASR, etc.)
%   3) GUI mode: Tools > Compute entropy
% or
%   EEG = get_entropy(EEG);     % launch GUI mode
% or
%   EEG = get_entropy(EEG, 'Fuzzy entropy','Fpz',[],[],'Variance'); 
%                                       % compute fuzzy entropy only on Fpz 
%                                       % channel with default tau and m but using
%                                       % variance for the coarse-graining
% or
%   EEG = get_entropy(EEG,'Multiscale fuzzy entropy',[],[],[],[],50,1,[],0);
%                                       % compute multiscale fuzzy entropy
%                                       % on all channels with default
%                                       % parameters but on 50 time scales,
%                                       % controlling for spectral bias, and
%                                       % turning plotting OFF.
% 
%
% Copyright - Cedric Cannard, 2022

function [EEG, com] = get_entropy(EEG, entropyType, chanlist, tau, m, coarseType, nScales, filtData, n, vis)

com = '';

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
% if exist(vis,'var') && ~isempty(vis)
if ~isfield(EEG.chanlocs, 'X') || isempty(EEG.chanlocs(1).X) 
    error("Electrode locations are required. " + ...
        "Go to 'Edit > Channel locations' and import the appropriate coordinates for your montage"); 
end
% end
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
    eTypes = {'Approximate entropy' 'Sample entropy' 'Fuzzy entropy' 'Multiscale entropy' 'Multiscale fuzzy entropy' 'Refined composite multiscale fuzzy entropy'};
    uigeom = { [.5 .9] .5 [.5 .4 .2] .5 [.5 .1] .5 [.5 .1] .5 .5};
    uilist = {
        {'style' 'text' 'string' 'Entropy type:' 'fontweight' 'bold'} ...
        {'style' 'popupmenu' 'string' eTypes 'tag' 'etype' 'value' 5} ...
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
    entropyType = 'Refined composite multiscale fuzzy entropy';
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
        disp('Number of scale factors not selected: selecting nScales = 30 (default).')
        nScales = 30;
    end
    if ~exist('filtData','var') || isempty(filtData)
        filtData = false;
    end
else
    coarseType = []; 
    nScales = [];
    filtData = [];
end
if contains(lower(entropyType), 'fuzzy')
    if ~exist('n','var') || isempty(n)
        disp('No fuzzy power selected: selecting n = 2 (default).')
        n = 2;
    end
else
    n = [];
end

%% Compute entropy depending on choices

% Hardcode r to .15 because data are scaled to have SD of 1
r = .15;

% index with channels of interest
nchan = length(chanlist);
if nchan > 1 && nchan < EEG.nbchan
    [~, chanIdx] = intersect({EEG.chanlocs.labels}, split(chanlist));
else
    chanIdx = 1:EEG.nbchan;
end

% preallocate memory for the entropy variable
if contains(lower(entropyType), 'multiscale')
    entropy = nan(nchan, nScales);
else
    entropy = nan(nchan,1);
end

% COMPUTE ENTROPY 
switch entropyType

    case 'Approximate entropy'
        disp('Computing approximate entropy...')
        entropy = nan(nchan,1);
        progressbar('Channels')
        for ichan = 1:nchan
%             fprintf('Channel %d \n', ichan)
            entropy(ichan,:) = compute_ae(EEG.data(chanIdx(ichan),:), m, r);
            fprintf('   %s: %6.3f \n', EEG.chanlocs(ichan).labels, entropy(ichan,:))
            progressbar(ichan/nchan)
        end
        if vis, plot_entropy(entropy, EEG.chanlocs, chanIdx); end

    case 'Sample entropy'
        disp('Computing sample entropy...')
        entropy = nan(nchan,1);
        if continuous && EEG.pnts <= 34000  % CHECK THRESHOLD
            disp('Computing sample entropy on continuous data (standard method)...')
            disp('If this takes too long, try the fast method (see main_script code)')
            progressbar('Channels')
            for ichan = 1:nchan
%                 fprintf('Channel %d \n', ichan)
                entropy(ichan,:) = compute_se(EEG.data(chanIdx(ichan),:),m,r,tau);  % standard method
                fprintf('   %s: %6.3f \n', EEG.chanlocs(chanIdx(ichan)).labels, entropy(ichan,:))
                progressbar(ichan/nchan)
            end
            
            % plot
            if vis, plot_entropy(entropy, EEG.chanlocs, chanIdx); end

        else
            disp('Large continuous data detected, computing sample entropy using the fast method...')
            progressbar('Channels')
            for ichan = 1:nchan
%                 fprintf('Channel %d \n', ichan)
                entropy(ichan,:) = compute_se_fast(EEG.data(chanIdx(ichan),:),m,r); % fast method
                fprintf('   %s: %6.3f \n', EEG.chanlocs(chanIdx(ichan)).labels, entropy(ichan,:))
                progressbar(ichan/nchan)
            end
        end

        % plot
        if vis, plot_entropy(entropy, EEG.chanlocs, chanIdx); end

    case 'Fuzzy entropy'
        disp('Computing fuzzy entropy...')
        entropy = nan(nchan,1);
        progressbar('Channels')
        for ichan = 1:nchan
%             fprintf('Channel %d \n', ichan)
            entropy(ichan,:) = compute_fe(EEG.data(chanIdx(ichan),:),m,r,n,tau);
            fprintf('   %s: %6.3f \n', EEG.chanlocs(chanIdx(ichan)).labels, entropy(ichan,:))
            progressbar(ichan/nchan)
        end    

        % Plot
        if vis, plot_entropy(entropy, EEG.chanlocs, chanIdx); end

    case 'Multiscale entropy'
        disp('Computing multiscale entropy...')
        progressbar('Channels')
        for ichan = 1:nchan
            fprintf('Channel %d: \n', ichan)
            [entropy(ichan,:), scales] = compute_mse(EEG.data(chanIdx(ichan),:), ...
                m,r,tau,coarseType,nScales,filtData,EEG.srate);
            progressbar(ichan/nchan)
        end
        
        % Remove NaN scales
        idx = isnan(entropy(1,:));
        entropy(:,idx) = []; scales(idx) = [];
        
        % Plot
        if vis, plot_entropy(entropy, EEG.chanlocs, chanIdx); end

    case 'Multiscale fuzzy entropy'
        disp('Computing multiscale fuzzy entropy...')
        progressbar('Channels')
        t1 = tic;
        for ichan = 1:nchan
            fprintf('Channel %d: \n', ichan)
            [entropy(ichan,:), scales] = compute_mfe(EEG.data(chanIdx(ichan),:),m,r, ...
                tau,coarseType,nScales,filtData,EEG.srate,n);
            progressbar(ichan/nchan)
        end
        toc(t1)

        % Remove NaN scales
        idx = isnan(entropy(1,:));
        entropy(:,idx) = []; scales(idx) = [];
        
        % Plot
        if vis, plot_entropy(entropy, EEG.chanlocs, chanIdx); end

    case 'Refined composite multiscale fuzzy entropy'

%      [rcmfe(ichan,:), scales] = compute_rcmfe(EEG.data(ichan,:),m,r,tau,coarseType,nScales,filtData,EEG.srate);

    otherwise
        error('Unknown entropy type. Please select one of the options (see help get_entropy).')
end

% save outputs in EEG structure
EEG.entropy = entropy;
if contains(lower(entropyType), 'multiscale')
    EEG.scales = scales;
end

% Command history
chanLabels = strjoin(chanlist);
chanLabels = insertBefore(chanLabels," ", "'");
chanLabels = insertAfter(chanLabels," ", "'");

% TO FIX   
com = sprintf('EEG = get_entropy(''%s'', {''%s''}, %d, %d, %s, %d, %d, %s, %d);', ...
    entropyType,chanLabels,tau,m,coarseType,nScales,filtData,'[]',vis);

gong
disp('Done!')

end

