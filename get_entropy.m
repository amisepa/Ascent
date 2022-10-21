%% Computes entropy on EEGLAB-formatted data.
%
% Cedric Cannard, August 2022

function [ae, se, rcmfe, f] = get_entropy(EEG,varargin)

ae = [];
se = [];
mse = [];
rcmfe = [];

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
% if ~isfield(EEG.chanlocs, 'X') || isempty(EEG.chanlocs(1).X), error("Electrode locations are required. " + ...
%         "Go to 'Edit > Channel locations' and import the appropriate coordinates for your montage"); end
if isempty(EEG.ref)
    warning(['EEG data not referenced! Referencing is highly recommended ' ...
        '(e.g., CSD-transformation, infinity-, or average- reference)!']); 
end

% Continuous/epoched data
if length(size(EEG.data)) == 2
    continuous = true;
else
%     continuous = false;
%     EEG = epoch2continuous(EEG);  %NOT CORRECT!! need to find another solution
end

%% GUI
if nargin < 2
    drawnow;
    eTypes = {'Approximate entropy' 'Sample entropy' 'Multiscale entropy' 'Refined composite multiscale fuzzy entropy (default)'};
    cTypes = {'Mean' 'Standard deviation (default)' 'Variance'};
    mTypes = {'1' '2 (default)' '3'};
    tTypes = {'0.5' '1 (default)' '2'};
    uigeom = { [.5 .6 .2] .5 [.5 .8] .5 [.5 .8] .5 [.5 .2] .5 [.5 .2] .5 .5};
    uilist = {
        {'style' 'text' 'string' 'Channel selection:'} ...
        {'style' 'edit' 'tag' 'chanlist'} ...
        {'style' 'pushbutton' 'string'  '...', 'enable' 'on' ...
        'callback' "tmpEEG = get(gcbf, 'userdata'); tmpchanlocs = tmpEEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels},'withindex','on'); set(findobj(gcbf,'tag','chanlist'),'string',tmpval); clear tmp tmpEEG tmpchanlocs tmpval" } ...
        {} ...
        {'style' 'text' 'string' 'Entropy type:'} ...
        {'style' 'popupmenu' 'string' eTypes 'tag' 'etype' 'value' 3} ...
        {} ...
        {'style' 'text' 'string' 'Coarse graining method:'} ...
        {'style' 'popupmenu' 'string' cTypes 'tag' 'stype' 'value' 2} ...
        {} ...
        {'style' 'text' 'string' 'Tau (time lag):'} ...
        {'style' 'popupmenu' 'tag' 'tau' 'string' tTypes  'value' 2}  ...
        {} ...
        {'style' 'text' 'string' 'Embedding dimension:'} ...
        {'style' 'popupmenu' 'tag' 'tau' 'string' mTypes  'value' 2}  ...
        {} ...
        {'style' 'checkbox' 'string' 'Bandpass filter each scale to control for spectral bias (recommended)?','tag' 'filter','value',1}  ...
        };
    userchoices = inputgui(uigeom,uilist,'pophelp(''pop_entropy'')','entropy EEGLAB plugin',EEG);
end

% Decode inputs
if nargin < 2   % GUI
    opt = struct;
    opt.chanlist = split(userchoices{1});
    opt.entropyType = eTypes(userchoices{2});
    % if ~isempty(userchoices{3})
    opt.coarseType = cTypes(userchoices{3});
    % end

    tau  = tTypes{userchoices{4}};
    if contains(tau, 'default'), tau = tau(1); end
    opt.tau = str2double(tau);

    m  = tTypes{userchoices{5}};
    if contains(m, 'default'),  m = m(1); end
    opt.m = str2double(m);

    opt.filter = logical(userchoices{6});

else % command line inputs
    opt = varargin;
end

%% Defaults
if ~isfield(opt, 'chanlist') || isempty(opt.chanlist)
    disp('No channels were selected: selecting all available channels.')
    chanlist = {EEG.chanlocs.labels}';
else
    chanlist = opt.chanlist;
end
if ~isfield(opt, 'entropyType') || isempty(opt.entropyType)
    disp('No entropy type selected: selecting Refined composite multiscale fuzzy entropy (default).')
    entropyType = 'Refined composite multiscale fuzzy entropy (default)';
else
    entropyType = opt.entropyType{:};
end
if ~isfield(opt, 'coarseType') || isempty(opt.coarseType)
    disp('No coarse graining method selected: selecting standard deviation (default).')
    coarseType = 'Standard deviation';
else
    coarseType = opt.coarseType{:};
end
if ~isfield(opt, 'tau') || isempty(opt.tau)
    disp('No time lag selected: selecting tau = 1 (default).')
    tau = 1;
else
    tau = opt.tau;
end
if ~isfield(opt, 'm') || isempty(opt.m)
    disp('No embedding dimension selected: selecting 2 (default).')
    m = 2;
else
    m = opt.m;
end
if ~isfield(opt, 'filter') || isempty(opt.filter)
    disp('Selecting bandpass filtering at each time scale to control for the spectral bias (default).')
    filtData = true;
else
    filtData = opt.filter;
end

% hardcode r because data are scaled to have SD of 1
r = .15;    % threshold (default = .15 of standard deviation)
nchan = length(chanlist);

switch entropyType

    case 'Approximate entropy'
        ae = nan(nchan,1);
        for ichan = 1:nchan
            ae(ichan,:) = approx_entropy(EEG.data(ichan,:), m, r);
            fprintf('   Channel %s: %s \n', EEG.chanlocs(ichan).labels, ae(ichan,:))
        end

    case 'Sample entropy'
        se = nan(nchan,1);
        % t1 = tic;
        if continuous && EEG.pnts <= 34000  % CHECK THRESHOLD
            disp('Computing sample entropy on continuous data (standard method)...')
            for ichan = 1:nchan
                se(ichan,:) = sample_entropy(EEG.data(ichan,:),m,r,tau);  % standard method
                fprintf('   Channel %d: %s \n', ichan, se(ichan,:))
            end

        else
            disp('Large continuous data detected, computing sample entropy using the fast method...')
            for ichan = 1:nchan
                se(ichan,:) = sample_entropy_fast(EEG.data(ichan,:),m,r); % fast method
                fprintf('   Channel %d: %s \n', ichan, se(ichan,:))
            end
        end
        % t2 = toc(t1)

    case 'Refined composite multiscale fuzzy entropy (default)'

        % number of scale factors to compute (starting at 2)
        nScales = inpdlg('Select fuzzy power: ');
        if isempty(nScales)
            nScales = 30;
        end

        fuzzypower = 2;

        % rcmfe_NN = RCMFE_std(data, 2, .15, 2, 1, nscalesNN);
        [rcmfe, f] = get_rcmfe(data, m, r, fuzzypower, tau, nScales, EEG.srate);

% [1] H. Azami and J. Escudero, "Refined Multiscale Fuzzy Entropy based on Standard Deviation 
% for Biomedical Signal Analysis", Medical & Biological Engineering & Computing, 2016.


end

%% PLot results

% open convert_3Dto2D for plotting
