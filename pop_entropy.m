%% Computes entropy on EEGLAB-formatted data.
%
% Cedric Cannard, August 2022

function [sampEn, rcmfe, f] = pop_entropy(EEG,varargin)

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
    eTypes = {'Sample entropy' 'Multiscale entropy' 'Refined composite multiscale fuzzy entropy (default)'};
    cTypes = {'Mean' 'Standard deviation (default)' 'Variance'};
    uigeom = { [.5 .6 .2] .5 [.5 .7] .5 [.5 .2] .5 [.5 .7] .5 .5};
    uilist = {
        {'style' 'text' 'string' 'Channel selection:'} ...
        {'style' 'edit' 'tag' 'chanlist'} ...
        {'style' 'pushbutton' 'string'  '...', 'enable' 'on' ...
        'callback' "tmpEEG = get(gcbf, 'userdata'); tmpchanlocs = tmpEEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels},'withindex','on'); set(findobj(gcbf,'tag','chanlist'),'string',tmpval); clear tmp tmpEEG tmpchanlocs tmpval" } ...
        {} ...
        {'style' 'text' 'string' 'Entropy type:'} ...
        {'style' 'popupmenu' 'string' eTypes 'tag' 'etype' 'value' 3} ...
        {} ...
        {'style' 'text' 'string' 'Tau (time lag):'} ...
        {'style' 'edit' 'tag' 'tau' 'value' 1}  ...
        {} ...
        {'style' 'text' 'string' 'Coarse graining method:'} ...
        {'style' 'popupmenu' 'string' cTypes 'tag' 'stype' 'value' 2} ...
        {} ...
        {'style' 'checkbox' 'string' 'Bandpass filter each scale to control for spectral bias (recommended)?','tag' 'filter','value',1}  ...
        };
    userchoices = inputgui(uigeom,uilist,'pophelp(''pop_entropy'')','entropy EEGLAB plugin',EEG);
   
    % Decode user inputs
    if ~isempty(userchoices)
        opt = struct;
        if ~isempty(userchoices{1})
            chanlist = split(userchoices{1});
            opt.chanlist = chanlist;
        end
        if ~isempty(userchoices{2})
            opt.etype = eTypes(userchoices{2});
        end
        if ~isempty(userchoices{3})
            opt.ctype = cTypes(userchoices{3});
        end
        if ~isempty(userchoices{4})
            opt.filter = logical(userchoices{4});
        end
    else
        return
    end
else % command line inputs
    opt = varargin;
end

%% Defaults
if ~isfield(opt, 'chanlist') || isempty(opt.chanlist)
    opt.chanlist = {EEG.chanlocs.labels}';
end
if ~isfield(opt, 'etype') || isempty(opt.etype)
    opt.etype = 'Refined composite multiscale fuzzy entropy (RCMFE; default)';
end
if ~isfield(opt, 'ctype') || isempty(opt.ctype)
    opt.ctype = 'Standard deviation (default)';
end
if ~isfield(opt, 'filter') || isempty(opt.filter)
    opt.filter = true;
end

tau = 1;    % time lag (default = 1)
m = 2;      % embedding dimension (default = 2)
r = .15;    % threshold (default = .15 of standard deviation)

% Sample entropy
if strcmp(opt.etype, 'Sample entropy')
%     t1 = tic;
    if continuous && EEG.pnts > 34000
        disp('Long continuous data detected --> computing sample entropy using fast method.')
        for iChan = 1:length([opt.chanlist])
            disp([' Channel ' num2str(iChan)])
            sampEn(iChan,:) = fastSampen(EEG.data(iChan,:),m,r);       % fast method
        end
    elseif continuous && EEG.pnts <= 34000
        disp('Computing standard sample entropy on continuous data. This may take a while...')
        for iChan = 1:length([opt.chanlist])
            disp([' Channel ' num2str(iChan)])
            sampEn(iChan,:) = sampEn(EEG.data(iChan,:),m,r,tau);     % standard method
        end
    end
% t2 = toc(t1)
% open convert_3Dto2D for plotting
end

%% Refined composite multiscale fuzzy entropy
fuzzypower = 2;     % fuzzy power
nscales = 10;       % number of scale factors to compute (starting at 2)

% rcmfe_NN = RCMFE_std(data, 2, .15, 2, 1, nscalesNN);
[rcmfe, f] = get_rcmfe(data,m,r,fuzzypower,tau,nscales,EEG.srate);


% [1] H. Azami and J. Escudero, "Refined Multiscale Fuzzy Entropy based on Standard Deviation 
% for Biomedical Signal Analysis", Medical & Biological Engineering & Computing, 2016.
