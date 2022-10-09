%% Computes entropy on EEGLAB-formatted data.
%
% Cedric Cannard, August 2022

function se = pop_entropy(EEG,varargin)

% Basic checks and warnings
if nargin < 1, help pop_entropy; return; end
if isempty(EEG.data), error('Cannot process empty dataset.'); end
if isempty(EEG.chanlocs(1).labels), error('Cannot process without channel labels.'); end
% if ~isfield(EEG.chanlocs, 'X') || isempty(EEG.chanlocs(1).X), error("Electrode locations are required. " + ...
%         "Go to 'Edit > Channel locations' and import the appropriate coordinates for your montage"); end
if isempty(EEG.ref), warning(['EEG data not referenced! Referencing is highly recommended ' ...
        '(e.g., CSD-transformation, infinity-, or average- reference)!']); end
if length(size(EEG.data)) == 2
    continuous = true;
else
%     continuous = false;
%     EEG = epoch2continuous(EEG);  %NOT CORRECT!! need to find another solution
end

% GUI
if nargin < 2
    drawnow;
    eTypes = {'Sample entropy' 'Multiscale entropy' 'Refined composite multiscale fuzzy entropy (default)'};
    cTypes = {'Mean' 'Standard deviation (default)' 'Variance'};
    uigeom = { [.5 .5 .5] .5 [.5 .5 .5] .5 [.5 .5 .5] .5 .5};
    uilist = {
        {'style' 'text' 'string' 'Channel selection:'}, ...
        {'style' 'edit' 'string' ' ' 'tag' 'chanlist'}, ...
        {'style' 'pushbutton' 'string'  '...', 'enable' 'on' ...
        'callback' "tmpEEG = get(gcbf, 'userdata'); tmpchanlocs = tmpEEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels},'withindex','on'); set(findobj(gcbf,'tag','chanlist'),'string',tmpval); clear tmp tmpEEG tmpchanlocs tmpval" }, ...
        {} ...
        {'style' 'text' 'string' 'Entropy type:'} ...
        {'style' 'popupmenu' 'string' eTypes 'tag' 'etype'} {} ...
        {} ...
        {'style' 'text' 'string' 'Coarse graining method:'} ...
        {'style' 'popupmenu' 'string' cTypes 'tag' 'stype'} {} ...
        {} ...
        {'style' 'checkbox' 'string' 'Bandpass filter each scale to control for spectral bias (recommended)?','tag' 'filter','value',0}  ...
        };
    result = inputgui(uigeom,uilist,'pophelp(''pop_entropy'')','entropy EEGLAB plugin',EEG);
   
    if isempty(result), return; end
    
    %decode user inputs
    args = struct;
    if ~isempty(result{1})
        chanlist = split(result{1})';
        args.chanlist = chanlist';
    end
    args.etype = eTypes(result{2});
    args.ctype = cTypes(result{3});
    args.filter = logical(result{4});
else
    args = varargin;
end

% Defaults
% args = struct(args{:});
if ~isfield(args, 'chanlist') || isempty(args.chanlist)
    args.chanlist = {EEG.chanlocs.labels}';
end
if ~isfield(args, 'etype') || isempty(args.etype)
    args.etype = 'Refined composite multiscale fuzzy entropy (default)';
end
if ~isfield(args, 'ctype') || isempty(args.ctype)
    args.ctype = 'Standard deviation (default)';
end
if ~isfield(args, 'filter') || isempty(args.filter)
    args.filter = false;
end

tau = 1;    % time lag (default = 1)
dim = 2;      % embedding dimension (default = 2)
r = .15;    % threshold (default = .15 of standard deviation)

% Sample entropy
if strcmp(args.etype, 'Sample entropy')
%     t1 = tic;
    if continuous && EEG.pnts > 34000
        disp('Long continuous data detected --> computing sample entropy using fast method.')
        for iChan = 1:length([args.chanlist])
            disp([' Channel ' num2str(iChan)])
%             se2(iChan,:) = fastSampen(EEG.data(iChan,:),m,r);       %fast method
            se(iChan,:) = sampEnFast(EEG.data(iChan,:),dim,r);
        end
    elseif continuous && EEG.pnts <= 34000
        disp('Computing standard sample entropy on continuous data. This may take a while...')
        for iChan = 1:length([args.chanlist])
            disp([' Channel ' num2str(iChan)])
            se(iChan,:) = sampEn(EEG.data(iChan,:),dim,r,tau);     %standard method
        end
    end
% t2 = toc(t1)
% open convert_3Dto2D for plotting
end

fuzzypower = 2;      % fuzzy power
nscales = 10; % number of scale factors to compute (starting at 2)
% rcmfe_NN = RCMFE_std(data, 2, .15, 2, 1, nscalesNN);
[rcmfe, f] = get_rcmfe(data,dim,r,fuzzypower,tau,nscales,EEG.srate);
    
% [1] H. Azami and J. Escudero, "Refined Multiscale Fuzzy Entropy based on Standard Deviation 
% for Biomedical Signal Analysis", Medical & Biological Engineering & Computing, 2016.
