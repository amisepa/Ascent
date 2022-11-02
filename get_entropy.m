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
%   nScales - number of scale factors (default = 15)
%   filtData - apply band pass filters to each time scale to control for
%   broadband spectral bias (see Kosciessa et al. 2020 for more detail). 
%
% Cedric Cannard, August 2022

function [ae,se,fe,p,mse,rcmfe,scales] = get_entropy(EEG, entropyType, chanlist, tau, m, coarseType, nScales, filtData, n, vis)

ae = [];
se = [];
fe = [];
p = [];
mse = [];
rcmfe = [];
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
    continuous = false; %%%%%%%%%%%%% ADD OPTION TO RUN ON EPOCHED DATA %%%%%%%%%
end

%% 1st GUI to select channels and type of entropy

if nargin == 1
    eTypes = {'Approximate entropy' 'Sample entropy' 'Fuzzy entropy' 'Multiscale entropy' 'Multiscale fuzzy entropy' 'Refined composite multiscale fuzzy entropy (default)'};
    uigeom = { [.5 .9] .5 [.5 .4 .2] .5 [.5 .1] .5 [.5 .1] };
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
        {'style' 'edit' 'string' '15' 'tag' 'n'}  ...
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
if contains(lower(entropyType), 'multiscale')
    if ~exist('coarseType','var') || isempty(coarseType)
        disp('No coarse graining method selected: selecting standard deviation (default).')
        coarseType = 'Standard deviation';
    end
    if ~exist('nScales','var') || isempty(nScales)
        disp('Number of scale factors not selected: selecting nScales = 15 (default).')
        nScales = 15;
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
        disp('Computing approximate entropy...')
        ae = nan(nchan,1);
        for ichan = 1:nchan
            ae(ichan,:) = compute_ae(EEG.data(ichan,:), m, r);
            fprintf('   %s: %6.3f \n', EEG.chanlocs(ichan).labels, ae(ichan,:))
        end

    case 'SE'
        disp('Computing sample entropy...')
        se = nan(nchan,1);
        if continuous && EEG.pnts <= 34000  % CHECK THRESHOLD
            disp('Computing sample entropy on continuous data (standard method)...')
            disp('If this takes too long, try the fast method (see main_script code)')
            for ichan = 1:nchan
                se(ichan,:) = compute_se(EEG.data(ichan,:),m,r,tau);  % standard method
                fprintf('   %s: %6.3f \n', EEG.chanlocs(ichan).labels, se(ichan,:))
            end

        else
            disp('Large continuous data detected, computing sample entropy using the fast method...')
            for ichan = 1:nchan
                se(ichan,:) = compute_se_fast(EEG.data(ichan,:),m,r); % fast method
                fprintf('   %s: %6.3f \n', EEG.chanlocs(ichan).labels, se(ichan,:))
            end
        end
    
    case 'FE'
        disp('Computing fuzzy entropy...')
        fe = nan(nchan,1);
        for ichan = 1:nchan
            [fe(ichan,:), p(ichan,:)] = compute_fe(EEG.data(ichan,:),m,r,n,tau);
            fprintf('   %s: %6.3f \n', EEG.chanlocs(ichan).labels, fe(ichan,:))
        end

    case 'MSE'
        disp('Computing multiscale entropy...')
        mse = nan(nchan,nScales);
        scales = nan(2,nScales);
%         if vis
%             figure('color','w')
%         end
        for ichan = 1:nchan
            [mse(ichan,:), scales] = compute_mse(EEG.data(ichan,:),m,r,tau,coarseType,nScales,filtData,EEG.srate);
%             if vis
%                 plot(1:nScales, mse, 'linewidth',2); hold on;
%                 nans = isnan(mse(1,:));
%                 if nans(1)
%                     xticks(2:nScales); xticklabels(join(string(scales(:,2:end)),1)); xtickangle(45)
%                     xlim([2 nScales]);
%                 else
%                     xticks(1:nScales); xticklabels(join(string(scales),1)); xtickangle(45)
%                 end
%                 legend({EEG.chanlocs.labels})
%             end
        end

    case 'RFCMFE'

%         % number of scale factors to compute (starting at 2)
%         nScales = inpdlg('Select fuzzy power: ');
%         if isempty(nScales)
%             nScales = 30;
%         end
% 
%         % rcmfe_NN = RCMFE_std(data, 2, .15, 2, 1, nscalesNN);
%         [rcmfe, scales] = get_rcmfe(data, m, r, fuzzypower, tau, nScales, EEG.srate);


end

% Cite references
disp('Please cite: ')
switch coarseType
    case 'Mean'
        disp('Costa, Goldberger, and Peng (2002) - Multiscale entropy analysis of complex physiologic time series. Physical review letters.')
    case 'SD'
        disp('   [1] Costa, Goldberger, and Peng (2002) - Multiscale entropy analysis of complex physiologic time series. Physical review letters.')
        disp('   [2] Azami and Escudero (2016) - Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis. Medical & Biological Engineering & Computing')
    case 'Variance'
        disp('   [1] Costa, Goldberger, and Peng (2002) - Multiscale entropy analysis of complex physiologic time series. Physical review letters.')
        disp('   [2] Azami and Escudero (2016) - Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis. Medical & Biological Engineering & Computing')
end
if filtData
    disp('Bandpass filters were applied to each scale factor to control for spectral bias, following recommendations by: ');
    disp('   [3] Kosciessa, Kloosterman, and Garrett (2020) - Standard multiscale entropy reflects neural dynamics at mismatched temporal scales: What''s signal irregularity got to do with it? Plos Computational Biology.')
end
