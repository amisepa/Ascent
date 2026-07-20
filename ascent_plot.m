function ascent_plot(entropyData, chanlocs, entropyType, scales, varargin)
% ascent_plot  Visualize Ascent measures on channel or IC data.
%
% Dispatcher only: parses inputs and routes to the figure that matches the
% measure. The aperiodic and multiscale figures share no layout, so they live
% in separate files:
%   ascent_plot_aperiodic   - entropyType containing 'aperiodic'
%   ascent_plot_multiscale  - everything else (heatmap/cluster, topo, curve, scalar)
%
% Neither of those defines any analysis: band power and the 1/f model come from
% compute_AperiodicBandPower, the same function ascent_compute uses, so a figure
% always shows what was computed and saved rather than its own re-derivation.
%
% entropyData : [sig x scale] or [sig x scale x time]; rows are channels, or ICs
%               when 'ICA' is true
% chanlocs    : EEGLAB channel locations. In ICA mode these are the scalp
%               channels the decomposition ran on (EEG.chanlocs(EEG.icachansind)),
%               i.e. the ones icawinv back-projects onto.
% entropyType : label, e.g., 'RCMFEsigma' or 'MSE (std)'
% scales      : cellstr of scale labels or numeric 1:S
% varargin    : context-dependent optional inputs:
%   Name-value pairs (can appear anywhere in varargin):
%       'ICA',           true/false
%       'icawinv',       EEG.icawinv  [nChans x nICs]
%       'ClusterThresh', scalar in (0,1)  percentile threshold (default 0.75)
%       'dipfit',        EEG.dipfit   (ICA only; adds the region/dipole figure)
%       'alphaPower',    struct('raw',[sig x 1],'osc',[sig x 1],'band',[f1 f2])
%                        Precomputed by ascent_compute. Aperiodic only. When
%                        omitted it is derived with compute_AperiodicBandPower.
%   Aperiodic positional (entropyType contains 'aperiodic'):
%       varargin{1} = offset        [sig x 1] or [sig x nTimes]
%       varargin{2} = freqs         [1 x freq]
%       varargin{3} = psd           [sig x freq] or [sig x freq x nTimes]
%       varargin{4} = psd_corrected [sig x freq] or [sig x freq x nTimes] or []
%       varargin{5} = times         [1 x nTimes] (time-varying only)
%   Time-resolved positional:
%       varargin{1} = time_sec      [1 x time]
%
% Copyright (C) - ASCENT EEGLAB PLUGIN - Cedric Cannard, 2021-2025

% Strip named pairs before positional parsing
opt.isICA         = false;
opt.icawinv       = [];
opt.clusterThresh = 0.75;
opt.dipfit        = [];
opt.alphaPower    = [];
rmIdx             = [];
for vi = 1:(numel(varargin)-1)
    if ischar(varargin{vi}) || isstring(varargin{vi})
        switch lower(char(varargin{vi}))
            case 'ica'
                opt.isICA = logical(varargin{vi+1});
                rmIdx = [rmIdx vi vi+1]; %#ok<AGROW>
            case 'icawinv'
                opt.icawinv = varargin{vi+1};
                rmIdx       = [rmIdx vi vi+1]; %#ok<AGROW>
            case 'clusterthresh'
                opt.clusterThresh = varargin{vi+1};
                rmIdx             = [rmIdx vi vi+1]; %#ok<AGROW>
            case 'dipfit'
                opt.dipfit = varargin{vi+1};
                rmIdx      = [rmIdx vi vi+1]; %#ok<AGROW>
            case 'alphapower'
                opt.alphaPower = varargin{vi+1};
                rmIdx          = [rmIdx vi vi+1]; %#ok<AGROW>
        end
    end
end
varargin(rmIdx) = [];

isAperiodic = contains(lower(entropyType),'aperiodic');

if opt.isICA
    % icawinv is needed for back-projection (non-aperiodic and static aperiodic).
    % Time-varying aperiodic only uses isICA to relabel axes, no back-projection.
    isTimeVaryingAperiodic = isAperiodic && size(entropyData,2) > 1;
    if isempty(opt.icawinv) && ~isTimeVaryingAperiodic
        error('ascent_plot: ICA=true requires ''icawinv'' (EEG.icawinv).');
    end
    if ~isempty(opt.icawinv)
        if size(opt.icawinv,2) ~= size(entropyData,1)
            error('ascent_plot: icawinv has %d ICs but data has %d rows.', ...
                  size(opt.icawinv,2), size(entropyData,1));
        end
        if size(opt.icawinv,1) ~= numel(chanlocs)
            error(['ascent_plot: icawinv maps onto %d channels but %d chanlocs were ' ...
                   'given. In ICA mode pass EEG.chanlocs(EEG.icachansind), not the ' ...
                   'full EEG.chanlocs.'], size(opt.icawinv,1), numel(chanlocs));
        end
    end
end

% Positional varargin
opt.time_sec = []; opt.offset = []; opt.freqs = [];
opt.psd      = []; opt.psd_corrected = []; opt.times = [];

if isAperiodic
    if numel(varargin)>=1, opt.offset        = varargin{1}; end
    if numel(varargin)>=2, opt.freqs         = varargin{2}; end
    if numel(varargin)>=3, opt.psd           = varargin{3}; end
    if numel(varargin)>=4, opt.psd_corrected = varargin{4}; end
    if numel(varargin)>=5, opt.times         = varargin{5}; end
    ascent_plot_aperiodic(entropyData, chanlocs, opt);
else
    if ~isempty(varargin), opt.time_sec = varargin{1}; end
    ascent_plot_multiscale(entropyData, chanlocs, entropyType, scales, opt);
end
end
