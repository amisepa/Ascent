function p = ascent_get_params(EEG, varargin)
% ascent_get_params - Parse inputs and apply defaults for ascent_compute.
%
% Called internally by ascent_compute. Not intended for direct use.
% Handles three input paths:
%   1) GUI mode   : launches ascent_compute_gui, unpacks outputs and
%                   aperiodic GUI struct into p
%   2) Name-value : parses varargin key-value pairs into p
%   3) Defaults   : fills any remaining empty fields with the values below
%
% INPUTS:
%   EEG      - EEGLAB EEG structure (required for GUI and channel defaults)
%   varargin - name-value pairs (empty triggers GUI mode)
%
% OUTPUT:
%   p - parameter struct. Returns [] if the user cancels the GUI.
%
%   General:
%     p.measure       'RCMFE' | 'SampEn' | 'FuzzEn' | 'ExSEnt' | 'FracDim' |
%                     'HigFracDim' | 'Aperiodic' | 'MSE' | 'mMSE' |
%                     'MFE' | 'RCMFE' | 'RCmvMFE'             (default: 'RCMFE')
%     p.domain        'channel' (scalp channels) | 'ica'     (default: 'channel')
%                     'ica' runs the measure on EEG.icaact; the channel
%                     selection is ignored and all ICs are used.
%     p.chanlist       cell array of channel labels           (default: all)
%     p.tau            time lag for embedding                 (default: 1)
%     p.m              embedding dimension                    (default: 2)
%     p.r              similarity bound (fraction of SD)      (default: 0.15)
%     p.vis            plot outputs                           (default: true)
%     p.parallel       use parallel computing                 (default: true)
%     p.progress       print progress to console              (default: true)
%
%   Multiscale (MSE, mMSE, MFE, RCMFE, RCmvMFE):
%     p.coarsing         coarse-graining: 'mean'|'median'|'std'|'var'   (default: 'mean')
%     p.num_scales       number of scale factors                        (default: 30)
%     p.filter_mode      'none' | 'narrowband' (mMSE only)             (default: 'narrowband' if mMSE)
%     p.TimeWin          window length for time-resolved mMSE          (default: [])
%     p.TimeStep         step size for time-resolved mMSE              (default: [])
%     p.TimeOnly         compute only time-resolved mMSE               (default: false)
%
%   Fuzzy (FuzzEn, MFE, RCMFE, RCmvMFE):
%     p.n              fuzzy power                            (default: 2)
%     p.kernel_meth    kernel type: 'exponential'|...         (default: 'exponential')
%     p.blocksize      block size for fuzzy computation       (default: 256)
%     p.fuzzy_mode     'local' | ...                          (default: 'local')
%
%   Aperiodic:
%     p.freqRange      PSD frequency range [fMin fMax] Hz     (default: [1 40])
%     p.winSec         Welch segment length (s)               (default: 4)
%     p.psdOverlap     Welch segment overlap fraction         (default: 0.5)
%     p.windowType     taper window type                      (default: 'hann')
%     p.aperiodicMode  'fixed' | 'knee'                       (default: 'fixed')
%     p.fitFreqRange   specparam fitting range [fMin fMax] Hz (default: p.freqRange)
%     p.maxPeaks       max number of spectral peaks to fit    (default: 6)
%     p.minPeakHeight  min peak height above aperiodic (log)  (default: 0.05)
%     p.peakThreshold  peak detection threshold (SDs)         (default: 2.0)
%     p.peakWidthLimits [min max] peak width in Hz            (default: [1 12])
%     p.correctAperiodic subtract aperiodic model from PSD    (default: true)
%     p.alphaBand      [fMin fMax] Hz for saved alpha power   (default: [8 13])
%     p.timeResolved   sliding-window aperiodic fit           (default: false)

% ----------------------------
% Initialize all fields to []
% ----------------------------
p.measure          = [];
p.domain           = [];
p.chanlist         = [];
p.tau              = [];
p.m                = [];
p.r                = [];
p.vis              = [];
p.parallel         = [];
p.progress         = [];

% Multiscale
p.coarsing         = [];
p.num_scales       = [];
p.filter_mode      = [];
p.TimeWin          = [];
p.TimeStep         = [];
p.TimeOnly         = [];

% Fuzzy
p.n                = [];
p.kernel_meth      = [];
p.blocksize        = [];
p.fuzzy_mode       = [];

% Aperiodic
p.freqRange        = [];
p.winSec           = [];
p.psdOverlap       = [];
p.windowType       = [];
p.aperiodicMode    = [];
p.fitFreqRange     = [];
p.maxPeaks         = [];
p.minPeakHeight    = [];
p.peakThreshold    = [];
p.peakWidthLimits  = [];
p.correctAperiodic = [];
p.alphaBand        = [];
p.timeResolved     = [];
p.slidWinSec       = [];
p.slidOverlap      = [];
p.slidNAvg         = [];
p.slidPeakWidthLimits = [];

% ----------------------------
% GUI or name-value parsing
% ----------------------------
if isempty(varargin) || (numel(varargin) == 1 && isempty(varargin{1}))

    % --- GUI mode ---
    [p.measure, p.chanlist, p.tau, p.m, p.coarsing, p.num_scales, extraParams, ...
     p.n, p.vis, p.parallel, p.progress, p.domain] = ascent_compute_gui(EEG);

    if isempty(p.measure)
        p = [];     % signal abort to caller
        return
    end

    % Unpack aperiodic GUI struct
    if strcmpi(p.measure, 'Aperiodic') && isstruct(extraParams)
        p.freqRange        = extraParams.freqRange;
        p.winSec           = extraParams.winSec;
        p.psdOverlap       = extraParams.overlap;
        p.windowType       = extraParams.window;
        p.aperiodicMode    = extraParams.aperiodicMode;
        p.fitFreqRange     = extraParams.fitFreqRange;
        p.maxPeaks         = extraParams.maxPeaks;
        p.minPeakHeight    = extraParams.minPeakHeight;
        p.peakThreshold    = extraParams.peakThreshold;
        p.peakWidthLimits  = extraParams.peakWidthLimits;
        p.correctAperiodic = extraParams.correctAperiodic;
        p.timeResolved     = extraParams.timeResolved;
        p.slidWinSec       = extraParams.slidWinSec;
        p.slidOverlap      = extraParams.slidOverlap;
        p.slidNAvg         = extraParams.slidNAvg;
        p.slidPeakWidthLimits = extraParams.slidPeakWidthLimits;
    end

else

    % --- Name-value mode ---
    if mod(numel(varargin), 2) ~= 0
        error('Options must be provided as name-value pairs.');
    end
    for k = 1:2:numel(varargin)
        key = lower(string(varargin{k}));
        val = varargin{k+1};
        switch key
            case 'measure',          p.measure          = val;
            case 'domain'
                dm = lower(char(val));
                switch dm
                    case {'channel','channels','chan','scalp'}, p.domain = 'channel';
                    case {'ica','ic','comp','components'},      p.domain = 'ica';
                    otherwise
                        error('domain must be ''channel'' or ''ica'' (got ''%s'').', dm);
                end
            case 'chanlist'
                if ischar(val) || isstring(val)
                    if strcmpi(string(val), 'all')
                        p.chanlist = {EEG.chanlocs.labels}';
                    else
                        p.chanlist = regexp(char(val), '[^,\s]+', 'match')';
                    end
                elseif iscell(val),  p.chanlist         = val(:);
                else, error('chanlist must be ''all'', a string, or a cellstr.');
                end
            case 'tau',              p.tau              = double(val);
            case 'm',                p.m                = double(val);
            case 'r',                p.r                = double(val);
            case 'vis',              p.vis              = logical(val);
            case 'parallel',         p.parallel         = logical(val);
            case 'progress',         p.progress         = logical(val);
            case 'coarsing',         p.coarsing         = val;
            case 'num_scales',       p.num_scales       = double(val);
            case 'filter_mode',      p.filter_mode      = val;
            case 'timewin',          p.TimeWin          = double(val);
            case 'timestep',         p.TimeStep         = double(val);
            case 'timeonly',         p.TimeOnly         = logical(val);
            case 'n',                p.n                = double(val);
            case 'kernel',           p.kernel_meth      = val;
            case 'fuzzy_mode',       p.fuzzy_mode       = val;
            case 'blocksize',        p.blocksize        = double(val);
            case 'freqrange',        p.freqRange        = double(val);
            case 'winsec',           p.winSec           = double(val);
            case 'overlap',          p.psdOverlap       = double(val);
            case 'window',           p.windowType       = val;
            case 'aperiodicmode',    p.aperiodicMode    = val;
            case 'fitfreqrange',     p.fitFreqRange     = double(val);
            case 'maxpeaks',         p.maxPeaks         = double(val);
            case 'minpeakheight',    p.minPeakHeight    = double(val);
            case 'peakthreshold',    p.peakThreshold    = double(val);
            case 'peakwidthlimits',  p.peakWidthLimits  = double(val);
            case 'correctaperiodic', p.correctAperiodic = logical(val);
            case 'alphaband',        p.alphaBand        = double(val);
            case 'timeresolved',     p.timeResolved     = logical(val);
            case 'slidwinsec',       p.slidWinSec       = double(val);
            case 'slidoverlap',      p.slidOverlap      = double(val);
            case 'slidnavg',         p.slidNAvg         = double(val);
            case 'slidPeakWidthLimits', p.slidPeakWidthLimits         = double(val);
                
            otherwise
                error('Unknown option: %s', key);
        end
    end
end

% ----------------------------
% Apply defaults
% ----------------------------

if isempty(p.domain)
    p.domain = 'channel';
end
if isempty(p.chanlist)
    if strcmpi(p.domain,'channel')
        disp('No channels selected: using all channels (default).');
    end
    p.chanlist = {EEG.chanlocs.labels}';
end
if isempty(p.measure)
    disp('No measure selected: using RCMFE (default).');
    p.measure = 'RCMFE';
end
if isempty(p.tau)
    disp('No time lag selected: using tau = 1 (default).');
    p.tau = 1;
end
if isempty(p.m)
    disp('No embedding dimension selected: using m = 2 (default).');
    p.m = 2;
end
if isempty(p.r)
    disp('No similarity bound selected: using r = 0.15 (default).');
    p.r = 0.15;
end
if isempty(p.vis),       p.vis       = true;  end
if isempty(p.progress)
    disp('Progress tracking not set: ON (default).');
    p.progress = true;
end
if isempty(p.parallel)
    disp('Parallel computing not set: ON (default).');
    p.parallel = true;
end

% Multiscale defaults
if contains(lower(p.measure), {'mse','mmse','mfe','cmfe','rcmfe','rcmvmfe'})
    if isempty(p.coarsing)
        disp('No coarse-graining method selected: using Standard Deviation (default).');
        p.coarsing = 'mean';
    end
    if isempty(p.num_scales)
        disp('Number of scales not set: using 30 (default).');
        p.num_scales = 30;
    end
    if isempty(p.filter_mode)
        if strcmpi(p.measure, 'mmse')
            p.filter_mode = 'narrowband';
            disp('Scale filtering: narrowband (Kosciessa et al. 2020).');
        else
            p.filter_mode = 'none';
        end
    end
    if isempty(p.TimeWin),  p.TimeWin  = []; end
    if isempty(p.TimeStep), p.TimeStep = []; end
    if isempty(p.TimeOnly), p.TimeOnly = false; end
end

% Fuzzy defaults
if contains(lower(p.measure), {'fuzzen','mfe','cmfe','rcmfe','rcmvmfe'})
    if isempty(p.n)
        disp('No fuzzy power selected: using n = 2 (default).');
        p.n = 2;
    end
    if isempty(p.kernel_meth)
        disp('Kernel method not set: using exponential (default).');
        p.kernel_meth = 'exponential';
    end
    if isempty(p.blocksize)
        disp('BlockSize not set: using 256 (default).');
        p.blocksize = 256;
    end
    if isempty(p.fuzzy_mode)
        disp('Fuzzy mode not set: using local (default).');
        p.fuzzy_mode = 'local';
    end
end

% Aperiodic defaults
if strcmpi(p.measure, 'aperiodic')
    if isempty(p.freqRange),        p.freqRange        = [1 40];       end
    if isempty(p.winSec),           p.winSec           = 4;            end
    if isempty(p.psdOverlap),       p.psdOverlap       = 0.5;          end
    if isempty(p.windowType),       p.windowType       = 'hann';       end
    if isempty(p.aperiodicMode),    p.aperiodicMode    = 'fixed';      end
    if isempty(p.fitFreqRange),     p.fitFreqRange     = p.freqRange;  end
    if isempty(p.maxPeaks),         p.maxPeaks         = 6;            end
    if isempty(p.minPeakHeight),    p.minPeakHeight    = 0.05;         end
    if isempty(p.peakThreshold),    p.peakThreshold    = 2.0;          end
    if isempty(p.peakWidthLimits),  p.peakWidthLimits  = [1 12];       end
    if isempty(p.correctAperiodic), p.correctAperiodic = true;         end
    if isempty(p.alphaBand),        p.alphaBand        = [8 13];       end
    if isempty(p.timeResolved),     p.timeResolved     = false;        end
    if isempty(p.slidWinSec),       p.slidWinSec       = 2;            end
    if isempty(p.slidOverlap),      p.slidOverlap      = 0.5;          end
    if isempty(p.slidNAvg),         p.slidNAvg         = 5;            end
    if isempty(p.slidPeakWidthLimits), p.slidPeakWidthLimits = [2 8]; end

end