function FuzzEn = compute_FuzzEn(data, varargin)
% compute_FuzzEn  Compute Fuzzy Entropy across channels.
%
%   FuzzEn = compute_FuzzEn(data, 'm', 2, 'n', 2, 'tau', 1, 'r', 0.15, ...
%                           'Mode', 'local', 'Kernel', 'exponential')
%
% Inputs
%   data      : numeric matrix [n_channels x n_samples] or [n_samples x n_channels]
%
% Name-value parameters
%   'm'         : embedding dimension (default = 2)
%   'n'         : fuzzy exponent for exponential kernel (default = 2)
%   'tau'       : time lag / embedding delay (default = 1)
%   'r'         : similarity bound (default = 0.15)
%   'Mode'      : 'local' | 'global'
%                 - 'local'  : subtract mean of each embedded vector before distance
%                              (original FuzEn-style local detrending)
%                 - 'global' : compare embedded vectors directly
%                              (SampEn-like / FuzEn(Glb)-style)
%                 default = 'local'
%   'Kernel'    : 'exponential' | 'gaussian' (default = 'exponential')
%   'Normalize' : true/false, z-score each channel before FuzzEn (default = true)
%   'BlockSize' : block size for exact blocked fallback (default = 2000)
%   'pdistMaxGB': max temporary GB allowed for pdist path (default = 2.0)
%   'Parallel'  : true/false, parfor over channels when available (default = true)
%   'Progress'  : true/false, print progress (default = true)
%
% Output
%   FuzzEn    : [n_channels x 1] Fuzzy Entropy values
%
% Notes
%   • 'local' matches the original FuzEn-style practice of removing the mean
%     of each embedded vector before computing Chebyshev distances.
%   • 'global' keeps absolute/global offsets and is closer in spirit to SampEn.
%   • When Normalize=true and r=0.15, this corresponds to the common practice
%     of using r = 0.15 * SD(original signal).
%
% See also: fuzz_engine_raw, compute_CMFE, compute_RCMFE

p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ndims(x) == 2);
p.addParameter('m', 2,                    @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('n', 2,                    @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('tau', 1,                  @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('r', 0.15,                 @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 2);
p.addParameter('Mode', 'local',           @(s) any(strcmpi(s, {'local','global'})));
p.addParameter('Kernel', 'exponential',   @(s) any(strcmpi(s, {'exponential','gaussian'})));
p.addParameter('Normalize', true,         @(x) islogical(x) && isscalar(x));
p.addParameter('BlockSize', 2000,         @(x) isnumeric(x) && isscalar(x) && x >= 100);
p.addParameter('pdistMaxGB', 2.0,         @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Parallel', true,          @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,          @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
n_exp        = p.Results.n;
tau          = p.Results.tau;
r            = p.Results.r;
modeType     = lower(p.Results.Mode);
kernelType   = lower(p.Results.Kernel);
doNormalize  = p.Results.Normalize;
blockSize    = p.Results.BlockSize;
pdistMaxGB   = p.Results.pdistMaxGB;
parallelMode = p.Results.Parallel;
showProgress = p.Results.Progress;

if size(data,1) > size(data,2)
    data = data.';
end
[nchan, ~] = size(data);

if showProgress
    if parallelMode && ~isempty(ver('parallel'))
        fprintf('FuzzEn: %d ch | m=%g, tau=%g, r=%g, n=%g | mode=%s | kernel=%s | parallel=on\n', ...
            nchan, m, tau, r, n_exp, upper(modeType), kernelType);
    else
        fprintf('FuzzEn: %d ch | m=%g, tau=%g, r=%g, n=%g | mode=%s | kernel=%s | parallel=off\n', ...
            nchan, m, tau, r, n_exp, upper(modeType), kernelType);
    end
end

FuzzEn = nan(nchan, 1);

useWB = ~parallelMode && usejava('desktop') && showProgress;
hWB = [];
if useWB
    try
        hWB = waitbar(0, 'Computing Fuzzy Entropy...', 'Name', 'compute_FuzzEn');
    catch
        hWB = [];
    end
end

useDQ = parallelMode && ~isempty(ver('parallel')) && showProgress;
if useDQ
    dq = parallel.pool.DataQueue;
    nDone = 0;
    afterEach(dq, @notifyProgress);
end

if parallelMode && ~isempty(ver('parallel'))
    parfor iChan = 1:nchan
        [fe, ~, ~] = fuzz_engine_raw(data(iChan,:), m, r, n_exp, tau, ...
            kernelType, doNormalize, blockSize, pdistMaxGB, modeType);
        FuzzEn(iChan) = fe;
        if useDQ, send(dq, 1); end
    end
else
    for iChan = 1:nchan
        [fe, ~, ~] = fuzz_engine_raw(data(iChan,:), m, r, n_exp, tau, ...
            kernelType, doNormalize, blockSize, pdistMaxGB, modeType);
        FuzzEn(iChan) = fe;

        if ~isempty(hWB) && isvalid(hWB)
            try
                waitbar(iChan / nchan, hWB, sprintf('Computing Fuzzy Entropy... (%d/%d)', iChan, nchan));
            catch
            end
        end
    end

    if ~isempty(hWB) && isvalid(hWB)
        try
            close(hWB);
        catch
        end
    end
end

    function notifyProgress(~)
        nDone = nDone + 1;
        step = max(1, round(0.05 * nchan));
        if nDone == 1 || nDone == nchan || mod(nDone, step) == 0
            fprintf('  progress: ch %d/%d\n', nDone, nchan);
        end
    end
end