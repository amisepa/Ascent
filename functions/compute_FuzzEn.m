function FuzzEn = compute_FuzzEn(data, varargin)
% compute_FuzzEn  Computes Fuzzy Entropy across multichannel data.

p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ndims(x) == 2);
p.addParameter('m', 2,                 @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('n', 2,                 @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('tau', 1,               @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('r', 0.15,              @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('Kernel','exponential', @(s) any(strcmpi(s,{'exponential','gaussian'})));
p.addParameter('BlockSize', 2000,      @(x) isnumeric(x) && isscalar(x) && x >= 100);
p.addParameter('pdistMaxGB', 2.0,      @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Parallel', true,       @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,       @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
n_exp        = p.Results.n;
tau          = p.Results.tau;
r            = p.Results.r;
kernelType   = lower(p.Results.Kernel);
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
        fprintf('FuzzEn: %d channel(s) | m=%g, tau=%g, r=%g, n=%g | kernel=%s | parallel=on\n', ...
            nchan, m, tau, r, n_exp, kernelType);
    else
        fprintf('FuzzEn: %d channel(s) | m=%g, tau=%g, r=%g, n=%g | kernel=%s | parallel=off\n', ...
            nchan, m, tau, r, n_exp, kernelType);
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
        [fe, ~, ~] = fuzz_engine_raw(data(iChan,:), m, r, n_exp, tau, kernelType, true, blockSize, pdistMaxGB);
        FuzzEn(iChan) = fe;
        if useDQ
            send(dq, 1);
        end
    end
else
    for iChan = 1:nchan
        [fe, ~, ~] = fuzz_engine_raw(data(iChan,:), m, r, n_exp, tau, kernelType, true, blockSize, pdistMaxGB);
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