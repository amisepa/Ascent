function [S, info] = compute_surrogates(X, varargin)
% compute_surrogates  Spectrum-preserving surrogates of multichannel data.
%
%   [S, info] = compute_surrogates(X, 'Method', 'gaussphase', 'N', 19, 'Seed', 1)
%
% A surrogate keeps the linear properties of a signal (its power spectrum)
% and destroys any nonlinear structure. An entropy that differs between a
% signal and its surrogates therefore reflects more than the spectrum
% (Theiler et al. 1992; Kosciessa et al. 2020; Courtiol et al. 2016).
%
% Inputs:
%   X          : [channels x samples]
%   'Method'   : 'gaussphase' (default) each channel is first rank-transformed to
%                a Gaussian distribution, then phase-randomized as for 'phase'.
%                Null hypothesis: a linear Gaussian process seen through a fixed
%                monotonic distortion, i.e. the amplitude distribution is not
%                tested. The surrogates belong to the transformed signal
%                (info.tested): compare a measure of info.tested with them.
%                'phase' Fourier phase randomization (Theiler et al. 1992), with
%                the same random phases for all channels (Prichard & Theiler
%                1994): keeps each power spectrum and the cross-spectra exactly;
%                the amplitudes become Gaussian, so a non-Gaussian amplitude
%                distribution alone makes a measure differ. Null: a linear
%                Gaussian process.
%                'iaaft' iterative amplitude-adjusted Fourier transform
%                (Schreiber & Schmitz 1996): keeps each channel's amplitude
%                distribution exactly and its spectrum closely; channels are
%                randomized independently, so cross-channel correlations are
%                not kept.
%                Measured (SampEn, 100 signals of 3000 samples per case, 19
%                surrogates, ends matched; % of tests with p <= .05):
%                                         gaussphase   phase   iaaft
%                  linear Gaussian            7          9       24
%                  linear, then exp()         8        100      100
%                  noisy logistic map       100        100      100
%   'N'        : number of surrogates (default 1)
%   'Seed'     : seed of a private random stream (default: the global stream)
%   'Stream'   : a RandStream to draw from (overrides 'Seed'), e.g. to make
%                one surrogate at a time from one reproducible sequence
%   'MaxIter'  : IAAFT iterations at most (default 100); it stops earlier when
%                no rank changes
%
% Outputs:
%   S          : [channels x samples x N] surrogates
%   info       : .method, .tested (the signal the surrogates belong to: X, or
%                X rank-transformed to Gaussian for 'gaussphase'), .iterations
%                (IAAFT, per surrogate), .spectrum_error (IAAFT: relative RMS
%                difference of the amplitude spectra)
%
% The FFT treats the signal as periodic: a jump between the first and last
% samples adds broadband power to the surrogates. Trim the data first with
% surrogate_segment, as ascent_compute does.
%
% Missing samples are filled by linear interpolation; channels with no
% finite samples give NaN surrogates.
%
% References:
%   Theiler, J., Eubank, S., Longtin, A., Galdrikian, B., & Farmer, J. D. (1992).
%     Testing for nonlinearity in time series: The method of surrogate data.
%     Physica D, 58, 77-94.
%   Prichard, D., & Theiler, J. (1994). Generating surrogate data for time series
%     with several simultaneously measured variables. Phys Rev Lett, 73, 951-954.
%   Schreiber, T., & Schmitz, A. (1996). Improved surrogate data for nonlinearity
%     tests. Phys Rev Lett, 77, 635-638.
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).

p = inputParser;
p.addParameter('Method', 'gaussphase', @(s) any(strcmpi(s, {'gaussphase','phase','iaaft'})));
p.addParameter('N', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == round(x));
p.addParameter('Seed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('Stream', [], @(x) isempty(x) || isa(x, 'RandStream'));
p.addParameter('MaxIter', 100, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.parse(varargin{:});
o = p.Results;
method = lower(o.Method);

if isvector(X), X = X(:).'; end
X = double(X);
[nCh, nS] = size(X);
if nS < 4, error('compute_surrogates: at least 4 samples are needed.'); end

if ~isempty(o.Stream), rs = o.Stream;
elseif ~isempty(o.Seed), rs = RandStream('mt19937ar', 'Seed', o.Seed);
else, rs = RandStream.getGlobalStream;
end

% Missing samples
dead = all(~isfinite(X), 2);
for c = find(any(~isfinite(X), 2) & ~dead).'
    x = X(c,:); x(~isfinite(x)) = NaN;
    X(c,:) = fillmissing(x, 'linear', 'EndValues', 'nearest');
end
X(dead,:) = 0;
if strcmp(method, 'gaussphase')
    X = gaussianize(X);
    X(dead,:) = 0;
end

S    = zeros(nCh, nS, o.N);
info = struct('method', method, 'tested', X, 'iterations', nan(1, o.N), 'spectrum_error', nan(1, o.N));
info.tested(dead,:) = NaN;
F    = fft(X, [], 2);                 % once, reused by every surrogate
amp  = abs(F);

switch method
    case {'phase', 'gaussphase'}
        nh = floor((nS-1)/2);         % positive frequencies, without DC and Nyquist
        pos = 2:nh+1; neg = nS:-1:nS-nh+1;
        for k = 1:o.N
            ph = exp(1i * 2*pi * rand(rs, 1, nh));   % same phases for every channel
            Fs = F;
            Fs(:,pos) = F(:,pos) .* ph;
            Fs(:,neg) = conj(Fs(:,pos));
            S(:,:,k) = ifft(Fs, [], 2, 'symmetric');
        end

    case 'iaaft'
        xs   = sort(X, 2);            % amplitude distribution of each channel
        rows = repmat((1:nCh).', 1, nS);
        for k = 1:o.N
            [~, ord] = sort(rand(rs, nCh, nS), 2);      % random shuffle per channel
            s = X(sub2ind([nCh nS], rows, ord));
            prevRank = [];
            for it = 1:o.MaxIter
                Fs = fft(s, [], 2);
                s  = ifft(amp .* exp(1i * angle(Fs)), [], 2, 'symmetric');   % impose the spectrum
                [~, rk] = sort(s, 2);
                s(sub2ind([nCh nS], rows, rk)) = xs;                        % impose the values
                if isequal(rk, prevRank), break, end
                prevRank = rk;
            end
            S(:,:,k) = s;
            info.iterations(k) = it;
            info.spectrum_error(k) = sqrt(sum((abs(fft(s, [], 2)) - amp).^2, 'all') / sum(amp.^2, 'all'));
        end
end
S(dead,:,:) = NaN;
end

function G = gaussianize(X)
% Rank-based inverse normal transform of each row (ties broken by order)
[nCh, nS] = size(X);
[~, ord] = sort(X, 2);
R = zeros(nCh, nS);
R(sub2ind([nCh nS], repmat((1:nCh).', 1, nS), ord)) = repmat(1:nS, nCh, 1);
G = sqrt(2) * erfinv(2 * (R - .5) / nS - 1);
end
