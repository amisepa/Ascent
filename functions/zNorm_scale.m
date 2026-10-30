function f = zNorm_scale(x, zNorm, n)
% zNorm_scale  Per-scale tolerance rescaling factor for varying-tolerance
%   multiscale entropy. Shared helper for the ASCENT measures.
%     0  fixed (factor 1) | 1  population std | 2  population var
%     3  mad(mean) | 4  mad(median)
%
%   n (optional): power of the spread in the tolerance. The exponential fuzzy
%   kernel exp(-d^n/r) has r in units of x^n, so the std must enter as std^n to
%   be equivalent to z-scoring each series. Pass 1 for the Gaussian kernel, whose
%   r is in units of x. Without n: the Chebyshev bound of sample entropy, in units
%   of x (zNorm 1 is then exact). zNorm 2 is always r*var: exact only at n = 2.
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).
if nargin < 3, n = []; end
switch zNorm
    case 1, f = std(x, 1);      % population standard deviation
    case 2, f = var(x, 1);      % population variance
    case 3, f = mad(x, 0);      % mean absolute deviation
    case 4, f = mad(x, 1);      % median absolute deviation
    otherwise, f = 1; return
end
if ~isempty(n) && zNorm ~= 2
    f = f^n;
end
end
