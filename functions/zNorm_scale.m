function f = zNorm_scale(x, zNorm)
% zNorm_scale  Per-scale tolerance rescaling factor for varying-tolerance
%   multiscale entropy. Shared helper for the ASCENT measures.
%     0  fixed (factor 1) | 1  population std | 2  population var
%     3  mad(mean) | 4  mad(median)
switch zNorm
    case 1, f = std(x, 1);      % population standard deviation
    case 2, f = var(x, 1);      % population variance
    case 3, f = mad(x, 0);      % mean absolute deviation
    case 4, f = mad(x, 1);      % median absolute deviation
    otherwise, f = 1;
end
end