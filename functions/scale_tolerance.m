function r_s = scale_tolerance(sig, s, coarseType, r, zNorm, minBinsHard, m)
% scale_tolerance  One tolerance per scale, shared across offsets, for the
%   composite multiscale fuzzy entropy measures (CMFE, RCMFE). The rescaling
%   factor is taken from the unshifted (offset-1) coarse-grained series so that
%   composite averaging (CMFE) and phi-pooling (RCMFE) stay coherent: every
%   offset at a given scale uses the same r. zNorm=0 returns the fixed r.
if zNorm == 0
    r_s = r;
    return
end
L = floor(numel(sig)/s) * s;
if L/s >= max(minBinsHard, m+1)
    cg_ref = coarsegrain(reshape(sig(1:L), s, []), coarseType);
    f = zNorm_scale(cg_ref, zNorm);
    if isfinite(f) && f > 0
        r_s = r * f;
    else
        r_s = r;
    end
else
    r_s = r;
end
end