function lbl = coarse_label(coarseType)
% coarse_label  Short uppercase label for a coarse-graining type (display only).
%   Shared helper for the ASCENT multiscale entropy measures.
cl = lower(strtrim(coarseType));
if any(strcmp(cl, {'sd','std','standard deviation'}))
    lbl = 'STD';
elseif any(strcmp(cl, {'var','variance'}))
    lbl = 'VAR';
elseif strcmp(cl, 'mean')
    lbl = 'MEAN';
elseif strcmp(cl, 'median')
    lbl = 'MEDIAN';
else
    lbl = upper(coarseType);
end
end