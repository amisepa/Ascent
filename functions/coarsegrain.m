function cg = coarsegrain(Y, coarseType)
% coarsegrain  Coarse-grain a windowed matrix down dim 1.
%   Y is [s x nBins] (s samples per non-overlapping window). Returns a
%   1 x nBins row vector, one value per window, using the requested statistic.
%   Shared helper for the ASCENT multiscale entropy measures
%   (MSE, MFE, mMSE, CMFE, RCMFE).
switch lower(strtrim(coarseType))
    case 'mean'
        cg = mean(Y, 1, 'omitnan');
    case 'median'
        cg = median(Y, 1, 'omitnan');
    case {'sd','std','standard deviation'}
        cg = std(Y, 0, 1, 'omitnan');
    case {'var','variance'}
        cg = var(Y, 0, 1, 'omitnan');
    otherwise
        error('ascent:coarsegrain:unknownType', 'Unknown coarsing "%s".', coarseType);
end
cg = cg(:).';
end