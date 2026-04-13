function cg = coarsegrain(Y, ct)
% coarsegrain  Column-wise coarse-graining.
%
%   Y  : [scale x nBins]
%   cg : [1 x nBins]
%
% Original-compatible modes:
%   'mean' : plain arithmetic mean over rows
%   'var'  : plain sample variance over rows (normalization N-1)
%
% Additional toolbox modes:
%   'median'
%   'trimmed mean' / 'trimmed' / 'tmean' / 'trim20'
%   'std' / 'sd' / 'standard deviation'

ct = lower(strtrim(ct));

switch ct
    case 'mean'
        % Match original coarse-graining behavior
        cg = mean(Y, 1);

    case {'variance','var'}
        % Match original variance coarse-graining behavior
        cg = var(Y, 0, 1);

    case 'median'
        % Robust optional mode
        cg = median(Y, 1, 'omitnan');

    case {'trimmed mean','trimmed','tmean','trim20'}
        % Robust optional mode: 20% total trimming (10% each tail)
        pct = 20;
        if exist('trimmean', 'file') == 2 && all(isfinite(Y(:)))
            cg = trimmean(Y, pct, 1);
        else
            nSeg = size(Y, 2);
            cg = nan(1, nSeg);
            kfrac = pct / 200;

            for j = 1:nSeg
                col = Y(:, j);
                col = col(isfinite(col));

                if isempty(col)
                    cg(j) = NaN;
                    continue
                end

                k = floor(kfrac * numel(col));
                if 2*k >= numel(col)
                    cg(j) = NaN;
                else
                    col = sort(col);
                    cg(j) = mean(col(k+1:end-k));
                end
            end
        end

    case {'sd','std','standard deviation'}
        % Optional spread-based coarse-graining
        cg = std(Y, 0, 1, 'omitnan');

    otherwise
        error('coarsegrain:BadMethod', ...
            'Unknown coarse-graining method "%s".', ct);
end

cg = cg(:).';
end