function RCmvMFE = compute_RCmvMFE(X,sm,sr,n,stau,nScales, coarsing, showProg)
% To improve the stabily and reliabilty of multivaraite multiscale fuzzy entropy (mvMFE) based on mean (mvMFE_mu), especially for short signals, we proposed refined composite mvMFE_mu (RCmvMFE_mu). In RCmvMFE_mu, for scale
%factor tau, tau different multivaraite time series, corresponding to different starting points of the coarse-graining process are created and the RCmvMFE_mu value is defined based on
%the averages of the total number of m- and m+1- dimensional matched vector pairs of those shifted sequences.
%
% Inputs:
% X: multivariate signal - a matrix of size nvar (the number of channels) x nsamp (the number of sample points for each channel)
% sm: scalar embedding value
% sr: scalar threshold value (it is usually equal to 0.15)
% n: fuzzy power (it is usually equal to 2)
% stau: scalar time lag  value (it is usually equal to 1)
% nScales: the number of scale factors
%
% Output:
% RCmvMFE: a scalar quantity
%
% Ref:
% [1] H. Azami and J. Escudero, "Refined Composite Multivariate Generalized Multiscale Fuzzy Entropy:
% A Tool for Complexity Analysis of Multichannel Signals", Physica A, 2016.
%
% If you use the code, please make sure that you cite reference [1].
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  10-June-16

% Because multi-channel signals may have different amplitude ranges, the distances calculated on embedded vectors may be biased
%toward the largest amplitude ranges variates. For this reason, we scale all the data channels to the same amplitude range and
% we normalize each data channel to unit standard deviation so that the total variation becomes equal to the number of channels
% or variables [1].
X = zscore(X')';
r   = sr * sum(std(X,0,2));
M   = sm   * ones(1,size(X,1));     % size(X,1) = nvar
tau = stau * ones(1,size(X,1));     % size(X,1) = nvar

RCmvMFE = NaN(1,nScales);



for iScale = 2:nScales
    % fprintf('scale %g/%g...\n', iScale, nScales)
    if showProg
        fprintf('  progress: starting scale %d/%d\n', s, S);
        tScale = tic;
    end


    % Preallocate scalars per shift iii
    fi_m1_vec = zeros(1, iScale);
    fi_m2_vec = zeros(1, iScale);

    for iii = 1:iScale
        Xs = coarsegrain(X(:, iii:end), iScale, coarsing);

        % [~, fi_m1, fi_m2] = compute_mvFE(Xs, M, r, n, tau);
        [~, fi_m1, fi_m2] = compute_mvFE_fast(Xs, M, r, n, tau, ...
            'blockSize', 2000, 'useSingle', true, 'progress', false);

        fi_m1_vec(iii) = fi_m1;
        fi_m2_vec(iii) = fi_m2;
    end

    AA = sum(fi_m1_vec);
    BB = sum(fi_m2_vec);

    RCmvMFE(iScale) = log(AA / BB);

    if showProg
        fprintf('  progress: finished scale %d/%d (%.1fs)\n', s, S, toc(tScale));
    end

end
end



function Y = coarsegrain(Data, S, method)
% 1. Coarse graining of multichannel data by mean or variance
% Inputs:
%   Data   : nvar x nsamp
%   S      : scale factor (positive integer)
%   method : 'mean' or 'var'
% Output:
%   Y      : nvar x floor(nsamp/S)

[nvar, L] = size(Data);
J = floor(L / S);

if J == 0
    Y = zeros(nvar, 0);
    return
end

Data = Data(:, 1:J*S);
Data = reshape(Data, nvar, S, J);   % nvar x S x J

switch lower(method)
    case 'mean'
        Y = squeeze(mean(Data, 2, 'omitnan'));
    case 'var'
        % If you want omitnan variance and your MATLAB supports it, use:
        % Y = squeeze(var(Data, 0, 2, 'omitnan'));
        Y = squeeze(var(Data, 0, 2));
    otherwise
        error('method must be ''mean'' or ''var''.');
end
end
