function [Out_SampEn,P] = SampEn(x,m,r,tau)
%
% This function calculates sample entropy (SampEn) of a univariate signal x
%
% Inputs:
%
% x: univariate signal - a vector of size 1 x N (the number of sample points)
% m: embedding dimension
% r: threshold (it is usually equal to 0.15 of the standard deviation of a signal - because we normalize signals to have a standard deviation of 1, here, r is usually equal to 0.15)
% tau: time lag (it is usually equal to 1)
%
% Outputs:
%
% Out_SampEn: scalar quantity - the SampEn of x
% P: a vector of length 2 : [the total number of template matches of length m, the total number of forward matches of length m+1]
%
%
% Ref:
% [1] H. Azami and J. Escudero, "Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis", Medical & Biological Engineering &
% Computing, 2016.
% [2] J. S. Richman and J. R. Moorman, "Physiological time-series analysis using approximate entropy and sample entropy"
% American Journal of Physiology-Heart and Circulatory Physiology, vol. 278, no. 6, pp.H2039-H2049, 2000.
%
%
% If you use the code, please make sure that you cite references [1] and [2].
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  7-September-16
%%
 
if nargin < 4, tau = 1; end
if tau > 1, x = downsample(x, tau); end

N = length(x);
P = zeros(1,2);
xMat = zeros(m+1,N-m);
for i = 1:m+1
    xMat(i,:) = x(i:N-m+i-1);
end

for k = m:m+1
    count = zeros(1,N-m);
    tempMat = xMat(1:k,:);
    
    for i = 1:N-k
        % calculate Chebyshev distance without counting self-matches
        dist = max(abs(tempMat(:,i+1:N-m) - repmat(tempMat(:,i),1,N-m-i)));
        % calculate the number of distance values that are less than the threshold r
        D = (dist < r);
        count(i) = sum(D)/(N-m);
    end
    
    P(k-m+1) = sum(count)/(N-m);
end
Out_SampEn = log(P(1)/P(2));
end