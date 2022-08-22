%% Computes sample entropy of a univariate signal. 
% Sample entropy quantifies the likelihood that a sequence of m consecutive 
% data points that matches another sequence of the same length (match within
% a tolerance of r) will still match the other sequence when their length 
% is increased of one sample (sequences of length m + 1).
% 
% Inputs:
%   x: univariate signal - a vector of size 1 x n (the number of sample points)
%   m: embedding dimension
%   r: threshold (it is usually equal to 0.15 of the standard deviation of a signal - because we normalize signals to have a standard deviation of 1, here, r is usually equal to 0.15)
%   tau: time lag (it is usually equal to 1)
%
% Outputs:
%   sampEn: scalar quantity - the SampEn of x
%   p: a vector of length 2 : [the total number of template matches of length m, the total number of forward matches of length m+1]
%
% Please cite:
%   Azami & Escudero, "Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis", 
%       Medical & Biological Engineering & Computing, 2016.
%   Richman & Moorman, "Physiological time-series analysis using approximate entropy and sample entropy"
%       American Journal of Physiology-Heart and Circulatory Physiology, vol. 278, no. 6, pp.H2039-H2049, 2000.
%
% Cedric Cannard, August 2022

function [sampEn,p] = sampEn(signal,m,r,tau)

if nargin < 4, tau = 1; end
if tau > 1 
    signal = downsample(signal, tau); 
end

n = length(signal);
p = zeros(1,2);
sMat = zeros(m+1,n-m);

for i = 1:m+1
    sMat(i,:) = signal(i:n-m+i-1);
end

for k = m:m+1
    count = zeros(1,n-m);
    tempMat = sMat(1:k,:);
    
    for i = 1:n-k
        % calculate Chebyshev distance without counting self-matches
        dist = max(abs(tempMat(:,i+1:n-m) - repmat(tempMat(:,i),1,n-m-i)));
        % calculate the number of distance values that are less than the threshold r
        D = (dist < r);
        count(i) = sum(D)/(n-m);
    end
    
    p(k-m+1) = sum(count)/(n-m);
end
sampEn = log(p(1)/p(2));

