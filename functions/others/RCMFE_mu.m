function Out_RCMFE = RCMFE_mu(x,m,r,n,tau,Scale)
%
%  This function calculates the refined composite multiscale fuzzy entropy (RCMFE) whose coarse-graining uses mean (RCMFE_mu)
%
%
% Inputs:
%
% x: univariate signal - a vector of size 1 x N (the number of sample points)
% m: embedding dimension
% r: threshold (it is usually equal to 0.15 of the standard deviation of a signal - because we normalize signals to have a standard deviation of 1, here, r is usually equal to 0.15)
% n: fuzzy power (it is usually equal to 2)
% tau: time lag (it is usually equal to 1)
% Scale: the number of scale factors
%
%
% Outputs:
%
% Out_RCMFE: a vector showing the RCMFE_mu of x
%
% Ref:
% [1] H. Azami and J. Escudero, "Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis", Medical & Biological Engineering &
% Computing, 2016.
%
% If you use the code, please make sure that you cite reference [1].
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  7-September-16
%%

% Signal is centered and normalised to standard deviation 1
x = x-mean(x);
x = x./std(x);

Out_RCMFE=NaN*ones(1,Scale);
Out_RCMFE(1)=FuzEn(x,m,r,n,tau);

for i=2:Scale
    temp_A=[];
    temp_B=[];
    for ii=1:i
        
        xs = Multi_mu(x(ii:end),i);
        [FE,P] = FuzEn(xs,m,r,n,tau);
        temp_A=[temp_A P(1)];
        temp_B=[temp_B P(2)];
        
    end
    A=sum(temp_A);
    B=sum(temp_B);
    Out_RCMFE(i)=log(A/B);
    
end



function M_Data = Multi_mu(Data,S)

%  the coarse-graining process based on mean
%  Input:   Data: time series;
%           S: the scale factor

% Output:
%           M_Data: the coarse-grained time series at the scale factor S

L = length(Data);
J = fix(L/S);

for i=1:J
    M_Data(i) = mean(Data((i-1)*S+1:i*S));
end
