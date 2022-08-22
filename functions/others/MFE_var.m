function Out_MFE=MFE_var(x,m,r,n,tau,Scale )
%
%  This function calculates the multiscale fuzzy entropy (MFE) whose coarse-graining uses variance (MFE_var)
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
% Out_MFE: a vector showing the MFE_var of x
%
% Ref:
% [1] H. Azami and J. Escudero, "Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis", Medical & Biological Engineering &
% Computing, 2016.
%
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

Out_MFE=NaN*ones(1,Scale);

for N_S=2:Scale
    xs=Multi_var(x,N_S);
    Out_MFE(N_S)=FuzEn(xs,m,r,n,tau);
end

end

function M_Data = Multi_var(Data,S)

%  the coarse-graining process based on variance
%  Input:   Data: time series;
%           S: the scale factor

% Output: 
%           M_Data: the new coarse-grained time series at the scale factor S
L = length(Data);

y=reshape(Data(1:floor(L/S)*S),S,[]);

M_Data=var(y);

end