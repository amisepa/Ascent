function Out_RCMSE = RCMSE_mu(x,m,r,tau,Scale)
%
%  This function calculates the refined composite multiscale sample entropy (RCMSE) whose coarse-graining uses mean (RCMSE_mu)
%
%
% Inputs:
%
% x: univariate signal - a vector of size 1 x N (the number of sample points)
% m: embedding dimension
% r: threshold (it is usually equal to 0.15 of the standard deviation of a signal - because we normalize signals to have a standard deviation of 1, here, r is usually equal to 0.15)
% tau: time lag (it is usually equal to 1)
% Scale: the number of scale factors
%
%
% Outputs:
%
% Out_RCMSE: a vector showing the RCMSE_mu of x
%
% Ref:
% [1] H. Azami and J. Escudero, "Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis", Medical & Biological Engineering &
% Computing, 2016.
% [2] S.D Wu, C.W. Wu, S.G. Lin, K.Y. Lee, and C.K. Peng, "Analysis of complex time series using refined composite multiscale entropy", Physics Letters A, vol.
% 378, no. 20, pp. 1369-1374, 2014
%
% If you use the code, please make sure that you cite references [1] and [2].
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  7-September-16
%%

% Signal is centered and normalised to standard deviation 1
x = x-mean(x);
x = x./std(x);

Out_RCMSE=NaN*ones(1,Scale);
Out_RCMSE(1)=SampEn(x,m,r,tau);

for i=2:Scale
    temp_A=[];
    temp_B=[];
    for ii=1:i
        
        xs = Multi_mu(x(ii:end),i);
        [SE,P] = SampEn(xs,m,r,tau);
        temp_A=[temp_A P(1)];
        temp_B=[temp_B P(2)];
        
    end
    A=sum(temp_A);
    B=sum(temp_B);
    Out_RCMSE(i)=log(A/B);
    
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
