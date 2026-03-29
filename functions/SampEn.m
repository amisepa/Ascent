function [Samp, A, B, Vcp_Ka_Kb] = SampEn(Sig, varargin)
% SampEn  estimates the sample entropy of a univariate data sequence.
%
%   [Samp, A, B] = SampEn(Sig) 
% 
%   Returns the sample entropy estimates (``Samp``) and the number of matched state 
%   vectors (``m: B``, ``m+1: A``) for ``m`` = [0,1,2] estimated from the data 
%   sequence (``Sig``) using the default parameters: embedding dimension = 2, 
%   time delay = 1,  radius threshold = 0.2*SD(``Sig``), logarithm = natural
% 
%   [Samp, A, B, (Vcp, Ka, Kb)] = SampEn(Sig, ..., Vcp = true) 
% 
%   If ``Vcp == True``, an additional vector ``(Vcp, Ka, Kb)`` is returned with    
%   the sample entropy estimates (``Samp``) and the number of matched state
%   vectors (``m: B``, ``m+1: A``). ``(Vcp, Ka, Kb)``  contains the variance of the conditional
%   probabilities (``Vcp``), i.e. CP = A/B,  and the number of **overlapping**
%   matching vector pairs of lengths m+1 (``Ka``) and m (``Kb``),
%   respectively.  Note ``Vcp`` is undefined for the zeroth embedding dimension (m = 0) 
%   and due to the computational demand, **will take substantially more time to return function outputs.**
%   See Appendix B in [2] for more info.
%
%   [Samp, A, B] = SampEn(Sig, name, value, ...)
% 
%   Returns the sample entropy estimates (``Samp``) for dimensions = [0,1,..., ``m``]
%   estimated from the data sequence (``Sig``) using the specified name/value pair
%   arguments:
% 
%       * ``m``     - Embedding Dimension, a positive integer
%       * ``tau``   - Time Delay, a positive integer
%       * ``r``     - Radius Distance Threshold, a positive scalar  
%       * ``Logx``  - Logarithm base, a positive scalar  
%       * ``Vcp``   - Option to return variance of conditional probabilities and the number of overlapping matching vector pairs, a boolean [default: false]
%
%   See also:
%       ApEn, FuzzEn, PermEn, CondEn, XSampEn, SampEn2D, MSEn.
%   
%   References:
%      [1] Joshua S Richman and J. Randall Moorman. 
%           "Physiological time-series analysis using approximate entropy
%           and sample entropy." 
%           American Journal of Physiology-Heart and Circulatory Physiology (2000).
% 
%      [2]  Douglas E Lake, Joshua S Richman, M.P. Griffin, J. Randall Moorman
%           "Sample entropy analysis of neonatal heart rate variability."
%           American Journal of Physiology-Regulatory, Integrative and Comparative Physiology
%           283, no. 3 (2002): R789-R797.
% 
narginchk(1,11)
if size(Sig,1) == 1
    Sig = Sig';
end
Sig = squeeze(Sig);
p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isscalar(x) && (x > 0);
addRequired(p,'Sig',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addParameter(p,'m',2,Chk);
addParameter(p,'tau',1,Chk);
addParameter(p,'r',.2*std(Sig,1),Chk2);
addParameter(p,'Logx',exp(1),Chk2);
addParameter(p,'Vcp',false,@(x) islogical(x));
parse(p,Sig,varargin{:})
m = p.Results.m; tau = p.Results.tau; 
r = p.Results.r; Logx = p.Results.Logx; 
Vcp = p.Results.Vcp;
N = length(Sig);
Counter = (abs(Sig - Sig') <= r).*(triu(ones(N),1));
M = [m*ones(1,N-(m*tau)) repelem((m-1):-1:1,tau)];
A(1) = sum(sum(Counter));  B(1) = N*(N-1)/2;
for n = 1:N - tau
    ix = find(Counter(n,:)==1);
    for k = 1:M(n)
        ix(ix + (k*tau) > N) = [];  
        if isempty(ix)
            break
        end
        p1 = repmat(Sig(n:tau:n+(k*tau))',length(ix),1);
        p2 = Sig(ix+(0:tau:tau*k)')';        
        ix = ix(max(abs(p1 - p2),[],2) <= r);
        Counter(n, ix) = Counter(n, ix)+1;
    end
end
for k = 1:m
    A(k+1) = sum(sum(Counter>k));
    B(k+1) = sum(sum(Counter(:,1:N-(k*tau))>=k));
end
Samp = -log(A./B)/log(Logx);
Vcp_Ka_Kb = [];
if Vcp
    [Temp1,Temp2] = find(Counter>m);
    if numel(Temp1)>1
        Temp1 = cast(Temp1,'int32'); Temp2 = cast(Temp2,'int32');
        Temp = [Temp1 Temp2];
        Ka = zeros(numel(Temp1)-1,1);
        for k = 1:length(Ka)
            TF = (abs(Temp(k+1:end,:) - Temp1(k)) <= m*tau) + (abs(Temp(k+1:end,:) - Temp2(k)) <= m*tau);
            Ka(k) = sum(any(TF,2));
        end
    else
        Ka = 0;
    end
    
    [Temp1,Temp2] = find(Counter(:,1:end-m*tau)>=m);
    if numel(Temp1)>1
        Temp1 = cast(Temp1,'int32'); Temp2 = cast(Temp2,'int32');
        Temp = [Temp1 Temp2];
        Kb = zeros(numel(Temp1)-1,1);
        for k = 1:length(Kb)
            TF = (abs(Temp(k+1:end,:) - Temp1(k)) <= (m-1)*tau) + (abs(Temp(k+1:end,:) - Temp2(k)) <= (m-1)*tau);
            Kb(k) = sum(any(TF,2));      
        end
    else
        Kb = 0;
    end
        
    Ka = sum(Ka(:));
    Kb = sum(Kb(:));
    CP = A(end)/B(end);
    Vcp = (CP*(1-CP)/B(end)) + (Ka - Kb*(CP^2))/(B(end)^2); 
    
    Vcp_Ka_Kb = [Vcp, Ka, Kb];
end
end
%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub