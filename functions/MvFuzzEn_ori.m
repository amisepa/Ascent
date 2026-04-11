function [MFuzz, B0, Bt, B1] = MvFuzzEn(Data, varargin)
% MvFuzzEn  estimates the multivariate fuzzy entropy of a multivariate dataset.
%
%   [MFuzz, B0, Bt, B1] = MvFuzzEn(Data) 
% 
%   Returns the multivariate fuzzy entropy estimate (``MFuzz``) and the 
%   average vector distances (``m``: ``B0``; joint total 
%   ``m+1`` subspace: ``Bt``; all possible ``m+1`` subspaces: ``B1``),
%   from the M multivariate sequences in ``Data`` using the default parameters: 
%   embedding dimension = 2*ones(M,1), time delay = ones(M,1), 
%   fuzzy membership function = "default", fuzzy function parameters= [0.2, 2],
%   logarithm = natural, data normalization = false, 
% 
%   .. attention::
%      The entropy value returned as ``MFuzz`` is estimated using the "full" 
%      method [i.e.  -log(Bt/B0)] which compares delay vectors across all possible ``m+1`` 
%      expansions of the embedding space as applied in [1][3]. Contrary to
%      conventional definitions of fuzzy entropy, this method does not provide a
%      lower bound of 0!!
%      Thus, it is possible to obtain negative entropy values for multivariate 
%      fuzzy entropy, even for stochastic processes...
%
%      Alternatively, one can calculate ``MFuzz`` via the "naive" method, 
%      which ensures a lower bound of 0, by using the average vector distances
%      for an individual ``m+1`` subspace (B1) [e.g. -log(B1(1)/B0)],
%      or the average for all ``m+1`` subspaces [i.e. -log(mean(B1)/B0)].
%
%   .. note::
%      To maximize the number of points in the embedding process, this algorithm 
%      uses N-max(m*tau) delay vectors and **not** N-max(m)*max(tau) as employed 
%      in [1] and [3].
% 
%   [MFuzz, B0, Bt, B1] = MvFuzzEn(Data, name, value, ...)
% 
%   Returns the multivariate fuzzy entropy estimates (``MFuzz``) estimated 
%   from the M multivariate data sequences in ``Data`` using the specified
%   name/value pair arguments:
% 
%       :Data:  - Multivariate dataset, NxM matrix of N (>10) observations (rows) and M (cols) univariate data sequences 
%       :m:     - Embedding Dimension, a vector of M positive integers
%       :tau:   - Time Delay, a vector of M positive integers
%       :Fx:    - Fuzzy function name, one of the following strings, {``'sigmoid'``, ``'modsampen'``, ``'default'``, ``'gudermannian'``, ``'bell'``, ``'triangular'``, ``'trapezoidal1'``, ``'trapezoidal2'``, ``'z_shaped'``, ``'gaussian'``, ``'constgaussian'``}
%       :r:     - Fuzzy function parameters, a 1 element scalar or a 2 element vector of positive values. The ``r`` parameters for each  fuzzy function are defined as follows:  (default: [.2 2])
%                 :sigmoid:     
%                             * r(1) = divisor of the exponential argument
%                             * r(2) = value subtracted from argument (pre-division)
%                 :modsampen:    
%                             *  r(1) = divisor of the exponential argument
%                             *  r(2) = value subtracted from argument (pre-division)
%                 :default:      
%                             *  r(1) = divisor of the exponential argument
%                             *  r(2) = argument exponent (pre-division)
%                 :gudermannian: 
%                             *  r  = a scalar whose value is the numerator of
%                               argument to gudermannian function:
%                               GD(x) = atan(tanh(r/x)). GD(x) is 
%                               normalised to have a maximum value of 1.                       
%                 :triangular:  
%                             * r = a positive scalar whose value is the threshold 
%                               (corner point) of the triangular function.
%                 :trapezoidal1:  
%                             * r = a positive scalar whose value corresponds
%                               to the upper (2r) and lower (r) corner points of the trapezoid.
%                 :trapezoidal2:  
%                             * r(1) = a value corresponding to the upper corner point of the trapezoid.
%                             * r(2) = a value corresponding to the lower corner point of the trapezoid.
%                 :z_shaped:  
%                             * r = a scalar whose value corresponds to the
%                               upper (2r) and lower (r) corner points of the z-shape.
%                 :bell:  
%                             * r(1) = divisor of the distance value
%                             * r(2) = exponent of generalized bell-shaped function
%                 :gaussian:  
%                             * r = a scalar whose value scales the slope of the Gaussian curve.
%                 :constgaussian:  
%                             * r = a scalar whose value defines the lower 
%                               threshod and shape of the Gaussian curve.
%       :Norm:  - Normalisation of all M sequences to unit variance, a boolean
%       :Logx:  - Logarithm base, a positive scalar  [default: natural]
%
%   For further information on the name/value paire arguments, see
%   the `EntropyHub guide <https://github.com/MattWillFlood/EntropyHub/blob/main/EntropyHub%20Guide.pdf>`_
%   
%   See also:
%       MvSampEn, FuzzEn, XFuzzEn, FuzzEn2D, MSEn, MvPermEn.
%   
%   References:
%      [1] Ahmed, Mosabber U., et al. 
%           "A multivariate multiscale fuzzy entropy algorithm with application
%           to uterine EMG complexity analysis." 
%           Entropy 19.1 (2016): 2.
% 
%      [2] Azami, Alberto Fernández, Javier Escudero. 
%           "Refined multiscale fuzzy entropy based on standard deviation for 
%           biomedical signal analysis." 
%           Medical & biological engineering & computing 55 (2017): 2037-2052.
% 
%      [3] Ahmed Mosabber Uddin, Danilo P. Mandic
%           "Multivariate multiscale entropy analysis."
%           IEEE signal processing letters 19.2 (2011): 91-94.
% 
narginchk(1,13)
p = inputParser;
Data = squeeze(Data);
Dn = size(Data,2);
N = size(Data,1);
Chk = @(x) isnumeric(x) && isvector(x) && (length(x)==Dn) && (min(x)>0) && all(mod(x,1)==0);
Chk2 = @(x) isnumeric(x) && isscalar(x) && (x > 0);
Chk3 = ["sigmoid","modsampen","default","gudermannian", "bell", "triangular", "trapezoidal1", "trapezoidal2", "z_shaped", "gaussian", "constgaussian"];
addRequired(p,'Data',@(x) isnumeric(x) && ismatrix(x) && (Dn>1) && N>10);
addParameter(p,'m',2*ones(Dn,1),Chk);
addParameter(p,'tau',ones(Dn,1),Chk);
addParameter(p,'r',[.2, 2.0],@(x) isnumeric(x) && isvector(x) && (min(x)>=0) && length(x)<=2);
addParameter(p,'Fx','default',@(x) (ischar(x) || isstring(x)) && any(strcmpi(string(x),Chk3)));
addParameter(p,'Logx',exp(1),Chk2);
addParameter(p,'Norm',false,@(x) islogical(x))
parse(p,Data,varargin{:})
m = p.Results.m(:); tau = p.Results.tau(:); 
Fun = str2func(p.Results.Fx);
r = p.Results.r(:); Logx = p.Results.Logx; 
Norm = p.Results.Norm;
if Norm, Data = Data./std(Data,1); end
Nx = N - max((m-1).*tau);
Ny = N - max(m.*tau);
Vex = zeros(Nx,sum(m));
q = 1;
for k = 1:Dn
    for p=1:m(k)
        Vex(:,q) = Data(1+(p-1)*tau(k):Nx+(p-1)*tau(k),  k);
        q = q+ 1;
    end    
end
Count0 = Distx(Vex-mean(Vex, 2), r, Fun); 
B0 = sum(Count0(:))/(Nx*(Nx-1)/2);  
B1 = zeros(1,Dn);
Temp = cumsum(m);
Vez = [];
for k = 1:Dn
    Sig = Data(1+m(k)*tau(k):Ny+m(k)*tau(k), k);
    Vey = [Vex(1:Ny, 1:Temp(k)) Sig  Vex(1:Ny, Temp(k)+1:end)];
    Vez = [Vez; Vey];
    Count1 = Distx(Vey - mean(Vey,2), r, Fun);
    B1(k) = sum(Count1(:))/(Ny*(Ny-1)/2);
end
Count1 = Distx(Vez - mean(Vez,2), r, Fun);
Bt = sum(Count1(:))/(Dn*Ny*(Dn*Ny-1)/2);
MFuzz = -log(Bt/B0)/log(Logx);
end
function [Counter] = Distx(Vex, r, Fun)
    Nt = size(Vex,1);
    Counter = zeros(Nt-1);
    for x=1:Nt-1
        Counter(x,x:end) = Fun(max(abs(Vex(x+1:end,:) - Vex(x,:)),[],2), r);
    end
end
function [y] = sigmoid(x,r)
if numel(r) == 1
    error('When Fx = "Sigmoid", r must be a two-element vector.')
end
y = 1./(1+exp((x-r(2))/r(1)));
end
function [y] = modsampen(x,r)
if numel(r) == 1
    error('When Fx = "Modsampen", r must be a two-element vector.')
end
y = 1./(1+exp((x-r(2))/r(1)));
end
function [y] = default(x,r)
if numel(r) == 1
    error('When Fx = "Default", r must be a two-element vector.')
end
y = exp(-(x.^r(2))/r(1));
end
function [y] = gudermannian(x,r)
if r <= 0
    error('When Fx = "Gudermannian", r must be a scalar > 0.')
end
y = atan(tanh(r(1)./x));
y = y/max(y);
end
function [y] = triangular(x,r)
    if numel(r) > 1
        error('When Fx = "Triangular", r must be a scalar > 0.')
    end
    y = 1 - (x/r);
    y(x > r) = 0;
end
function [y] = trapezoidal1(x,r)
    if numel(r) > 1
        error('When Fx = "Trapezoidal1", r must be a scalar > 0.')
    end
    y = zeros(size(x));
    y(x <= r*2) = 2 - (x(x <= r*2)/r);
    y(x <= r) = 1;
end
function [y] = trapezoidal2(x,r)
    if numel(r) ~=2 
        error('When Fx = "Trapezoidal2", r must be a two-element vector.')
    end
    y = zeros(size(x));
    y(x <= max(r)) = 1 - (x(x <= max(r))-min(r))/range(r);
    y(x <= min(r)) = 1;
end
function [y] = z_shaped(x, r)
    if numel(r) > 1
        error('When Fx = "Z_shaped", r must be a scalar > 0.')
    end
    y = zeros(size(x));
    y(x <= 2*r) = 2*(((x(x <= 2*r) - 2*r)/r).^2);
    y(x <= 1.5*r) = 1 - (2*(((x(x <= 1.5*r) - r)/r).^2));
    y(x <= r) = 1;
end
function [y] = bell(x, r)
    if numel(r) ~=2 
        error('When Fx = "Bell", r must be a two-element vector.')
    end
    y = 1./(1 + abs(x/r(1)).^(2*r(2)));
end
function [y] = gaussian(x, r)
    if numel(r) > 1
        error('When Fx = "Gaussian", r must be a scalar > 0.')
    end
    y = exp(-((x.^2)/(2*(r^2))));
end
function [y] = constgaussian(x, r)
    if numel(r) > 1
        error('When Fx = "ConstGaussian", r must be a scalar > 0.')
    end
    y = ones(size(x));
    y(x > r) = exp(-log(2)*((x(x > r) - r)/r).^2);
end
%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub