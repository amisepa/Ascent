function [MSx, CI] = MvMSEn(Data, Mobj, varargin) 
% MvMSEn  returns the multivariate multiscale entropy of a multivariate dataset.
%
%   [MSx, CI] = MvMSEn(Data, Mobj) 
% 
%   Returns a vector of multivariate multiscale entropy values (``MSx``) and the complexity 
%   index (``CI``) of the data sequences ``Data`` using the parameters specified 
%   by the multiscale object (``Mobj``) over 3 temporal scales with coarse-
%   graining (default). 
% 
%   .. caution::
%          By default, the ``MvSampEn`` and ``MvFuzzEn`` multivariate entropy algorithms
%          estimate entropy values using the "full"  method by comparing delay vectors 
%          across all possible ``m+1`` expansions of the embedding space as applied in [1].
%          These methods are not lower-bounded to 0, like most entropy algorithms,
%          so ``MvMSEn`` may return negative entropy values if the base multivariate 
%          entropy function is ``MvSampEn`` and ``MvFuzzEn``, even for stochastic processes...
% 
% 
%   [MSx,CI] = MvMSEn(Data, Mobj, name, value, ...)
% 
%   Returns a vector of multivariate multiscale entropy values (``MSx``) and the complexity 
%   index (``CI``) of the data sequences ``Data`` using the parameters specified by
%   the multiscale object (``Mobj``) and the following name/value pair arguments:
%   
%       :Scales:   - Number of temporal scales, an integer > 1   (default = 3)
%       :Methodx:  - Graining method, one of the following, [default = ``'coarse'``]  {``'coarse'``,``'generalized'``, ``'modified'``} 
%       :Plotx:    - When ``Plotx == true``, returns a plot of the entropy value at each time scale (i.e. the multiscale entropy curve) [default: false]
% 
%   For further information on multiscale graining procedures, see
%   the `EntropyHub guide <https://github.com/MattWillFlood/EntropyHub/blob/main/EntropyHub%20Guide.pdf>`_
% 
%   See also:
%       MSobject, MSEn, cMvMSEn, MvSampEn, MvFuzzEn, MvPermEn, MvDispEn, MvCoSiEn
% 
%   References:
%      [1] Ahmed Mosabber Uddin, Danilo P. Mandic
%           "Multivariate multiscale entropy analysis."
%           IEEE signal processing letters 19.2 (2011): 91-94.
% 
%      [2] Madalena Costa, Ary Goldberger, and C-K. Peng,
%           "Multiscale entropy analysis of complex physiologic time series."
%           Physical review letters
%           89.6 (2002): 068102.
% 
%      [3] Vadim V. Nikulin, and Tom Brismar,
%           "Comment on "Multiscale entropy analysis of complex physiologic
%           time series"." 
%           Physical Review Letters 
%           92.8 (2004): 089803.
% 
%      [4] Madalena Costa, Ary L. Goldberger, and C-K. Peng. 
%           "Costa, Goldberger, and Peng reply." 
%           Physical Review Letters
%           92.8 (2004): 089804.
% 
%      [5] Madalena Costa, Ary L. Goldberger and C-K. Peng,
%           "Multiscale entropy analysis of biological signals." 
%           Physical review E 
%           71.2 (2005): 021906.
% 
%      [6] Ranjit A. Thuraisingham and Georg A. Gottwald,
%           "On multiscale entropy analysis for physiological data."
%           Physica A: Statistical Mechanics and its Applications
%           366 (2006): 323-332.
% 
%      [7] Ahmed Mosabber Uddin, Danilo P. Mandic
%           "Multivariate multiscale entropy: A tool for complexity
%           analysis of multichannel data."
%           Physical Review E 84.6 (2011): 061918.
%
narginchk(2,8)
Data = squeeze(Data);
Dn = size(Data,2);
p = inputParser;
Chk = ["coarse";"generalized";"modified";"timeshift"];
addRequired(p,'Data',@(x) isnumeric(x) && ismatrix(x) && (size(x,1)>10) && Dn>1);
addRequired(p,'Mobj',@(x) isstruct(x));
addParameter(p,'Methodx','coarse',@(x) any(strcmpi(string(x),Chk)));
addParameter(p,'Scales', 3, @(x) isscalar(x) && isnumeric(x) && (mod(x,1)==0) && (x>1));
addParameter(p,'Plotx', false, @(x) islogical(x));
parse(p, Data, Mobj, varargin{:})
Methodx = str2func(lower(p.Results.Methodx));
MSx = zeros(1,p.Results.Scales);
Fields = fieldnames(Mobj);
Y = struct2cell(Mobj);
C = [Fields(2:end),Y(2:end)].';
assert(startsWith(func2str(Mobj.Func), "Mv"),...
    "Multivariate entropy estimator is a not a multivariate method!")
for T = 1:p.Results.Scales       
    fprintf(' .')
    Temp = Methodx(Data, T, Dn);  
    Temp2 = Mobj.Func(Temp,C{:});    
    MSx(T) = mean(Temp2);  
end
if numel(Temp2)>1
    warning("More than one multivariate entropy value was found at each scale. "+ ...
        "The multivariate multiscale entropy returned at each scale represents the average of these values.")
end
CI = sum(MSx);
fprintf(' .\n')
if any(isnan(MSx))
    fprintf('\nSome entropy values may be undefined.')
end
if p.Results.Plotx
   figure, hold on
   plot(1:p.Results.Scales,MSx,'color',[8 63 77]/255,'LineWidth',3)   
   scatter(1:p.Results.Scales,MSx,60,[1 0 1],'filled')   
   xlabel('Scale Factor','FontSize',12,'FontWeight','bold','Color',[7 54 66]/255)
   ylabel('Entropy Value','FontSize',12,'FontWeight','bold','Color',[7 54 66]/255)
   Temp = func2str(Y{1});
   title(sprintf('Multivariate Multiscale %s (%s-graining method)',Temp(3:end),...
       func2str(Methodx)),'FontSize',16,'FontWeight','bold','Color',[7 54 66]/255)
   ylim([0 max(MSx)+.2])
end
end
function [Y] = generalized(Z, sx, Dn)
    Ns = floor(size(Z,1)/sx);
    Y = zeros(Ns, Dn);
    for k = 1:Dn
        Y(:,k) = var(reshape(Z(1:sx*Ns,k),sx,Ns),1,1);
    end
end
function [Y] = coarse(Z,sx, Dn)
    Ns = floor(length(Z)/sx);
    Y = zeros(Ns, Dn);
    for k = 1:Dn
        Y(:,k) = mean(reshape(Z(1:sx*Ns,k),sx,Ns),1);
    end
end
function [Y] = modified(Z,sx, Dn)
    Y = zeros(size(Z,1)-sx+1, Dn);
    for k = 1:Dn
        Y(:,k) = conv(Z(:,k),ones(1,sx),'valid')/sx;
    end
end
%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub