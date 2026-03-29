function [MSx, CI] = compute_MSE_ori(Sig, Mobj, varargin) 
% MSEn  returns the multiscale entropy of a univariate data sequence.
%
%   [MSx,CI] = MSEn(Sig, Mobj) 
% 
%   Returns a vector of multiscale entropy values (``MSx``) and the complexity 
%   index (``CI``) of the data sequence ``Sig`` using the parameters specified 
%   by the multiscale object (``Mobj``) over 3 temporal scales with coarse-
%   graining (default). 
% 
%   [MSx,CI] = MSEn(Sig, Mobj, name, value, ...)
% 
%   Returns a vector of multiscale entropy values (``MSx``) and the complexity 
%   index (``CI``) of the data sequence ``Sig`` using the parameters specified by
%   the multiscale object (``Mobj``) and the following name/value pair arguments:
%   
%       * ``Scales``   - Number of temporal scales, an integer > 1   (default = 3)
%       * ``Methodx``  - Graining method, one of the following: [default = ``'coarse'``]
%         {``'coarse'``,``'generalized'``, ``'modified'``, ``'imf'``, ``'timeshift'``} 
%       * ``RadNew``   - Radius rescaling method, an integer in the range [1 4].
%         When the entropy specified by ``Mobj`` is ``SampEn`` or ``ApEn``, 
%         RadNew rescales the radius threshold in each sub-sequence
%         at each time scale (Xt). If a radius value (``r``) is specified 
%         by ``Mobj``, this becomes the rescaling coefficient, otherwise
%         it is set to 0.2 (default). The value of ``RadNew`` specifies
%         one of the following methods:
%                   * [1]    Standard Deviation          - r*std(Xt)
%                   * [2]    Variance                    - r*var(Xt)
%                   * [3]    Mean Absolute Deviation     - r*mad(Xt)
%                   * [4]    Median Absolute Deviation   - r*mad(Xt,1)
%       * ``Plotx``    - When ``Plotx == true``, returns a plot of the entropy value at
%         each time scale (i.e. the multiscale entropy curve) [default: false]
% 
%   For further information on multiscale graining procedures, see
%   the `EntropyHub guide <https://github.com/MattWillFlood/EntropyHub/blob/main/EntropyHub%20Guide.pdf>`_
% 
%   See also:
%       MSobject, rMSEn, cMSEn, hMSEn, XMSEn, cXMSEn, hXMSEn, SampEn
% 
%   References:
%      [1] Madalena Costa, Ary Goldberger, and C-K. Peng,
%           "Multiscale entropy analysis of complex physiologic time series."
%           Physical review letters
%           89.6 (2002): 068102.
% 
%      [2] Vadim V. Nikulin, and Tom Brismar,
%           "Comment on “Multiscale entropy analysis of complex physiologic
%           time series”." 
%           Physical Review Letters 
%           92.8 (2004): 089803.
% 
%      [3] Madalena Costa, Ary L. Goldberger, and C-K. Peng. 
%           "Costa, Goldberger, and Peng reply." 
%           Physical Review Letters
%           92.8 (2004): 089804.
% 
%      [4] Madalena Costa, Ary L. Goldberger and C-K. Peng,
%           "Multiscale entropy analysis of biological signals." 
%           Physical review E 
%           71.2 (2005): 021906.
% 
%      [5] Ranjit A. Thuraisingham and Georg A. Gottwald,
%           "On multiscale entropy analysis for physiological data."
%           Physica A: Statistical Mechanics and its Applications
%           366 (2006): 323-332.
% 
%      [6] Meng Hu and Hualou Liang,
%           "Intrinsic mode entropy based on multivariate empirical mode
%           decomposition and its application to neural data analysis." 
%           Cognitive neurodynamics
%           5.3 (2011): 277-284.
% 
%      [7] Anne Humeau-Heurtier 
%           "The multiscale entropy algorithm and its variants: A review." 
%           Entropy 
%           17.5 (2015): 3110-3123.
% 
%      [8] Jianbo Gao, et al.,
%           "Multiscale entropy analysis of biological signals: a 
%           fundamental bi-scaling law." 
%           Frontiers in computational neuroscience 
%           9 (2015): 64.
% 
%      [9] Paolo Castiglioni, et al.,
%           "Multiscale Sample Entropy of cardiovascular signals: Does the
%           choice between fixed-or varying-tolerance among scales 
%           influence its evaluation and interpretation?." 
%           Entropy
%           19.11 (2017): 590.
% 
%      [10] Tuan D Pham,
%           "Time-shift multiscale entropy analysis of physiological 
%           signals." 
%           Entropy 
%           19.6 (2017): 257.
% 
%      [11] Hamed Azami and Javier Escudero,
%           "Coarse-graining approaches in univariate multiscale sample 
%           and dispersion entropy." 
%           Entropy 
%           20.2 (2018): 138.
%
%      [12] Magdalena Costa and Ary Goldberger, 
%           "Generalized multiscale entropy analysis: Application to
%           quantifying the complex volatility of human heartbeat time
%           series"
%           Entropy
%           17 (2015): 1197–1203
% 
narginchk(2,10)
Sig = squeeze(Sig);
p = inputParser;
Chk = {'coarse';'generalized';'modified';'imf';'timeshift'};
addRequired(p,'Sig',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addRequired(p,'Mobj',@(x) isstruct(x));
addParameter(p,'Methodx','coarse',@(x) ...
    any(validatestring(string(lower(x)),Chk)));
addParameter(p,'Scales',3,@(x) isnumeric(x) && (length(x)==1) && (x>1));
addParameter(p,'Plotx',false,@(x) islogical(x));
addParameter(p,'RadNew',0,@(x) x==0 || (ismember(x,1:4) && ...
    any(validatestring(func2str(Mobj.Func),{'SampEn';'ApEn'}))));
parse(p,Sig, Mobj, varargin{:})
Methodx = str2func(lower(p.Results.Methodx));
MSx = zeros(1,p.Results.Scales);
RadNew = p.Results.RadNew;
if strcmp(func2str(Mobj.Func),'SampEn'),  Mobj.Vcp = false; end
Fields = fieldnames(Mobj);
Y = struct2cell(Mobj);
C = [Fields(2:end),Y(2:end)].';
assert(~startsWith(lower(func2str(Mobj.Func)), 'x'),...
    "Base entropy estimator is a cross-entropy method. To perform multiscale CROSS-entropy estimation, use XMSEn.")
if strcmpi(p.Results.Methodx,'imf')
    [Imfx,Resx] = emd(Sig,'MaxNumIMF',p.Results.Scales-1,'Display',0); 
    Sig = [Imfx Resx]';
    clear Imfx Resx
end
if RadNew    
    switch RadNew
        case 1
            Rnew = @(x) std(x,1);
        case 2
            Rnew = @(x) var(x,1);
        case 3            
            Rnew = @(x) mad(x);
        case 4
            Rnew = @(x) mad(x,1);
    end
    
    try
        C_Loc = find(strcmp(C(1,:),'r'));
        Cx = C{2,C_Loc};
    catch
        Cy = {'Standard Deviation';'Variance';...
            'Mean Abs Deviation';'Median Abs Deviation'};
        warning(['No radius value provided.\n' ...
            'Default set to 0.2*(%s) of each new time-series.'],Cy{RadNew})
        %  C_Loc = length(C(1,:)) + 1;
        C_Loc = size(C,1) + 1;
        
        C{1,C_Loc} = 'r';
        Cx = .2;
    end
end
for T = 1:p.Results.Scales       
    fprintf(' .')
    Temp = Methodx(Sig,T);  
    if strcmpi(p.Results.Methodx,'timeshift')
        Tempx = zeros(1,T);
        for k = 1:T
            if RadNew
                C{2,C_Loc} = Cx*Rnew(Temp(k,:));
            end
            Tempy = Mobj.Func(Temp(k,:),C{:});  
            Tempx(k) = Tempy(end);
        end
        Temp2 = mean(Tempx);
        clear Tempx Tempy
    else
        if RadNew
            C{2,C_Loc} = Cx*Rnew(Temp);
        end
        Temp2 = Mobj.Func(Temp,C{:});
    end
    MSx(T) = Temp2(end);  
    clear Temp Temp2    
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
   title(sprintf('Multiscale %s (%s-graining method)',func2str(Y{1}),...
       func2str(Methodx)),'FontSize',16,'FontWeight','bold','Color',[7 54 66]/255)
   ylim([0 max(MSx)+.2])
end
end
function [Y] = generalized(Z,sx)
    Ns = floor(length(Z)/sx);
    Y = var(reshape(Z(1:sx*Ns),sx,Ns),1,1);
end
function [Y] = coarse(Z,sx)
    Ns = floor(length(Z)/sx);
    Y = mean(reshape(Z(1:sx*Ns),sx,Ns),1);
end
function [Y] = modified(Z,sx)
    % Y = movmean(Z,sx);
    % Y = Y(ceil((sx+1)/2):end-ceil(sx/2)+1);
    Y = conv(Z,ones(1,sx),'valid')/sx;
end
function [Y] = imf(Z,sx)
    Y = squeeze(sum(Z(1:sx,:),1));
end
function [Y] = timeshift(Z,sx)
    Y = reshape(Z(1:sx*floor(length(Z)/sx)),sx,floor(length(Z)/sx));
end
%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub