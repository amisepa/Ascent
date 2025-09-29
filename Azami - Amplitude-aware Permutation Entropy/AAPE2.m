function [Out_AAPE hist] = AAPE2(y,m,varargin)


% Permutation entropy (PE) has two key shortcomings. First, when a signal is symbolized using the Bandt-Pompe procedure,
% only the order of the amplitude  values are considered and information regarding the amplitudes is discarded.
% Second, in the PE, the effect of equal amplitude values in each embedded vector is not addressed.
% To address these issues, we proposed a new entropy measure based on PE: the amplitude-aware permutation entropy (AAPE)
% For more information, please see reference [1].
%
%
% AAPE1: address the first above mentioned problem
% AAPE2: address the second above mentioned problem
% AAPE: address both of the above mentioned problems
%
% Input:
% y: univariate signal
% m: order of permuation entropy
% t: delay time of permuation entropy,
%
% Output:
%Out_AAPE: amplitude-aware permutation entropy (AAPE)
%
%
% Ref: [1] H. Azami and J. Escudero, Amplitude-aware Permutation Entropy: Illustration in Spike Detection
% and Signal Segmentation Computer Methods and Programs in Biomedicine, 2016.     %
%
%
% If you use the code, please make sure that you cite reference [1].
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  10-February-16


if length(varargin) == 1;
    t = varargin{1};
else
    t = 1;
end

ly = length(y);
permlist = perms(1:m);
c(1:length(permlist))=0;


for j=1:ly-t*(m-1)
    [a,iv]=sort(y(j:t:j+t*(m-1)));
    
    if 0~=min(abs(diff(a)))  % in case the embedding vectors do not have any equal values
        for jj=1:length(permlist)
            if (abs(permlist(jj,:)-iv))==0
                c(jj) = c(jj) + 1;
            end
        end
    else % in case the embedding vectors have equal values (addressing the second above mentioned problem of PE)
        
        [u1 u2 u3]=unique(y(j:t:j+t*(m-1)));
        TF=1;
        for ip=1:length(u1)
            f1=find(u3==ip);
            TF=factorial(length(f1))*TF;
        end
        
        yy=NaN*ones(TF,m);
        
        for ip=1:length(u1)
            
            a=find(u3==ip);
            pa=perms(a);
            yy(:,a)=repmat(pa,TF/(factorial(length(a))),1);
        end
        
        for jjj=1:TF
            iv=yy(jjj,:);
            for jj=1:length(permlist)
                if (abs(permlist(jj,:)-iv))==0
                    c(jj) = c(jj) + 1;
                end
            end
        end
        
    end
    
end
hist = c;

c=c(find(c~=0));
p = c/sum(c);
Out_AAPE = -sum(p .* log(p));
