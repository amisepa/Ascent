
%% Helpers

function S = pairwise_fuzzy_sum_cheby(V, r, n, blockSize)
% S = sum_{i<j} exp(-(d_ij^n)/r), d_ij Chebyshev distance
m = size(V,1);
p = size(V,2);
S = 0;
for i0 = 1:blockSize:m-1
    i1 = min(i0 + blockSize - 1, m-1);
    Vi = V(i0:i1,:);
    nb = size(Vi,1);
    % within-block upper triangle
    for ii = 1:nb-1
        D = max(abs(Vi(ii+1:end,:) - Vi(ii,:)), [], 2);
        S = S + sum(exp(-(D.^n)/r), 'omitnan');
    end

    % cross blocks
    for j0 = i1+1:blockSize:m
        j1 = min(j0 + blockSize - 1, m);
        Vj = V(j0:j1,:);
        nj = size(Vj,1);

        Dmax = zeros(nb, nj, 'single');
        for k = 1:p
            Dmax = max(Dmax, abs(Vi(:,k) - Vj(:,k)'));
        end

        S = S + sum(exp(-(Dmax.^n)/r), 'all', 'omitnan');
    end
end
end


function A=embd(M,tau,ts)
% This function creates multivariate delay embedded vectors with embedding
% vector parameter M and time lag vector parameter tau.
% M is a row vector [m1 m2 ...mnvar] and tau is also a row vector [tau1 tau2....taunvar] where nvar is the
% number of channels;
% ts is the multivariate time series-a matrix of size nvarxnsamp;
% Ref: M. U. Ahmed and D. P. Mandic, "Multivariate multiscale entropy
% analysis", IEEE Signal Processing Letters, vol. 19, no. 2, pp.91-94.2012
[nvar,nsamp]=size(ts);
A=[];
for j=1:nvar
    for i=1:nsamp-max(M)
        temp1(i,:)=ts(j,i:tau(j):i+M(j)-1);
    end
    A=horzcat(A,temp1);
    temp1=[];
end
end