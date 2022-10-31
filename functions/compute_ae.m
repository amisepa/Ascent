% Computes the approximate entropy (AE)
%
% INPUTS
%   signal : time-series data
%   m    : Template length
%   r    : tolerance (typically 0.2 * std)
% OUTPUT
%   ae : approximate entropy
%
% Cedric Cannard


function ae = compute_ae(signal, m, r)

N = length(signal);
result = zeros(1,2);

for j = 1:2
    m = m+j-1;
    phi = zeros(1,N-m+1);
    dataMat = zeros(m,N-m+1);
    
    % setting up data matrix
    for i = 1:m
        dataMat(i,:) = signal(i:N-m+i);
    end
    
    % counting similar patterns using distance calculation
    for i = 1:N-m+1
        
        if sum( isnan( dataMat(:,i) ) ) == 0
        
            tempMat = abs(dataMat - repmat(dataMat(:,i),1,N-m+1));
            boolMat = any( (tempMat > r),1);
            phi(i) = sum(~boolMat)/(N-m+1);
        else

	        %discarding blocks with nan values
            phi(i) = nan;
        end
    end
    
    % summing over the counts
    result(j) = log(mean(phi, 'omitnan'));    
end

ae = result(1) - result(2) ;

end