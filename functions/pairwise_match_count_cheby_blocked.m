function hits = pairwise_match_count_cheby_blocked(V, r, blockSize)

m = size(V,1);
pdim = size(V,2);

hits = 0;

for i0 = 1:blockSize:m-1
    i1 = min(i0 + blockSize - 1, m-1);
    Vi = V(i0:i1,:);
    nb = size(Vi,1);

    % within-block
    for ii = 1:nb-1
        D = max(abs(Vi(ii+1:end,:) - Vi(ii,:)), [], 2);
        hits = hits + nnz(D <= r);
    end

    % across blocks
    for j0 = i1+1:blockSize:m
        j1 = min(j0 + blockSize - 1, m);
        Vj = V(j0:j1,:);
        nj = size(Vj,1);

        Dmax = zeros(nb, nj);
        for k = 1:pdim
            Dmax = max(Dmax, abs(Vi(:,k) - Vj(:,k)'));
        end

        hits = hits + nnz(Dmax <= r);
    end
end
end