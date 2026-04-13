function p = match_fraction_cheby_exact(V, r, blockSize, pdistMaxGB)
% Exact fraction of pairs with Chebyshev distance <= r

M = size(V,1);
if M < 2
    p = NaN;
    return
end

pairs = double(M) * double(M-1) / 2;
gbNeeded = (pairs * 8) / (1024^3);

if gbNeeded <= pdistMaxGB
    d = pdist(V, 'chebychev');
    p = mean(d <= r);
else
    hits = pairwise_match_count_cheby_blocked(V, r, blockSize);
    p = hits / pairs;
end
end