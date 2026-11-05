function [se, cnt] = sampen_exact(x, m, r, tau, route)
% sampen_exact  Sample entropy of one series from exact match counts.
%
%   [se, cnt] = sampen_exact(x, m, r, tau)
%
%   x   : vector (non-finite values are dropped)
%   m   : embedding dimension
%   r   : absolute Chebyshev tolerance (a pair matches when every coordinate
%         differs by <= r)
%   tau : time lag
%
%   cnt = [matches(m+1) pairs(m+1) matches(m) pairs(m)], over the N-m*tau
%   templates of length m+1 and the N-(m-1)*tau templates of length m, and
%   se = -log((cnt(1)/cnt(2)) / (cnt(3)/cnt(4))) (NaN without matches).
%
%   The counts are those of comparing every pair of templates, found without
%   forming the pairs. Templates are sorted by their first coordinate, so a
%   template can only match the w(p) templates that follow it in that order
%   (those whose first coordinate is within r; floating-point subtraction is
%   monotone, so they form a run). Each lag k in sorted order is then checked
%   for all templates at once on the other coordinates. Memory is O(N).
%   Comparing all pairs in compiled code (pdist) is faster when many pairs are
%   within r on the first coordinate: on one thread the two took the same time
%   at a share of .4-.5 (4000-12000 samples, M1 Max), and the sweep was 5x
%   faster at the usual .03-.1 (and faster from 100 samples up). pdist uses
%   every computational thread and the sweep one, so the share is divided by
%   maxNumCompThreads (parallel workers have one). pdist is used only when the
%   pairs fit in 2 GB. Counts are integers: both routes give the counts of
%   pdist(..., 'chebychev') <= r exactly.
%   route (optional, for tests): 'auto' (default), 'sweep' or 'pairs'.
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).

if nargin < 5 || isempty(route), route = 'auto'; end
cnt = zeros(1, 4);
se  = NaN;
x   = x(isfinite(x));
x   = x(:);
N   = numel(x);
if N <= m + 1
    return
end
Lm  = N - (m-1)*tau;            % templates of length m
Lm1 = N - m*tau;                % templates of length m+1 (the first Lm1 of them)
if Lm < 2 || Lm1 < 2
    return
end
pB = double(Lm)  * double(Lm - 1)  / 2;
pA = double(Lm1) * double(Lm1 - 1) / 2;
fits = pB * 8 <= 2 * 1024^3;

if strcmp(route, 'pairs')
    [hB, hA] = all_pairs(x, Lm, Lm1, m, r, tau);
    [se, cnt] = finish(hA, pA, hB, pB);
    return
end

[xs, ord] = sort(x(1:Lm));

% w(p): number of templates after p (in sorted order) whose first
% coordinate is within r of xs(p)
hi = zeros(Lm, 1);
q = 1;
for p = 1:Lm
    if q < p, q = p; end
    while q < Lm && xs(q+1) - xs(p) <= r
        q = q + 1;
    end
    hi(p) = q;
end
w = hi - (1:Lm).';

if strcmp(route, 'auto') && fits && sum(w) > .4 / maxNumCompThreads * pB
    [hB, hA] = all_pairs(x, Lm, Lm1, m, r, tau);
else
    [hB, hA] = sweep(x, ord, w, N, Lm1, m, r, tau);
end
[se, cnt] = finish(hA, pA, hB, pB);
end


function [se, cnt] = finish(hA, pA, hB, pB)
cnt = [hA pA hB pB];
se  = NaN;
Bm  = hB / pB;
Am1 = hA / pA;
if Bm > 0 && Am1 > 0 && isfinite(Bm) && isfinite(Am1)
    se = -log(Am1 / Bm);
end
end


function [hB, hA] = all_pairs(x, Lm, Lm1, m, r, tau)
hB = nnz(pdist(x((1:Lm).'  + (0:m-1)*tau), 'chebychev') <= r);
hA = nnz(pdist(x((1:Lm1).' + (0:m)*tau),   'chebychev') <= r);
end


function [hB, hA] = sweep(x, ord, w, N, Lm1, m, r, tau)
% Match counts over the pairs (p, p+k) of sorted positions, k <= w(p)
Lm = numel(ord);

% coordinates of each template in sorted order; the (m+1)-th coordinate of
% templates that have none in dimension m+1 is a placeholder, masked by ok1
Ys = zeros(Lm, m+1);
for c = 0:m
    idx = ord + c*tau;
    idx(idx > N) = N;
    Ys(:, c+1) = x(idx);
end
ok1 = ord <= Lm1;

% first and last sorted position with w >= k, for every lag k
kmax = max(w);
first = zeros(kmax, 1);
last  = zeros(kmax, 1);
cmL = cummax(w);
cmR = cummax(w, 'reverse');
p = 1;
for k = 1:kmax
    while cmL(p) < k, p = p + 1; end
    first(k) = p;
end
p = Lm;
for k = 1:kmax
    while cmR(p) < k, p = p - 1; end
    last(k) = p;
end

hB = 0;
hA = 0;
for k = 1:kmax
    P  = (first(k):last(k)).';
    mt = w(P) >= k;
    for c = 2:m
        mt = mt & abs(Ys(P, c) - Ys(P + k, c)) <= r;
    end
    hB = hB + nnz(mt);
    mt = mt & ok1(P) & ok1(P + k) & abs(Ys(P, m+1) - Ys(P + k, m+1)) <= r;
    hA = hA + nnz(mt);
end
end
