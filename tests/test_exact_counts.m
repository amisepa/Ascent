classdef test_exact_counts < matlab.unittest.TestCase
% Exact pair counts and sums against all-pairs references.
%   sampen_exact (templates sorted by their first coordinate, or all pairs for
%   short series and when many are within r on it) against counting
%   pdist(..., 'chebychev') <= r over all pairs, on each route. Counts are integers, so they must be equal, and so
%   must SampEn. The signals include ties (integer data,
%   a constant series), distances exactly equal to r, decimal steps whose
%   differences round on either side of r, non-finite values, heavy tails, lags
%   above 1 and series too short for two templates.
%   The blocked fuzzy sum (used above pdistMaxGB) takes its across-block tiles
%   from pdist2: each tile must equal the maximum over coordinates of
%   |Vi(:,k) - Vj(:,k)'| exactly, and the sum must equal the pdist path to
%   round-off (the order of summation differs).
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).

    methods (TestClassSetup)
        function addPaths(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'functions')));
        end
    end

    methods (Test)
        function countsMatchAllPairs(tc)
            rng(1);
            n = 600;
            g = randn(n, 1);
            ar = filter(1, [1 -.95], randn(n, 1));
            nf = g; nf([5 77 300]) = [NaN Inf -Inf];
            S = {
                'gauss',    g,                                  [.05 .15 .5 2]
                'ar1',      ar / std(ar),                       [.05 .15 .5]
                'sine',     sin(2*pi*.05*(1:n)') + .1*randn(n, 1), [.15 .5]
                'integer',  randi([-3 3], n, 1),                [1 2]
                'quarter',  round(4*randn(n, 1)) / 4,           [.25 .5]
                'decimal',  round(30*rand(n, 1)) / 10,          [.1 .2 .3]
                'constant', 3*ones(80, 1),                      [.15]
                'cauchy',   randn(n, 1) ./ randn(n, 1),         [.15 1]
                'nonfinite', nf,                                [.15]
                };
            nCase = 0;
            for i = 1:size(S, 1)
                x = S{i, 2};
                for r = S{i, 3}
                    for m = 1:4
                        for tau = [1 2 3 5]
                            [seRef, cntRef] = ref_sampen(x, m, r, tau);
                            for route = {'sweep', 'pairs', 'auto'}
                                [se, cnt] = sampen_exact(x, m, r, tau, route{1});
                                lbl = sprintf('%s m=%d tau=%d r=%g %s', S{i, 1}, m, tau, r, route{1});
                                tc.verifyEqual(cnt, cntRef, lbl);
                                tc.verifyTrue(isequaln(se, seRef), lbl);
                            end
                            nCase = nCase + 1;
                        end
                    end
                end
            end
            tc.verifyEqual(nCase, 320);
        end

        function shortSeries(tc)
            rng(2);
            for N = 1:9
                x = randn(N, 1);
                for m = 1:3
                    for tau = 1:3
                        [seRef, cntRef] = ref_sampen(x, m, 1, tau);
                        for route = {'sweep', 'pairs'}
                            [se, cnt] = sampen_exact(x, m, 1, tau, route{1});
                            lbl = sprintf('N=%d m=%d tau=%d %s', N, m, tau, route{1});
                            tc.verifyEqual(cnt, cntRef, lbl);
                            tc.verifyTrue(isequaln(se, seRef), lbl);
                        end
                    end
                end
            end
        end

        function longSeries(tc)
            % 8000 samples through pdist; 25000 through the blocked counter
            % (the old route above 2 GB of pairs)
            rng(3);
            x = filter(1, [1 -.9], randn(8000, 1)); x = x / std(x);
            [~, cnt] = sampen_exact(x, 2, .15, 1, 'sweep');
            [~, cntRef] = ref_sampen(x, 2, .15, 1);
            tc.verifyEqual(cnt, cntRef);
            x = randn(25000, 1);
            m = 3; tau = 2; r = .3;
            [~, cnt] = sampen_exact(x, m, r, tau, 'auto');   % beyond 2 GB of pairs: the sweep
            Xm  = x((1:numel(x)-(m-1)*tau)' + (0:m-1)*tau);
            Xm1 = x((1:numel(x)-m*tau)' + (0:m)*tau);
            tc.verifyEqual(cnt([3 1]), [pairwise_match_count_cheby_blocked(Xm, r, 2000) ...
                pairwise_match_count_cheby_blocked(Xm1, r, 2000)]);
        end

        function computeSampEnUsesExactCounts(tc)
            rng(4);
            X = randn(3, 900) .* [1; 5; .1] + [0; 2; -1];
            [se, cnt] = compute_SampEn(X, 'm', 2, 'r', .2, 'tau', 2, 'Parallel', false, 'Progress', false);
            for c = 1:3
                z = (X(c, :) - mean(X(c, :))) ./ std(X(c, :));
                [seRef, cntRef] = ref_sampen(z, 2, .2, 2);
                tc.verifyEqual(cnt(c, :), cntRef);
                tc.verifyTrue(isequaln(se(c), seRef));
            end
        end

        function blockedTileIsExact(tc)
            rng(5);
            for d = [1 2 3 8 17]
                Vi = round(randn(300, d) * 8) / 8;
                Vj = round(randn(250, d) * 8) / 8;
                D = zeros(300, 250);
                for k = 1:d
                    D = max(D, abs(Vi(:, k) - Vj(:, k)'));
                end
                tc.verifyEqual(pdist2(Vi, Vj, 'chebychev'), D, sprintf('dim %d', d));
            end
        end

        function blockedSumMatchesPdist(tc)
            rng(6);
            K = {'exponential', 2; 'exponential', 3; 'gaussian', 2};
            for d = [1 3 8]
                V = randn(1037, d);
                for i = 1:size(K, 1)
                    a = fuzzy_pairmean_cheby_exact(V, .2, K{i, 2}, 100, 1e-9, K{i, 1});
                    b = fuzzy_pairmean_cheby_exact(V, .2, K{i, 2}, 100, 2, K{i, 1});
                    tc.verifyLessThan(abs(a - b) / b, 1e-12, sprintf('dim %d %s n=%d', d, K{i, 1}, K{i, 2}));
                end
            end
        end
    end
end


function [se, cnt] = ref_sampen(x, m, r, tau)
% SampEn from all pairs (pdist), as compute_SampEn computed it before sampen_exact
x = x(isfinite(x)); x = x(:);
N = numel(x);
se = NaN; cnt = zeros(1, 4);
if N <= m + 1, return, end
Lm = N - (m-1)*tau; Lm1 = N - m*tau;
if Lm < 2 || Lm1 < 2, return, end
hB = nnz(pdist(x((1:Lm)' + (0:m-1)*tau), 'chebychev') <= r);
hA = nnz(pdist(x((1:Lm1)' + (0:m)*tau), 'chebychev') <= r);
pB = Lm*(Lm-1)/2; pA = Lm1*(Lm1-1)/2;
cnt = [hA pA hB pB];
B = hB / pB; A = hA / pA;
if B > 0 && A > 0 && isfinite(B) && isfinite(A)
    se = -log(A / B);
end
end
