function c_pct = mc_containment_X(sys, R0, U, W, Xsets, C)
% Monte-Carlo fidelity proxy: simulate true states and test if they are inside
% the predicted reachable sets Xsets using an APPROX method (interval hull).
%
% Tunables in C.shared (all optional):
%   .mc_trials        (default 50)    number of Monte Carlo rollouts
%   .mc_steps_max     (default  min(10, numel(Xsets))) steps per trial to evaluate
%   .mc_reduce_order  (default  40)   reduce() order cap before testing
%   .mc_tol           (default 1e-8)  interval slack
%   .mc_method        (default 'interval')  'interval' | 'none'
%
% Notes:
%   - We avoid CORA's contains() to bypass polytope conversion.
%   - 'interval' test: x in [c - sum|G|, c + sum|G|] (O(nGen), no combinatorics).
%   - This is a fidelity *proxy*; containment may be slightly optimistic.

    % Defaults
    trials   = getfielddef(C.shared,'mc_trials',50);
    Ktotal   = numel(Xsets);
    Kcheck   = min(getfielddef(C.shared,'mc_steps_max',min(10,Ktotal)), Ktotal);
    ord_cap  = getfielddef(C.shared,'mc_reduce_order',40);
    tol      = getfielddef(C.shared,'mc_tol',1e-8);
    method   = getfielddef(C.shared,'mc_method','interval');

    % Pre-reduce predicted sets to keep membership cheap
    Zcap = cell(1,Kcheck);
    for k = 1:Kcheck
        Zk = Xsets{k};
        if ~isa(Zk,'zonotope'); Zk = zonotope(Zk); end
        if size(Zk.generators,2) > ord_cap
            Zk = reduce(Zk,'girard',ord_cap);
        end
        Zcap{k} = Zk;
    end

    contain = 0; total = 0;
    for s = 1:trials
        % Sample initial state from R0
        x = randPoint(R0);
        for k = 1:Kcheck
            % Simulate one step of the true system with random admissible (u,w)
            u = randPoint(U);
            w = randPoint(W);
            x = sys.A * x + sys.B * u + w;

            Zk = Zcap{k};

            % ---- APPROXIMATE membership (interval hull) ----
            % interval(Z) = [c - sum|G|, c + sum|G|]
            I = interval(Zk);
            inside = all(x >= I.inf - tol & x <= I.sup + tol);

            contain = contain + double(inside);
            total   = total + 1;
        end
    end
    c_pct = 100 * (contain / max(total,1));
end
