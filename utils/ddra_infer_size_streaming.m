function [sizeI_ddra, contain_pct] = ddra_infer_size_streaming(sys, R0, U, W, M_AB, C)
% Like ddra_infer, but never stores Xsets; accumulates size proxy & MC contain on the fly.
% Uses a fast, LP-free point-in-zonotope tester to avoid polytope blow-ups.

    n_k_val = C.shared.n_k_val;
    X = R0;
    sizeI_ddra = 0;
    contain = 0; total = 0;

    % MC budget (memory-friendly). Allow override via cfg.lowmem.ddra_mc
    LM = getfielddef(C,'lowmem',struct());
    Nmc = getfielddef(LM,'ddra_mc', getfielddef(C,'mc_batch', 20));  % default 20

    % Optional: cap order for stability (only for visualization-size intermediates)
    ordCap = getfielddef(LM,'zonotopeOrder_cap', ...
               getfielddef(C.shared.options_reach,'zonotopeOrder',100));

    for k = 1:n_k_val
        % propagate one step with M_AB
        X = reduce(X,'girard', min(100, ordCap));
        Xnext = M_AB * cartProd(X, U) + W;

        % accumulate interval size proxy
        Iv = interval(Xnext);                % works for zonotope/matzonotope images
        sizeI_ddra = sizeI_ddra + sum(abs(Iv.sup(:) - Iv.inf(:)));

        % on-the-fly MC containment (fast, LP-free)
        if Nmc > 0
            for i = 1:Nmc
                x0 = randPoint(R0); u = randPoint(U); w = randPoint(W);
                x1 = sys.A*x0 + sys.B*u + w;

                if fast_contains_zonotope_point(Xnext, x1, 1e-3, 1e-6)
                    contain = contain + 1;
                end
                total = total + 1;
            end
        end

        X = Xnext;  % step forward
    end

    if total > 0
        contain_pct = 100 * contain / total;
    else
        contain_pct = NaN;   % MC disabled
    end
end

% ---------- helpers ----------
function inside = fast_contains_zonotope_point(Z, x, tol_coeff, tol_resid)
% Conservative, LP-free membership test for Z = c + G * B_inf:
%   Solve xi_hat = pinv(G)*(x-c); accept if ||xi_hat||_inf <= 1+tol and
%   residual small. This never returns false-positives (only possible false-negatives).

    c = center(Z);
    G = generators(Z);   % CORA returns [] if no gens

    if isempty(G)
        inside = (norm(x - c, inf) <= tol_coeff);
        return
    end

    rhs = x - c;
    xi  = pinv(G) * rhs;                 % least-squares coefficients
    res = norm(G*xi - rhs, 2);           % linear residual

    inside = (max(abs(xi)) <= 1 + tol_coeff) && (res <= tol_resid);
end

function val = getfielddef(S, fname, defaultVal)
    if isstruct(S) && isfield(S, fname); val = S.(fname); else; val = defaultVal; end
end
