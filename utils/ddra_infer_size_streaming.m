function [sizeI_ddra, contain_pct] = ddra_infer_size_streaming(sys, R0, U, W, M_AB, C)
% Like ddra_infer, but never stores Xsets; accumulates size proxy & MC contain on the fly.

    n_k_val = C.shared.n_k_val;
    X = R0;
    sizeI_ddra = 0;
    contain = 0; total = 0;

    % small MC batch per step
    Nmc = getfielddef(C,'mc_batch', 50);

    % cap order for memory if present
    zOrd = getfielddef(C.shared.options_reach,'zonotopeOrder',100);
    zOrd = min(100, zOrd);

    for k = 1:n_k_val
        % propagate one step with M_AB
        X = reduce(X, 'girard', zOrd);
        Xnext = M_AB * cartProd(X, U) + W;

        % --- interval size proxy (works for both matZonotope and zonotope)
        [lo, hi] = local_interval_bounds(Xnext);
        sizeI_ddra = sizeI_ddra + sum(abs(hi(:) - lo(:)));

        % --- on-the-fly MC containment (single-step)
        for i = 1:Nmc
            x0 = randPoint(R0); u = randPoint(U); w = randPoint(W);
            x1 = sys.A*x0 + sys.B*u + w;
            if contains(Xnext, x1, 'approx', 1e-6), contain = contain + 1; end
            total = total + 1;
        end

        X = Xnext;  % step forward
    end
    contain_pct = 100 * contain / max(1, total);
end

% --- helpers --------------------------------------------------------------
function [lo, hi] = local_interval_bounds(S)
    % matZonotope -> use intervalMatrix, then pull .int.inf/.sup
    if isa(S, 'matZonotope')
        BM = intervalMatrix(S);      % struct with .int.inf/.sup
        lo = BM.int.inf; hi = BM.int.sup;
        return
    end

    % zonotope / polyZonotope / interval -> use interval(S)
    try
        IV = interval(S);            % returns object/struct with .inf/.sup
        lo = IV.inf; hi = IV.sup;
    catch
        % some CORA versions also support [lo,hi] = interval(S)
        try
            [lo, hi] = interval(S);
        catch ME
            error('Could not obtain interval bounds for object of class %s: %s', class(S), ME.message);
        end
    end
end
