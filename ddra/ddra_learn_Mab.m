function M_AB = ddra_learn_Mab(Xminus, Uminus, Xplus, W, Zinfo, sys)
    % NEW: zero-noise fast path (skip Mw entirely)
    if ~hasGenerators(W)
        % X1W is exactly Xplus; no matrix-zonotope generators
        X1W_center = Xplus;
        X1W = matZonotope(X1W_center, zeros(size(Xplus,1), size(Xplus,2), 0));
    else
        % Original behavior
        Mw = build_Mw_matrix_zonotope(W, size(sys.A,1), size(Xplus,2));
        X1W_center = Xplus - Mw.center;
        X1W        = matZonotope(X1W_center, Mw.generator);
    end

    % Regress M_AB with ridge fallback if rank-deficient
    Z = [Xminus; Uminus];
    full_row = size(Z,1);
    if Zinfo.rankZ < min(full_row, size(Z,2))
        lambda = 1e-8;
        M_AB = X1W * ((Z'/(Z*Z' + lambda*eye(full_row)))');
    else
        M_AB = X1W * pinv(Z);
    end

    % Containment check 
    boxM   = intervalMatrix(M_AB);
    bounds = boxM.int;
    AB     = [sys.A, sys.B];
    inside = all(all(bounds.inf <= AB & AB <= bounds.sup));
    if ~inside
        warning('DDRA: [A B] not contained in M_AB interval.');
    end
end
