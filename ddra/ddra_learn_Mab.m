function M_AB = ddra_learn_Mab(Xminus, Uminus, Xplus, W, Zinfo, sys)
    % Build Mw exactly like original DDRA
    Mw = build_Mw_matrix_zonotope(W, size(sys.A,1), size(Xplus,2));

    % Original behavior: subtract center, pass generators as-is
    X1W_center = Xplus - Mw.center;
    X1W        = matZonotope(X1W_center, Mw.generator);

    % Regress M_AB with ridge fallback if rank-deficient
    Z = [Xminus; Uminus];
    full_row = size(Z,1);
    if Zinfo.rankZ < min(full_row, size(Z,2))
        lambda = 1e-8;
        M_AB = X1W * ((Z'/(Z*Z' + lambda*eye(full_row)))');
    else
        M_AB = X1W * pinv(Z);
    end

    % --- Containment check: avoid nested subsref on CORA objects ---
    boxM   = intervalMatrix(M_AB);
    bounds = boxM.int;                 % one level of subsref only
    AB     = [sys.A, sys.B];

    % If your MATLAB supports it:
    % inside = all(bounds.inf <= AB & AB <= bounds.sup,'all');
    % Backward-compatible:
    inside = all(all(bounds.inf <= AB & AB <= bounds.sup));

    if ~inside
        warning('DDRA: [A B] not contained in M_AB interval.');
    end
end
