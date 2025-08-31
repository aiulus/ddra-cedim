function [M_AB, ridgeInfo, W_out] = ddra_learn_Mab(Xminus, Uminus, Xplus, W_in, Zinfo, sys, C)
% DDRA_LEARN_MAB  Learn matrix-zonotope over [A B] with optional ridge guard.
% If Z is rank-deficient:
%   - allow_ridge = false  -> mark skip (return []), caller should skip this grid point
%   - allow_ridge = true   -> ridge with λ, then widen either M_AB (extra generator)
%                             or W (extra generator) proportional to λ*||Z^\dagger||_2.

    if nargin < 7, C = struct(); end
    allow_ridge = getfielddef(C,'ddra',struct()); 
    if isstruct(allow_ridge)
        allow_ridge_flag = getfielddef(C.ddra,'allow_ridge',false);
        lambda           = getfielddef(C.ddra,'lambda',1e-8);
        ridge_gamma      = getfielddef(C.ddra,'ridge_gamma',1.0);
        ridge_policy     = string(getfielddef(C.ddra,'ridge_policy',"MAB"));
    else
        allow_ridge_flag = false; lambda = 1e-8; ridge_gamma = 1.0; ridge_policy = "MAB";
    end

    % Build Mw (matrix zonotope of stacked disturbances) exactly like DDRA snippet
    Mw = build_Mw_matrix_zonotope(W_in, size(sys.A,1), size(Xplus,2));

    % Subtract center, keep generators
    X1W_center = Xplus - Mw.center;
    X1W        = matZonotope(X1W_center, Mw.G);

    % Regressor
    Z        = [Xminus; Uminus];
    n_rows   = size(Z,1);
    n_cols   = size(Z,2);
    fullRank = Zinfo.rankZ >= min(n_rows, n_cols);

    ridgeInfo = struct('used',false,'skipped',false,'lambda',0,'kappa',NaN,'policy',char(ridge_policy));
    W_out = W_in;  % default pass-through

    if fullRank
        % Pseudoinverse as in the ground truth
        M_AB = X1W * pinv(Z);

    else
        % Rank-deficient
        if ~allow_ridge_flag
            warning('DDRA:RankDef','Z rank-deficient (rank=%d, rows=%d, cols=%d). Skipping point (allow_ridge=false).', ...
                    Zinfo.rankZ, n_rows, n_cols);
            M_AB = [];    % signal "skip"
            ridgeInfo.skipped = true;
            return
        end

        % Ridge
        Zpinv_ridge = Z' / (Z*Z' + lambda * eye(n_rows));  % (n_cols x n_rows)
        M_AB = X1W * Zpinv_ridge;

        % Uncertainty inflation proportional to λ * ||Z^†||_2
        kappa = norm(Zpinv_ridge, 2);                      % proxy for ||Z^†||
        eps_scale = ridge_gamma * lambda * kappa;

        ridgeInfo.used    = true;
        ridgeInfo.lambda  = lambda;
        ridgeInfo.kappa   = kappa;

        switch upper(char(ridge_policy))
            case 'MAB'
                % Add one extra generator matrix of uniform magnitude eps_scale
                dim_x   = size(sys.A,1);
                dim_col = size(sys.A,2) + size(sys.B,2);
                Gr      = eps_scale * ones(dim_x, dim_col);
                if isempty(M_AB.G)
                    G3 = Gr;              % (dim_x x dim_col) -> treated as one slice
                else
                    G3 = cat(3, M_AB.G, Gr);
                end
                M_AB = matZonotope(M_AB.center, G3);

            case 'W'
                % Inflate process noise W by adding a generator of magnitude eps_scale
                if isa(W_in,'zonotope')
                    G_old = generators(W_in);
                    if isempty(G_old), G_old = zeros(size(center(W_in),1),0); end
                    G_new = [G_old, eps_scale * ones(size(center(W_in),1),1)];
                    W_out = zonotope(center(W_in), G_new);
                else
                    warning('DDRA:RidgePolicy','W inflation requested but W is not a zonotope; leaving W unchanged.');
                end

                % Rebuild Mw and X1W with inflated W
                Mw2 = build_Mw_matrix_zonotope(W_out, size(sys.A,1), size(Xplus,2));
                X1W2_center = Xplus - Mw2.center;
                X1W2        = matZonotope(X1W2_center, Mw2.G);
                M_AB        = X1W2 * Zpinv_ridge;

            otherwise
                warning('DDRA:RidgePolicy','Unknown ridge_policy "%s"; defaulting to MAB.', ridge_policy);
                dim_x  = size(sys.A,1);
                dim_col = size(sys.A,2) + size(sys.B,2);
                Gr     = eps_scale * ones(dim_x, dim_col);
                if isempty(M_AB.G)
                    G3 = Gr;
                else
                    G3 = cat(3, M_AB.G, Gr);
                end
                M_AB = matZonotope(M_AB.center, G3);
        end

    end

    % Optional sanity check: does [A B] lie within the interval hull?
    boxM   = intervalMatrix(M_AB);
    bounds = boxM.int;
    AB     = [sys.A, sys.B];
    inside = all(all(bounds.inf <= AB & AB <= bounds.sup));
    if ~inside
        warning('DDRA:ContainWarn','[A B] not contained in learned M_{AB} interval hull.');
    end
end
