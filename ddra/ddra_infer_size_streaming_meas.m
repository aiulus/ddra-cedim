function [sizeI_ddra, cval_ddra, wid_k] = ddra_infer_size_streaming_meas(sys_true_dt, R0, ~, W_used, ABc, AV_one, V_meas, C, VAL)
%DDRA_INFER_SIZE_STREAMING_MEAS  Per-step widths + containment using
% a measurement-noise aware update:
%   X_{k+1} = ABc * (X_k (+ V_meas for k>1)) × {u_k}  + AV_one + W_used
% Containment uses VAL.y and sys_true_dt.C,D like the standard path.

    % --- reduction cap (match standard path defaults) ---
    ordCap = getfielddef(getfielddef(C,'lowmem',struct()), 'zonotopeOrder_cap', 100);
    Kred   = max(1, min(100, round(ordCap)));

    nkv = C.shared.n_k_val;
    ny  = size(sys_true_dt.C,1);
    Bv  = numel(VAL.x0);

    wid_sums   = zeros(nkv,1);
    wid_counts = zeros(nkv,1);
    num_in = 0; num_all = 0;
    tol = getfielddef(C.metrics,'tol', 1e-6);

    for b = 1:Bv
        Xk = reduce(VAL.R0 + VAL.x0{b}, 'girard', Kred);   % PRE-update X_0
        Uk = VAL.u{b};
        Yb = VAL.y{b};

        for k = 1:nkv
            Xk = reduce(Xk, 'girard', Kred);

            % --- measurement noise enters the *regressor state* for k>1 ---
            if k > 1 && ~isempty(generators(V_meas))
                X_for_AB = Xk + V_meas;
            else
                X_for_AB = Xk;
            end

            uk   = Uk(k,:).';
            U_pt = zonotope(uk);

            % PRE-update output (for width & containment)
            Yset = sys_true_dt.C * Xk + sys_true_dt.D * uk;

            % widths (sum across outputs)
            Iv   = interval(Yset);
            lo   = try_get(@() infimum(Iv), @() Iv.inf);
            hi   = try_get(@() supremum(Iv), @() Iv.sup);
            wvec = max(hi - lo, 0);
            wid_sums(k)   = wid_sums(k)   + sum(double(wvec), 'omitnan');
            wid_counts(k) = wid_counts(k) + 1;

            % containment
            y_meas = Yb(k,:).';
            if contains_interval(y_meas, Yset, tol), num_in = num_in + 1; end
            num_all = num_all + 1;

            % --- POST update
            Xk = (ABc * cartProd(X_for_AB, U_pt)) + AV_one + W_used;
        end
    end

    % aggregate like the standard variant
    wid_k        = wid_sums ./ max(1, wid_counts);
    [sizeI_ddra, wid_k] = normalize_widths(wid_k, ny, "mean");
    cval_ddra    = (num_all>0) * (100 * num_in / max(1,num_all));
end

% ------- tiny utils (local) -------
function v = getfielddef(S, f, d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end
function out = try_get(f_ok, f_fallback)
    try, out = f_ok(); catch, out = f_fallback(); end
end
