function [sizeI_ddra, contain_pct] = ddra_infer_size_streaming(sys, R0, U_unused, W, M_AB, C, VAL)
% Like ddra_infer, but never stores Xsets; accumulates output-size proxy &
% synchronized containment on the exact validation points in VAL.

    % ---- guards ----
    assert(nargin>=7 && ~isempty(VAL) && isfield(VAL,'u') && isfield(VAL,'y'), ...
        'ddra_infer_size_streaming: VAL (y,u,...) is required for synchronized validation.');

    % dims and options
    nx = size(sys.A,1);
    m  = size(sys.B,2);

    LM    = getfielddef(C,'lowmem',struct());
    ordCap= getfielddef(LM,'zonotopeOrder_cap', ...
                 getfielddef(C.shared.options_reach,'zonotopeOrder',100));
    Kred  = min(100, ordCap);

    % ----- pick C, D robustly (avoid isprop on CORA classes) -----
    Cmat = [];
    Dmat = [];

    if isfield(VAL,'C') && ~isempty(VAL.C), Cmat = VAL.C; end
    if isempty(Cmat)
        Cmat = safe_get_prop(sys,'C',[]);
    end
    if isempty(Cmat)
        % no output matrix provided -> identity output
        Cmat = eye(nx);
    end
    ny = size(Cmat,1);

    if isfield(VAL,'D') && ~isempty(VAL.D), Dmat = VAL.D; end
    if isempty(Dmat)
        Dsys = safe_get_prop(sys,'D',[]);
        if isempty(Dsys), Dmat = zeros(ny,m); else, Dmat = Dsys; end
    end

    sizeI_ddra = 0;
    num_in = 0; num_all = 0;
    tol = 1e-6;

    B = numel(VAL.x0);
    for b = 1:B
        x0 = VAL.x0{b};                 % (nx×1)
        Uk = VAL.u{b};                  % (n_k×m) or (m×n_k)
        Yk = VAL.y{b};                  % (n_k×ny)

        % normalize input to (n_k×m)
        if size(Uk,2) ~= m && size(Uk,1) == m, Uk = Uk.'; end
        n_k = size(Uk,1);

        % start set at R0 shifted by this x0
        X = R0 + x0;

        for k = 1:n_k
            X  = reduce(X,'girard',Kred);
            u_k = Uk(k,:).';            % (m×1)
            U_point = zonotope(u_k);    % deterministic input
            Xnext = M_AB * cartProd(X, U_point) + W;

            % output set and its interval size
            Yset = Cmat*Xnext + Dmat*u_k;
            Iv   = interval(Yset);
            sizeI_ddra = sizeI_ddra + sum(abs(Iv.sup(:) - Iv.inf(:)));

            % synchronized containment
            y_meas = Yk(k,:).';
            if contains_interval(y_meas, Yset, tol), num_in = num_in + 1; end
            num_all = num_all + 1;

            X = Xnext;
        end
    end

    contain_pct = (num_all > 0) * (100 * num_in / max(1,num_all));
    if num_all == 0, contain_pct = NaN; end
end

% ---------- helpers ----------
function val = safe_get_prop(obj, pname, defaultVal)
    try
        val = obj.(pname);
    catch
        val = defaultVal;
    end
end

function val = getfielddef(S, fname, defaultVal)
    if isstruct(S) && isfield(S, fname); val = S.(fname); else; val = defaultVal; end
end
