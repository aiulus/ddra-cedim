function [sizeI_ddra, contain_pct, wid_ddra] = ...
    ddra_infer_size_streaming(sys, R0, U_unused, W_in, M_AB, C, VAL)

    assert(nargin>=7 && ~isempty(VAL) && isfield(VAL,'u') && isfield(VAL,'y'), ...
        'ddra_infer_size_streaming: VAL (y,u,...) is required for synchronized validation.');

    nx = size(sys.A,1);
    m  = size(sys.B,2);

    LM     = getfielddef(C,'lowmem',struct());
    ordCap = getfielddef(LM,'zonotopeOrder_cap', ...
                 getfielddef(C.shared.options_reach,'zonotopeOrder',100));
    Kred   = min(100, ordCap);

    cR = center(R0);
    if any(abs(cR) > 1e-12)
        R0 = zonotope(zeros(nx,1), R0.G);   % strip center; x0 already carries it
    end


    % Robust C, D (output map)
    Cmat = []; Dmat = [];
    if isfield(VAL,'C') && ~isempty(VAL.C), Cmat = VAL.C; end
    if isempty(Cmat), Cmat = safe_get_prop(sys,'C',[]); end
    if isempty(Cmat), Cmat = eye(nx); end
    ny = size(Cmat,1);

    if isfield(VAL,'D') && ~isempty(VAL.D), Dmat = VAL.D; end
    if isempty(Dmat)
        Dsys = safe_get_prop(sys,'D',[]);
        if isempty(Dsys), Dmat = zeros(ny,m); else, Dmat = Dsys; end
    end

    % Disturbance accessor
    useSeq = iscell(W_in);

    sizeI_ddra = 0; wid_ddra = [];
    num_in = 0; num_all = 0;
    tol = 1e-6;

    B = numel(VAL.x0);
    for b = 1:B
        x0 = VAL.x0{b}; if isrow(x0), x0 = x0.'; end
        Uk = VAL.u{b};  if size(Uk,2) ~= m && size(Uk,1) == m, Uk = Uk.'; end
        Yk = VAL.y{b};

        assert(size(Yk,2) == ny, 'VAL.y has %d outputs, model has %d.', size(Yk,2), ny);
        assert(size(Uk,2) == m,  'VAL.u has %d inputs, model has %d.',   size(Uk,2), m);
        n_k = size(Uk,1);

        if useSeq
            assert(numel(W_in) >= n_k, 'W sequence shorter than n_k');
        end

        % state set at k=0
        X = R0 + x0;

        for k = 1:n_k
            X  = reduce(X,'girard',Kred);           % keep sets compact

            % deterministic input at time k
            u_k   = Uk(k,:).';
            U_pt  = zonotope(u_k);

            % --- OUTPUT at time k uses PRE-UPDATE state X ---
            Yset = Cmat*X + Dmat*u_k;

            % interval-hull width as size proxy
            IY   = interval(Yset);
            w    = supremum(IY) - infimum(IY);                         % 1×ny widths
            v_k  = sum(w);                           % aggregate width
            sizeI_ddra      = sizeI_ddra + v_k;
            wid_ddra(end+1,1)= v_k;

            % containment in interval hull of Yset
            y_meas = Yk(k,:).';
            if contains_interval(y_meas, Yset, tol)
                num_in = num_in + 1;
            end
            num_all = num_all + 1;

            % --- propagate to next state (post-update) ---
            if useSeq
                Wk = W_in{k};
            else
                Wk = W_in;
            end

            % sanity check: disturbance dim must match state dim
            if ~isempty(Wk)
                assert(size(center(Wk),1) == nx, 'W dim %d ≠ nx %d', size(center(Wk),1), nx);
            end
            if iscell(W_in)
                warning('Algorithm 1 expects constant Zw; using W_in{1} for all k.');
                W_const = W_in{1};
            else
                W_const = W_in;
            end
            Xnext = M_AB * cartProd(X, U_pt) + W_const;     % DDRA-LTI: x_{k+1} ∈ M_AB * [x_k;u_k] \oplus W_k
            Xnext = reduce(Xnext, 'girard', Kred);
            X = Xnext;
        end
    end

    contain_pct = (num_all > 0) * (100 * num_in / max(1,num_all));
    if num_all == 0, contain_pct = NaN; end

    % optional: return average size 
    sizeI_ddra = sizeI_ddra / max(1,num_all);
end

% ---------- helpers ----------
function val = safe_get_prop(obj, pname, defaultVal)
    try
        val = obj.(pname);
        if isempty(val), val = defaultVal; end
    catch
        val = defaultVal;
    end
end

function val = getfielddef(S, fname, defaultVal)
    if isstruct(S) && isfield(S, fname); val = S.(fname); else; val = defaultVal; end
end
