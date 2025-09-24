function [sizeI_ddra, contain_pct, wid_k] = ...
    ddra_infer_size_streaming(sys, R0, U_unused, W_in, M_AB, C, VAL)
    
    % --- defensive defaults so dot-indexing never hits non-structs ---
    if ~exist('C','var') || ~isstruct(C), C = struct(); end
    if ~isfield(C,'shared') || ~isstruct(C.shared), C.shared = struct(); end
    if ~isfield(C,'metrics') || ~isstruct(C.metrics), C.metrics = struct(); end
    if ~isfield(C.metrics,'width_agg') || isempty(C.metrics.width_agg)
        C.metrics.width_agg = "sum";
    end
    if ~isfield(C,'lowmem') || ~isstruct(C.lowmem), C.lowmem = struct(); end
    if ~isfield(C.lowmem,'zonotopeOrder_cap') || isempty(C.lowmem.zonotopeOrder_cap)
        C.lowmem.zonotopeOrder_cap = 50;
    end


    assert(nargin>=7 && ~isempty(VAL) && isfield(VAL,'u') && isfield(VAL,'y'), ...
        'ddra_infer_size_streaming: VAL (y,u,...) is required for synchronized validation.');

    % --- accumulate per-time-step widths (align with Gray)
    nk_global = max(cellfun(@(Ucell) size(Ucell,1), VAL.u));
    wid_sums   = zeros(nk_global,1);
    wid_counts = zeros(nk_global,1);

    nx = size(sys.A,1);
    m  = size(sys.B,2);

    % Use the lowmem cap only; avoid dot-chaining through C.shared.options_reach
    LM      = getfielddef(C,'lowmem',struct());
    ordCap  = getfielddef(LM,'zonotopeOrder_cap', 100);
    Kred    = max(1, min(100, round(ordCap)));

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
    tol = getfielddef(getfielddef(C,'metrics',struct()),'tol',1e-6);

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

            % interval-hull width proxy: SUM across outputs (canonical internal form)
            IY   = interval(Yset);
            w    = supremum(IY) - infimum(IY);  % 1×ny widths
            v_k  = sum(w);

            % scalar width for this time step (sum of output interval widths)
            sizeI_ddra       = sizeI_ddra + v_k;   % for global scalar (kept)
            wid_sums(k)      = wid_sums(k)   + v_k;
            wid_counts(k)    = wid_counts(k) + 1;


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
            if ~isa(Wk, 'emptySet')
                assert(size(center(Wk),1) == nx, 'W dim %d ≠ nx %d', size(center(Wk),1), nx);
            end

            Xnext = M_AB * cartProd(X, U_pt) + Wk;

            Xnext = reduce(Xnext, 'girard', Kred);
            X = Xnext;
        end
    end

    contain_pct = (num_all > 0) * (100 * num_in / max(1,num_all));
    if num_all == 0, contain_pct = NaN; end
    
    sizeI_ddra = sizeI_ddra / max(1, num_all);         
    wid_k      = wid_sums ./ max(1, wid_counts);       

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

function tf = contains_interval(y, Yset, tol)
    I = interval(Yset);
    try
        lo = infimum(I); hi = supremum(I);
    catch
        lo = I.inf; hi = I.sup;  % very old CORA
    end
    tf = all(y <= hi + tol) && all(y >= lo - tol);
end

