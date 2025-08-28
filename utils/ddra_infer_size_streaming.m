function [sizeI_ddra, contain_pct] = ddra_infer_size_streaming(sys, R0, U, W, M_AB, C, VAL)
% Streaming DDRA: unified size + containment on the SAME validation points as Gray.
% VAL must be built by build_validation_bundle().

    assert(nargin>=7 && ~isempty(VAL) && isfield(VAL,'y') && isfield(VAL,'u'), ...
        'ddra_infer_size_streaming: VAL (y,u,...) is required for synchronized validation.');

    K   = C.shared.options_reach.zonotopeOrder;   % unified reduction
    X   = R0;
    S   = 0;   % size accumulator
    hit = 0; tot = 0;

    n_k_val = C.shared.n_k_val;
    if isfield(VAL,'n_k'), n_k_val = min(n_k_val, VAL.n_k); end

    % pull output maps; missing -> identity/zero
    dim_y = size(VAL.y,2);
    hasC = isfield(VAL,'C') && ~isempty(VAL.C);
    hasD = isfield(VAL,'D') && ~isempty(VAL.D);
    if ~hasC, VAL.C = eye(size(sys.A,1), size(sys.A,1)); end
    if ~hasD, VAL.D = zeros(dim_y, size(sys.B,2)); end

    for k = 1:n_k_val
        % state set step (pre-reduce -> image -> post-reduce)
        X = reduce(X,'girard',K);
        Xnext = M_AB * cartProd(X, U) + W;
        Xnext = reduce(Xnext,'girard',K);

        % unified size metric in OUTPUT space: size of Y = C*Xnext
        Yset_base = VAL.C * Xnext;       % matZonotope/zonotope linear map
        S = S + size_interval_sum(Yset_base);

        % unified containment: y_k^(s) ∈ C Xnext + D u_k^(s)
        yk = squeeze(VAL.y(k,:,:));      % (dim_y × Sval)
        uk = squeeze(VAL.u(k,:,:));      % (dim_u × Sval)
        if isvector(yk), yk = yk(:); end
        if isvector(uk), uk = uk(:); end
        Sval = size(yk,2);

        for s = 1:Sval
            ys = yk(:,s);
            dus = VAL.D * uk(:,s);
            % Shift measurement into the frame of Yset_base:
            % ys ∈ Yset_base + dus   <=>   ys - dus ∈ Yset_base
            inside = contains_interval(ys - dus, Yset_base, 1e-6);
            hit = hit + inside;
            tot = tot + 1;
        end

        X = Xnext;
    end

    sizeI_ddra  = S;
    contain_pct = 100 * hit / max(1,tot);
end
