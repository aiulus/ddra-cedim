function [Xsets, sizeI] = ddra_infer(sys, R0, U, W, M_AB, C, varargin)
% DDRA_INFER  Stored-sets propagation.
% If VAL is provided and contains u, uses deterministic inputs per step.
    K  = red_order(C);           nk = C.shared.n_k;
    Xsets = cell(nk,1); Xk = R0; sizeI = 0;

    % detect deterministic VAL.u
    VAL = []; if ~isempty(varargin), VAL = varargin{1}; end
    hasVAL = isstruct(VAL) && isfield(VAL,'u') && ~isempty(VAL.u);

    if hasVAL
        % single-block convenience: use first block in VAL
        if iscell(VAL.u), Useq = VAL.u{1}; else, Useq = VAL.u; end  % (n_k×m) or (m×n_k)
        if size(Useq,2) == size(sys.B,2) && size(Useq,1) ~= nk, Useq = Useq.'; end
    end

    for k = 1:nk
        Xk = reduce(Xk,'girard',K);
        if hasVAL
            u_k = Useq(min(k,size(Useq,1)),:).';
            Ueff = zonotope(u_k);                 % deterministic per step
        else
            Ueff = U;                             % legacy: uncertain input set
        end
        Xnext = M_AB * cartProd(Xk, Ueff) + W;
        Xnext = reduce(Xnext,'girard',K);
        Xsets{k} = Xnext;
        sizeI = sizeI + sum(abs(generators(Xnext)),'all');
        Xk = Xnext;
    end
end

