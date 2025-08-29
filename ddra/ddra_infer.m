function [Xsets, sizeI] = ddra_infer(sys, R0, U, W, M_AB, C, varargin) 
    K = red_order(C);
    nk = C.shared.n_k;

    Xsets = cell(nk,1);
    Xk = R0;
    sizeI = 0;

    for k = 1:nk
        Xk    = reduce(Xk,'girard',K);
        Xnext = M_AB * cartProd(Xk, U) + W;   
        Xnext = reduce(Xnext,'girard',K);

        Xsets{k} = Xnext;
        sizeI    = sizeI + sum(abs(generators(Xnext)),'all');
        Xk       = Xnext;
    end
end


% ------------------------ helpers (local) --------------------------------
function u_k = val_get_u_k(VAL, k, nk, dim_u)
    % Try the most explicit schema first: u_blocks{b} is [n_k x dim_u]
    if isfield(VAL,'u_blocks') && ~isempty(VAL.u_blocks) && ~isempty(VAL.u_blocks{1})
        Ublk = VAL.u_blocks{1};
        % allow longer blocks; clip if needed
        k_eff = min(k, size(Ublk,1));
        u_k = Ublk(k_eff,:).';  % column
        return
    end

    % Fallback: VAL.u either [n_k x dim_u] or [dim_u x n_k] or 3D
    if isfield(VAL,'u') && ~isempty(VAL.u)
        Uraw = VAL.u;
        sz = size(Uraw);
        if numel(sz) == 2
            if sz(1) == nk && sz(2) == dim_u
                u_k = Uraw(min(k,sz(1)),:).';
                return
            elseif sz(2) == nk && sz(1) == dim_u
                u_k = Uraw(:,min(k,sz(2)));
                return
            end
        elseif numel(sz) == 3
            % Common shapes: [n_k x dim_u x n_cases] or [dim_u x n_k x n_cases]
            if sz(1) == nk && sz(2) == dim_u
                u_k = Uraw(min(k,sz(1)),:,1).';
                return
            elseif sz(2) == nk && sz(1) == dim_u
                u_k = Uraw(:,min(k,sz(2)),1);
                return
            end
        end
    end

    % Last resort: use the *center* of U (keeps old behavior)
    warning('ddra_infer:VALshape',...
        'VAL provided but no recognizable u-sequence; falling back to center(U).');
    u_k = center(U);
end
