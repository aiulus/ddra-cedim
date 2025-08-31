function pe_eff = finalize_pe_order(pe, sys, C)
% Ensure persistency-of-excitation order is feasible for n_k and follows a policy.
% Backward compatible default: 'strict' -> L_min = nx + lag + 1 (your current code).

    % dims
    if isfield(sys,'nrOfDims'), nx = sys.nrOfDims; else, nx = size(sys.A,1); end
    if isfield(sys,'nrOfInputs'), nu = sys.nrOfInputs; else, nu = size(sys.B,2); end
    lag = getfielddef(sys,'n_p', 0);     % ARX lag if present
    nk  = getfielddef(C.shared,'n_k',  C.shared.n_k);

    % policy knobs (new; default keeps old behavior)
    pol = getfielddef(getfielddef(C,'shared',struct()), 'pe_min_policy', 'strict');  % 'strict' | 'rank' | 'none'
    verbose = getfielddef(getfielddef(C,'shared',struct()), 'pe_verbose', false);

    % pick floor
    switch lower(string(pol))
        case 'none'
            L_min = 1;                            % honor explicit L (PEness)
        case 'rank'
            % softer floor: enough to hope for full Hankel rank (heuristic)
            L_min = max(1, ceil((nx + lag) / max(1, nu)));
        otherwise % 'strict' (backwards-compatible)
            L_min = nx + lag + 1;
    end

    % requested vs effective
    L_req = max(1, getfielddef(pe, 'order', L_min));
    % hard feasibility cap: need at least a column in the block Hankel
    L_cap = max(1, nk - 1);
    % final
    L_eff = max( ifelse(strcmpi(pol,'none'), 1, L_min), min(L_req, L_cap) );

    if verbose
        fprintf('[finalize_pe_order] policy=%s, nx=%d, nu=%d, lag=%d, nk=%d: requested L=%d -> effective L=%d (min=%d, cap=%d)\n',...
            string(pol), nx, nu, lag, nk, L_req, L_eff, ifelse(strcmpi(pol,'none'),1,L_min), L_cap);
    end

    pe_eff = pe; pe_eff.order = L_eff;
end

function out = ifelse(cond, a, b)
    if cond, out = a; else, out = b; end
end
