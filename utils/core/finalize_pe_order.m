function pe_eff = finalize_pe_order(pe, sys, C)
% Ensure persistency-of-excitation is at least (dim + lag + 1), but feasible for n_k
% Works for linearSysDT / linearARX; falls back safely for others.
    if isfield(sys,'nrOfDims'),   nx = sys.nrOfDims; else, nx = size(sys.A,1); end
    lag = getfielddef(sys,'n_p', 0);  % ARX lag if present, else 0
    L_min = nx + lag + 1;

    nk    = getfielddef(C.shared,'n_k',  C.shared.n_k);
    L_req = max(L_min, getfielddef(pe,'order', L_min));

    % keep at least 1, at most nk-1 (need space to realize L)
    L_eff = max(1, min(L_req, max(1, nk-1)));

    pe_eff = pe; pe_eff.order = L_eff;
end
