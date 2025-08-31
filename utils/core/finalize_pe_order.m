function pe_eff = finalize_pe_order(pe, sys, C, varargin)
% FINALIZE_PE_ORDER  Choose an effective PE order.
% Backwards-compatible: pe_eff.order = max(1, min(max(L_min, pe.order), nk-1))
%
% Optional name-value:
%   'respect_explicit_order' (false) : if true, do NOT lift pe.order to L_min
%   'default_minimal'         (false) : if pe.order missing/empty, use L_min

    % --- parse optional flags ---
    p = inputParser;
    addParameter(p,'respect_explicit_order',false);
    addParameter(p,'default_minimal',false);
    parse(p,varargin{:});
    RESPECT = p.Results.respect_explicit_order;
    DEFMIN  = p.Results.default_minimal;

    % --- dims / lag ---
    if isfield(sys,'nrOfDims'), nx = sys.nrOfDims; else, nx = size(sys.A,1); end
    lag = getfielddef(sys,'n_p',0);
    L_min = nx + lag + 1;

    % --- requested L ---
    hasOrder = isstruct(pe) && isfield(pe,'order') && ~isempty(pe.order);
    if RESPECT && hasOrder
        L_req = pe.order;                % honor explicit order as-is
    elseif DEFMIN && (~hasOrder)
        L_req = L_min;                   % minimal sufficient if nothing given
    else
        L_req = max(L_min, getfielddef(pe,'order',L_min));   % legacy behavior
    end

    % --- clip to horizon ---
    nk  = getfielddef(C.shared,'n_k',C.shared.n_k);
    L_eff = max(1, min(L_req, max(1,nk-1)));

    % --- out ---
    pe_eff       = pe;
    pe_eff.order = L_eff;
end
