function pe_eff = pe_normalize(pe, U, sys_or_nx, n_k, varargin)
%PE_NORMALIZE  Canonicalize PE spec across callers.
% pe: struct with fields you already use: mode, order (optional), strength or strength_rel, deterministic, seed_base
% U : zonotope for input uncertainty (to scale relative amplitudes)
% sys_or_nx: CORA sys (with .nrOfDims) OR numeric nx
% n_k: horizon length
% opts:
%   'respect_explicit_order' (default true)
%   'default_minimal'        (default true)  % if order missing, use L_min = nx+lag+1
%
% Returns a struct pe_eff with concrete fields:
%   .mode, .order (clamped to [1, n_k-1]), .strength (scalar), .deterministic, .seed_base

    p = inputParser;
    addParameter(p,'respect_explicit_order',true);
    addParameter(p,'default_minimal',true);
    parse(p,varargin{:});
    RESPECT = p.Results.respect_explicit_order;
    DEFMIN  = p.Results.default_minimal;

    % --- dims / minimal L ---
    if isnumeric(sys_or_nx), nx = sys_or_nx; else, nx = sys_or_nx.nrOfDims; end
    lag = 0;                               % ARX lag if relevant; 0 for your MSD path
    Lmin = nx + lag + 1;

    % --- order ---
    hasOrder = isfield(pe,'order') && ~isempty(pe.order);
    if RESPECT && hasOrder
        Lreq = pe.order;
    elseif ~hasOrder && DEFMIN
        Lreq = Lmin;
    else
        Lreq = getfield_default(pe,'order', Lmin);
    end
    Leff = max(1, min(Lreq, max(1,n_k-1)));

    % --- strength (absolute from relative if requested) ---
    s = getfield_default(pe,'strength', []);
    if isempty(s)
        s_rel = getfield_default(pe,'strength_rel', 0.8);   % default: 80% of U half-width
        try
            G = generators(U);                   % nu × r
            if isempty(G), hw = ones(dim(U),1);
            else, hw = sum(abs(G),2);           % per-channel half-width
            end
            s = mean(s_rel .* hw);              % collapse to scalar for genPEInput
        catch
            s = 1; % safe fallback
        end
    end

    % --- pack ---
    pe_eff = struct();
    pe_eff.mode          = pe.mode;
    pe_eff.order         = Leff;
    pe_eff.strength      = s;
    pe_eff.deterministic = getfield_default(pe,'deterministic',true);
    pe_eff.seed_base     = getfield_default(pe,'seed_base', uint32(0));

    % Optional: show normalization
    if getfield_default(pe,'verbose',false)
        fprintf('[PE] mode=%s, L_req=%g -> L_eff=%g (n_k=%g), strength=%g\n', ...
            string(pe_eff.mode), Lreq, pe_eff.order, n_k, pe_eff.strength);
    end
end

function v = getfield_default(S,f,d)
    if isstruct(S) && isfield(S,f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end
