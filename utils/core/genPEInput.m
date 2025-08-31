function [U_nom, diagPE] = genPEInput(mode, L, n_u, n_k, dt, U_set, opts)
%GENPEINPUT  Build a nominal input U_nom (n_u x n_k) with target PE order L.
%   mode : 'randn' | 'sinWave'
%   L    : desired PE order (>=1)
%   n_u  : #inputs
%   n_k  : horizon length (timesteps)
%   dt   : sample time
%   U_set: zonotope for inputs (used to set directions/scale); can be []
%   opts : struct with fields (all optional)
%            .strength       (default 1)
%            .deterministic  (default true -> reproducible phases/noise)
%            .contInput      (default true -> cumsum for continuity)
%
% Returns:
%   U_nom : (n_u x n_k) nominal input
%   diagPE: struct with Hankel diagnostics (rank, full, fraction, message)

    if nargin < 7, opts = struct(); end
    Aamp = getfielddef(opts,'strength', 1);
    det  = getfielddef(opts,'deterministic', true);
    cont = getfielddef(opts,'contInput', true);

    % ---- feasibility & effective window length ---------------------------
    % Block Hankel H has size (n_u*L) x (n_k-L+1). To be able to reach
    % full row rank n_u*L, we must have n_k-L+1 >= n_u*L  
    % n_k >= (n_u+1)*L - 1. If not, we reduce L and flag it.
    L_req = L;
    L_eff = min(L_req, max(1, floor((n_k+1)/(n_u+1))));
    msg = "";
    if L_eff < L_req
        msg = sprintf('Requested L=%d is infeasible for (n_u=%d, n_k=%d); using L=%d.', ...
                      L_req, n_u, n_k, L_eff);
    end

    % ---- temporal bases S (L_eff x n_k) ---------------------------------
    switch lower(string(mode))
        case "sinwave"
            t  = (0:n_k-1)*dt;
            % Spread L_eff distinct frequencies over a band; fixed phases in det mode
            f  = linspace(0.08/dt, 0.45/dt, L_eff);      % [Hz] equiv (per-sample band)
            if det
                phi = linspace(0, pi, L_eff);
            else
                phi = 2*pi*rand(1, L_eff);
            end
            S = zeros(L_eff, n_k);
            for l = 1:L_eff
                S(l,:) = sin(2*pi*f(l)*t + phi(l));
            end
        otherwise % 'randn'
            S = randn(L_eff, n_k);
    end

    % ---- spatial mixing V (n_u x L_eff) ---------------------------------
    if isa(U_set,'zonotope') && ~isempty(U_set.G)
        GU = U_set.G;  % n_u x eta_u
        if size(GU,2) >= L_eff
            V = GU(:,1:L_eff);
        else
            V = [GU, randn(size(GU,1), L_eff-size(GU,2))];
        end
    else
        V = randn(n_u, L_eff);
    end
    vn = vecnorm(V); vn(vn==0) = 1; V = V ./ vn;

    % ---- center + scale --------------------------------------------------
    if isa(U_set,'zonotope'), cU = center(U_set); else, cU = zeros(n_u,1); end
    U_nom = cU + Aamp * (V * S);

    % Optional continuity (matches createTestSuite’s default contInput=true)
    if cont, U_nom = cumsum(U_nom, 2); end

    % ---- PE diagnostics on the built sequence ----------------------------
    [hrank, hfull, rfrac] = hankel_stats(U_nom, n_u, L_eff);
    diagPE = struct('hankel_rank', hrank, ...
                    'hankel_full', hfull, ...
                    'rank_frac',   rfrac, ...
                    'L_req',       L_req, ...
                    'L_eff',       L_eff, ...
                    'feasible',    hfull, ...
                    'message',     msg);
end

% ------- helpers (same signatures as in your file) -----------------------
function [hrank, hfull, rfrac] = hankel_stats(U, m, L)
    N = size(U,2);
    if N < L, hrank=0; hfull=false; rfrac=0; return; end
    H = block_hankel(U, L);  % (mL x (N-L+1))
    hrank = rank(H);
    rfrac = hrank / min(size(H));
    hfull = (hrank >= m*L) && (hrank >= min(size(H)));
end

function H = block_hankel(U, L)
    [m, N] = size(U);
    cols = N - L + 1;
    if cols <= 0, H = zeros(m*L,0); return; end
    H = zeros(m*L, cols);
    for i = 1:L
        H((i-1)*m+1:i*m, :) = U(:, i:i+cols-1);
    end
end

function v = getfielddef(S, f, d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end
