function [U_nom, diagPE] = genPEInput(mode, L, n_u, n_k, dt, U_set, opts)
    if nargin < 7, opts = struct(); end
    Aamp = getfielddef(opts,'strength', 1);
    det  = getfielddef(opts,'deterministic', true);
    cont = getfielddef(opts,'contInput', true);
    seed_base = getfielddef(opts,'seed_base', 14621);   % new: for reproducibility

    % -------- feasibility (unchanged) --------
    L_req = L;
    L_eff = min(L_req, max(1, floor((n_k+1)/(n_u+1))));
    msg = "";
    if L_eff < L_req
        msg = sprintf('Requested L=%d is infeasible for (n_u=%d, n_k=%d); using L=%d.', ...
                      L_req, n_u, n_k, L_eff);
    end

    % -------- local RNG keyed to (L, n_u) --------
    % Ensures reproducible, L-dependent draws even if caller restores global RNG.
    % (Include dt/n_k if you want those to affect the seed as well.)
    if det
        seed = uint32( seed_base + 131*uint32(L_eff) + 977*uint32(n_u) );
        rs   = RandStream('mt19937ar','Seed',seed);
    else
        rs   = RandStream.getGlobalStream();  % use current global RNG
    end

    % -------- temporal bases S (depends on L) --------
    switch lower(string(mode))
        case "sinwave"
            t  = (0:n_k-1)*dt;
            f  = linspace(0.08/dt, 0.45/dt, L_eff);   % spread below Nyquist
            % deterministic, L-dependent phases; random only if det=false
            phi = linspace(0, pi, L_eff);
            if ~det, phi = 2*pi*rand(rs,1,L_eff); end
            S = zeros(L_eff, n_k);
            for l = 1:L_eff
                S(l,:) = sin(2*pi*f(l)*t + phi(l));
            end
        otherwise  % 'randn'
            S = randn(rs, L_eff, n_k);  % use local stream
    end

    % -------- spatial mixing V (n_u x L_eff) --------
    if isa(U_set,'zonotope') && ~isempty(U_set.G)
        GU = U_set.G;  % n_u x eta_u
        if size(GU,2) >= L_eff
            V = GU(:,1:L_eff);
        else
            V = [GU, randn(rs, size(GU,1), L_eff-size(GU,2))];  % use local stream
        end
    else
        % deterministic path still uses local stream so it's L-dependent
        V = randn(rs, n_u, L_eff);
    end
    vn = vecnorm(V); vn(vn==0) = 1; V = V ./ vn;

    % -------- build U_nom --------
    if isa(U_set,'zonotope'), cU = center(U_set); else, cU = zeros(n_u,1); end
    U_nom = cU + Aamp * (V * S);
    if cont, U_nom = cumsum(U_nom, 2); end

    % -------- diagnostics (unchanged) --------
    if exist('compute_pe_stats','file') == 2
        PE = compute_pe_stats(U_nom, n_u, L_eff);
        diagPE = struct('hankel_rank', PE.hankel_rank, ...
                        'hankel_full', PE.hankel_full, ...
                        'rank_frac',   PE.rank_frac,   ...
                        'L_req',       L_req, ...
                        'L_eff',       L_eff, ...
                        'feasible',    PE.hankel_full, ...
                        'message',     msg);
    else
        [hrank, hfull, rfrac] = local_hankel_stats(U_nom, n_u, L_eff);
        diagPE = struct('hankel_rank', hrank, ...
                        'hankel_full', hfull, ...
                        'rank_frac',   rfrac, ...
                        'L_req',       L_req, ...
                        'L_eff',       L_eff, ...
                        'feasible',    hfull, ...
                        'message',     msg);
    end

end
