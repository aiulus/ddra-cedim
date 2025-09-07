function [U_nom, diagPE] = genPEInput(mode, L, n_u, n_k, dt, U_set, opts)
    if nargin < 7, opts = struct(); end
    Aamp = getfielddef(opts,'strength', 1);
    det  = getfielddef(opts,'deterministic', true);
    cont = getfielddef(opts,'contInput', true);
    seed_base = getfielddef(opts,'seed_base', 14621);

    % --- feasibility ---
    L_req = L;
    L_eff = min(L_req, max(1, floor((n_k+1)/(n_u+1))));
    msg = "";
    if L_eff < L_req
        msg = sprintf('Requested L=%d is infeasible for (n_u=%d, n_k=%d); using L=%d.', ...
                      L_req, n_u, n_k, L_eff);
    end

    % --- local RNG ---
    if det
        seed = uint32( seed_base + 131*uint32(L_eff) + 977*uint32(n_u) );
        rs   = RandStream('mt19937ar','Seed',seed);
    else
        rs   = RandStream.getGlobalStream();
    end

    built_direct = false;   % true for *CS modes that build U_nom directly*

    switch lower(string(mode))
        case "sinwave"
            t  = (0:n_k-1)*dt;
            f  = linspace(0.08/dt, 0.45/dt, L_eff);
            phi = linspace(0, pi, L_eff);
            if ~det, phi = 2*pi*rand(rs,1,L_eff); end
            S = zeros(L_eff, n_k);
            for l = 1:L_eff
                S(l,:) = sin(2*pi*f(l)*t + phi(l));
            end

        case "randn"
            S = randn(rs, L_eff, n_k);

        case "randncs"   % emulate createTestSuite "randn"
            % default parameter ranges from createTestSuite
            d_min = 0; d_max = (n_k*dt)/4;
            delay = d_min + (d_max-d_min)*rand(rs,n_u,1);
            A     = 2*rand(rs,n_u,1)-1;           % [-1,1]
            tZero = (3/4)*(n_k*dt) + (n_k*dt - (3/4)*(n_k*dt))*rand(rs,n_u,1);

            U_nom = zeros(n_u,n_k);
            for i=1:n_u
                t1 = max(ceil(delay(i)/dt),1);
                t2 = min(floor(tZero(i)/dt), n_k);
                if t2 >= t1
                    len = t2-t1+1;
                    U_nom(i,t1:t2) = A(i) * randn(rs, 1, len);
                end
            end
            if cont, U_nom = cumsum(U_nom,2); end
            U_nom = Aamp * U_nom;           % apply global strength
            built_direct = true;

        case "sinwavecs" % emulate createTestSuite "sinWave"
            d_min = 0; d_max = (n_k*dt)/4;
            delay = d_min + (d_max-d_min)*rand(rs,n_u,1);
            A     = 2*rand(rs,n_u,1)-1;
            % crude but OK proxy for T distribution used in example
            T     = max(dt, (n_k*dt)/10*rand(rs,n_u,1));  T = 10*T;
            Ts    = zeros(n_u,1);

            U_nom = zeros(n_u,n_k);
            tgrid = (1:n_k)*dt;
            for i=1:n_u
                idx1 = (tgrid >= delay(i)) & (tgrid <= min((T(i)/4 + delay(i)), n_k*dt));
                idx2 = (tgrid >= (T(i)/4+Ts(i)+delay(i))) & ...
                       (tgrid <= min((T(i)/2+Ts(i)+delay(i)), n_k*dt));
                U_nom(i,idx1) = A(i)*cos((2*pi/T(i))*(tgrid(idx1)-delay(i)));
                U_nom(i,idx2) = A(i)*cos((2*pi/T(i))*(tgrid(idx2)-Ts(i)-delay(i)));
            end
            if cont, U_nom = cumsum(U_nom,2); end
            U_nom = Aamp * U_nom;
            built_direct = true;

        otherwise
            % fall back to random bases
            S = randn(rs, L_eff, n_k);
    end

    if ~built_direct
        % spatial mixing using U_set's directions (optional)
        if isa(U_set,'zonotope') && ~isempty(U_set.G)
            GU = U_set.G;  % n_u × eta_u
            if size(GU,2) >= L_eff
                V = GU(:,1:L_eff);
            else
                V = [GU, randn(rs, size(GU,1), L_eff-size(GU,2))];
            end
        else
            V = randn(rs, n_u, L_eff);
        end
        vn = vecnorm(V); vn(vn==0) = 1; V = V ./ vn;

        % build U_nom (no cU to mirror createTestSuite nominal)
        U_nom = Aamp * (V * S);
        if cont, U_nom = cumsum(U_nom,2); end
    end

    % --- diagnostics ---
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
        [hrank, hfull, rfrac] = hankel_stats(U_nom, n_u, L_eff);  
        diagPE = struct('hankel_rank', hrank, ...
                        'hankel_full', hfull, ...
                        'rank_frac',   rfrac, ...
                        'L_req',       L_req, ...
                        'L_eff',       L_eff, ...
                        'feasible',    hfull, ...
                        'message',     msg);
    end
end
