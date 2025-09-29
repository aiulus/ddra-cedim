function [Xminus, Uminus, Xplus, W, Zinfo, DATASET] = ddra_generate_data(C, sys_ddra, dt, R0, U, pe, sys_sim)
%DDRA_GENERATE_DATA  Build (X^-, U^-, X^+) and disturbance set W.
% If sys_sim (CORA linearSysDT) is provided, generate trajectories via
%   simulate(sys_sim, params) exactly as CORA's createTestSuite does.
% Otherwise, fall back to manual A,B (+w) propagation (legacy).
%
% Inputs:
%   C.shared: n_m, n_s, n_k, p_extr, options_reach
%   C.ddra:   alpha_w, eta_w
%   sys_ddra: struct/object with fields A,B (for learning)
%   dt:       sample time (seconds)
%   R0,U:     zonotopes (init/input uncertainty)
%   pe:       struct('mode','...','order',L,'strength',A,'deterministic',tf,'seed_base',uint32)
%   sys_sim:  (optional) CORA linearSysDT used for simulate(...)
%
% Outputs:
%   Xminus (nx×N), Uminus (nu×N), Xplus (nx×N), W (zonotope), Zinfo, DATASET (blocks)

    % ------------------ unpack & basics ------------------
    n_m   = C.shared.n_m;
    n_s   = C.shared.n_s;
    n_k   = C.shared.n_k;
    nx    = size(sys_ddra.A,1);
    nu    = size(sys_ddra.B,2);
    useSim = (nargin >= 7) && ~isempty(sys_sim);

    if nargin < 7 || isempty(sys_sim)
        sys_sim = [];   % no simulator available -> legacy/manual path
    end

    % --- Disturbance set W: push forward input uncertainty through B  (W = B*U)
    if resolve_use_noise(C.shared)
        G_U = generators(U);                 % (nu × r)
        %G_W = sys_ddra.B * G_U;              % (nx × r)
        G_W = linearMap(sys_ddra.B, G_U);
        W   = zonotope(zeros(nx,1), G_W);    % zero center; growth only
    else
        W   = zonotope(zeros(nx,1));
    end

    % ------------------ PE nominal inputs (per trajectory m) ---------------
    U_nom_all = cell(1, n_m);
    % match createTestSuite default: continuous input
    pe_opts = struct('strength',      getfielddef(pe,'strength',1), ...
                     'deterministic', getfielddef(pe,'deterministic',true), ...
                     'contInput',     true, ...
                     'seed_base',     getfielddef(pe,'seed_base', uint32(0)));

    base_seed = getfielddef(pe,'seed_base', uint32(0));
    for m = 1:n_m
        % ensure different nominal per m (still reproducible per row)
        pe_opts_m = pe_opts;
        pe_opts_m.seed_base = uint32(mod(double(base_seed) + 44497*double(m), 2^32));
        [U_nom, diagPE]     = genPEInput(pe.mode, pe.order, nu, n_k, dt, U, pe_opts_m);
        if m == 1, diagPE_first = diagPE; end
        U_nom_all{m} = U_nom; % (nu × n_k)
    end

    % ------------------ allocate outputs ------------------
    N       = n_m * n_s * n_k;
    Xminus  = zeros(nx, N);
    Uminus  = zeros(nu, N);
    Xplus   = zeros(nx, N);

    Mtot       = n_m * n_s;
    x0_list    = zeros(nx, Mtot);
    U_blocks   = zeros(nu, n_k, Mtot);
    X_blocks   = zeros(nx, n_k+1, Mtot);   % includes x0 at index 1
    Y_blocks   = [];                        % filled only if we detect outputs later
    W_blocks   = zeros(nx, n_k, Mtot);      % kept for schema; zero if useSim==true

    % Optional output map (for DATASET.Y_blocks)
    [Cmat, Dmat, ny] = local_get_CD(sys_sim, sys_ddra, nx, nu);

    col = 0; msi = 0;
    for m = 1:n_m
        U_nom = U_nom_all{m}; % (nu × n_k)
        % one nominal initial state vector (standard uses zeros if no stateSet)
        % we keep nominal x0 = 0 and add uncertainty below (per simulation)
        x0_nom = zeros(nx,1);

        for s = 1:n_s
            msi = msi + 1;

            % ---- extreme sampling like createTestSuite -------------------
            if rand < getfielddef(C.shared,'p_extr',0)
                x0_s   = x0_nom + randPoint(R0, 1, 'extreme');
                U_pert = U_nom   + randPoint(U,  n_k, 'extreme');
            else
                x0_s   = x0_nom + randPoint(R0, 1);
                U_pert = U_nom   + randPoint(U,  n_k);
            end

            x0_list(:,msi)     = x0_s;

            % =========== PATH A: standard via simulate(sys_sim, params) ===========
            if useSim
                params = struct();
                params.R0    = zonotope(x0_s);  
                params.x0 = x0_s;
                params.u     = U_pert;            
                params.tFinal= dt*(n_k-1);
            
                try
                    [~, x_sim, u_sim, y_sim] = simulate(sys_sim, params);
                catch
                    [~, x_sim, u_sim] = simulate(sys_sim, params);
                    y_sim = [];
                end

  
                if isempty(u_sim)
                    u_sim = params.u.';           
                end
                            
                % normalize shapes to (n_k × nx/nu/ny)
                [norm_xSim, norm_uSim, norm_ySim] = local_norm_shapes(x_sim, u_sim, y_sim, n_k);

                % sanity checks
                assert(size(norm_xSim,1)==n_k, 'x_sim must be n_k×nx, got %s', mat2str(size(norm_xSim)));
                assert(size(norm_uSim,1)==n_k, 'u_sim must be n_k×nu, got %s', mat2str(size(norm_uSim)));


                % x_sim(1,:) == x0 for CORA linearSysDT. Build a clean state timeline X_tl.
                tol = 1e-10;
                X_tl = zeros(nx, n_k+1);
                X_tl(:,1) = x0_s;
                U_blocks(:,:,msi) = norm_uSim.';   % nu × n_k
                
                if norm(norm_xSim(1,:).' - x0_s) < tol
                    % norm_xSim rows are [x0; x1; ...; x_{n_k-1}] (length n_k)
                    for t = 1:n_k-1
                        X_tl(:,t+1) = norm_xSim(t+1,:).';      % x_t
                    end
                    % synthesize last state x_{n_k} with (A,B) and u(n_k)
                    X_tl(:,n_k+1) = sys_ddra.A*X_tl(:,n_k) + sys_ddra.B*norm_uSim(n_k,:).';
                else
                    % fallback: treat rows as x_t directly
                    for t = 1:n_k
                        X_tl(:,t+1) = norm_xSim(t,:).';
                    end
                end
                
                % store states (block)
                X_blocks(:,1,msi) = X_tl(:,1);
                for t = 1:n_k
                    X_blocks(:,t+1,msi) = X_tl(:,t+1);
                end
                
                % store triples
                for t = 1:n_k
                    col = col + 1;
                    Xminus(:,col) = X_tl(:,t);            % x_{t-1}
                    Uminus(:,col) = norm_uSim(t,:).';     % u_t
                    Xplus(:,col)  = X_tl(:,t+1);          % x_t
                end

                % outputs (keeps your previous semantics that y(1)=C x0 + D u1)
                if ny > 0
                    if isempty(Y_blocks), Y_blocks = zeros(ny, n_k, Mtot); end
                    if ~isempty(norm_ySim)
                        Y_blocks(:,:,msi) = norm_ySim.';          % ny × n_k
                    else
                        for t = 1:n_k
                            Y_blocks(:,t,msi) = Cmat*X_tl(:,t) + Dmat*U_blocks(:,t,msi);
                        end
                    end
                end

                % W_blocks stay zero (no process noise injection in standard)

            % =========== PATH B: legacy manual propagation ===========
            else
                X_blocks(:,1,msi) = x0_s;
                x = x0_s;
                for t = 1:n_k
                    col = col + 1;
                    Xminus(:,col) = x;
                    u             = U_pert(:,t);
                    Uminus(:,col) = u;
                    U_blocks(:,t,msi) = u;     % NEW
                
                    w  = randPoint(W);
                    W_blocks(:,t,msi) = w;
                
                    x = sys_ddra.A*x + sys_ddra.B*u + w;
                    Xplus(:,col) = x;
                    X_blocks(:,t+1,msi) = x;
                end


                if ny > 0
                    if isempty(Y_blocks), Y_blocks = zeros(ny, n_k, Mtot); end
                    for t = 1:n_k
                        Y_blocks(:,t,msi) = Cmat*X_blocks(:,t+1,msi) + Dmat*U_blocks(:,t,msi);
                    end
                end
            end
        end
    end

    % ------------------ PE / regressor diagnostics ------------------
    % use first nominal for Hankel stats (clearer than concat)
    U1 = U_nom_all{1};  % (nu × n_k)
    Lraw = max(1, round(getfielddef(pe,'order',2)));
    Leff = min(Lraw, max(1, n_k));
    [hrank, hfull, rfrac] = hankel_stats(U1, nu, Leff);

    Z = [Xminus; Uminus];
    Zinfo.rankZ        = rank(Z);
    Zinfo.condZ        = cond(Z*Z');
    Zinfo.pe_L_req     = Lraw;
    Zinfo.pe_L_eff     = getfielddef(diagPE_first,'L_eff', Leff);
    Zinfo.hankel_rank  = getfielddef(diagPE_first,'hankel_rank', hrank);
    Zinfo.hankel_full  = getfielddef(diagPE_first,'hankel_full', hfull);
    Zinfo.rank_frac    = getfielddef(diagPE_first,'rank_frac', rfrac);
    Zinfo.pe_message   = string(getfielddef(diagPE_first,'message', ""));

    % ------------------ DATASET (matches validateReach expectations) ------
    DATASET = struct();
    DATASET.dim_x    = nx;
    DATASET.n_u      = nu;
    DATASET.n_y      = max(0, ny);
    DATASET.n_k      = n_k;
    DATASET.n_blocks = Mtot;
    DATASET.x0_list  = x0_list;                   % (nx × Mtot)
    DATASET.U_blocks = U_blocks;                  % (nu × n_k × Mtot)
    DATASET.W_blocks = W_blocks;                  % (nx × n_k × Mtot)
    DATASET.X_blocks = X_blocks;                  % (nx × (n_k+1) × Mtot)
    DATASET.U_nom_all = U_nom_all;                
    if ny > 0
        DATASET.Y_blocks = Y_blocks;              % (ny × n_k × Mtot)
    end
end

% =================== helpers ===================
function v = getfielddef(S, f, d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end

function use_noise = resolve_use_noise(S)
    if isfield(S,'noise_for_ddra') || isfield(S,'noise_for_gray')
        g = ~isfield(S,'noise_for_gray') || logical(S.noise_for_gray);
        d = ~isfield(S,'noise_for_ddra') || logical(S.noise_for_ddra);
        use_noise = g && d;
    else
        use_noise = true;
    end
end

function [Cmat,Dmat,ny] = local_get_CD(sys_sim, sys_ddra, nx, nu)
    Cmat = []; Dmat = [];
    if nargin>=1 && ~isempty(sys_sim)
        try, Cmat = sys_sim.C; catch, end
        try, Dmat = sys_sim.D; catch, end
    end
    if isempty(Cmat)
        try, Cmat = sys_ddra.C; catch, end
    end
    if isempty(Cmat), Cmat = eye(nx); end
    ny = size(Cmat,1);

    if isempty(Dmat)
        try, Dmat = sys_ddra.D; catch, end
        if isempty(Dmat), Dmat = zeros(ny,nu); end
    end
end

function [xS, uS, yS] = local_norm_shapes(x_sim, u_sim, y_sim, n_k)
    xS = []; uS = []; yS = [];

    if ~isempty(x_sim)
        xs = squeeze(x_sim);
        % want n_k × n_x
        if size(xs,2) == n_k, xs = xs.'; end
        xS = xs;
    end

    if ~isempty(u_sim)
        us = squeeze(u_sim);
        % simulate often returns m × (steps+1) = m × n_k
        if size(us,2) == n_k
            us = us.';           % -> n_k × m
        elseif size(us,1) ~= n_k
            % last resort: if neither dim equals n_k, raise a clear error
            error('Unexpected u_sim size %s; cannot normalize to %d×m.', mat2str(size(us)), n_k);
        end
        uS = us;
    end

    if ~isempty(y_sim)
        ys = squeeze(y_sim);
        % want n_k × n_y
        if size(ys,2) == n_k, ys = ys.'; end
        yS = ys;
    end
end

function [hrank, hfull, rfrac] = hankel_stats(Uminus, m, L)
    N = size(Uminus,2);
    if N < L, hrank=0; hfull=false; rfrac=0; return; end
    H = block_hankel(Uminus, L);
    hrank = rank(H);
    rfrac = hrank / min(size(H));
    hfull = (hrank >= m*L) && (hrank >= min(size(H)));
end

function H = block_hankel(U, L)
    [m, N] = size(U);
    cols = N - L + 1; if cols <= 0, H = zeros(m*L,0); return; end
    H = zeros(m*L, cols);
    for i = 1:L
        H((i-1)*m+1:i*m, :) = U(:, i:i+cols-1);
    end
end
