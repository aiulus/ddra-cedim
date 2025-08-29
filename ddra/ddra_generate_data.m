function [Xminus, Uminus, Xplus, W, Zinfo, DATASET] = ddra_generate_data(sys, R0, U, C, pe)
%DDRA_GENERATE_DATA  Build (X^-, U^-, X^+) blocks and disturbance set W.
% Adds PE-driven nominal inputs (randn or sin) with optional randomness
% around the nominal (via U), and reports PE diagnostics in Zinfo.
%
% Inputs:
%   sys:   discrete LTI with fields A,B
%   R0,U:  zonotopes (dim_x and n_u respectively)
%   C:     config with C.shared.{n_m,n_s,n_k,p_extr}, C.ddra.{alpha_w,eta_w}
%   pe:    struct (optional)
%            .mode         : "randn" | "sinWave"    (default "randn")
%            .order        : PE order (>=1)        (default 2)
%            .strength     : amplitude scale       (default 1)
%            .deterministic: true/false            (default true)
%
% Outputs:
%   Xminus (dim_x×N), Uminus (n_u×N), Xplus (dim_x×N), W (zonotope), Zinfo

    % ------------------ unpack & basics ------------------
    n_m   = C.shared.n_m;
    n_s   = C.shared.n_s;
    n_k   = C.shared.n_k;
    dim_x = size(sys.A,1);
    n_u   = size(sys.B,2);

    % Disturbance set W
    eta_w   = getfielddef(C.ddra,'eta_w',1);
    alpha_w = getfielddef(C.ddra,'alpha_w',0);
    W = zonotope(zeros(dim_x,1), alpha_w * ones(dim_x, max(1,eta_w)));
    
    % NEW: global noise switch
    use_noise = resolve_use_noise(C.shared);
    if ~use_noise
        % hard zero: no generators, zero center
        W = zonotope(zeros(dim_x,1));  
    end

    % PE defaults
    if nargin < 5 || isempty(pe)
        pe = struct('mode','randn','order',2,'strength',1,'deterministic',true);
    end
    mode   = lower(char(getfielddef(pe,'mode','randn')));
    L      = max(1, round(getfielddef(pe,'order',2)));
    Aamp   = getfielddef(pe,'strength',1);
    det    = getfielddef(pe,'deterministic',true);
    p_extr = getfielddef(C.shared,'p_extr',0.0);

    % ------------------ nominal PE inputs (per trajectory) ------------------
    % ------------------ nominal PE inputs (per trajectory) ------------------
    U_nom_all = cell(1, n_m);
    % Honor CORA-like continuity default via options.contInput
    contInput = true;
    if isfield(C,'shared') && isfield(C.shared,'options_reach')
        % nothing needed; contInput is purely an input-generation choice
    end
    pe_opts = struct('strength', getfielddef(pe,'strength',1), ...
                     'deterministic', getfielddef(pe,'deterministic',true), ...
                     'contInput', contInput);

    for m = 1:n_m
        if getfielddef(pe,'deterministic',true)
            rng(10 + m,'twister');   % or any stable mapping of m
        end
        [U_nom, diagPE] = genPEInput(pe.mode, pe.order, n_u, n_k, sys.dt, U, pe_opts);
        if ~isempty(diagPE.message)
            fprintf('[PE] %s\n', diagPE.message);
        end
        U_nom_all{m} = U_nom;
    end

    % ------------------ assemble (X^-,U^-,X^+) ------------------
    % Preallocate exactly N columns and fill with an explicit counter.
    N = n_m * n_s * n_k;
    Xminus = zeros(dim_x, N);
    Uminus = zeros(n_u,  N);
    Xplus  = zeros(dim_x, N);

    % ------------------ NEW: DATASET ------------------
    Mtot = n_m * n_s;
    x0_list   = zeros(dim_x, Mtot);
    U_blocks  = zeros(n_u, n_k, Mtot);
    W_blocks  = zeros(dim_x, n_k, Mtot);
    X_blocks  = zeros(dim_x, n_k+1, Mtot);   % includes x0 in slot 1

    col = 0; msi = 0;
    for m = 1:n_m
        U_nom = U_nom_all{m}; % (n_u × n_k)

        for s = 1:n_s
            msi = msi + 1;

            % Perturb nominal with a sample from U (per step)
            if rand < p_extr
                U_pert = U_nom + randPoint(U, n_k, 'extreme');
            else
                U_pert = U_nom + randPoint(U, n_k);
            end

            % Simulate this block and capture everything
            x = randPoint(R0);
            x0_list(:,msi) = x;
            X_blocks(:,1,msi) = x;

            for t = 1:n_k
                col = col + 1;
                Xminus(:,col) = x;

                u = U_pert(:,t);
                Uminus(:,col) = u;

                w = randPoint(W);
                W_blocks(:,t,msi) = w;

                x = sys.A*x + sys.B*u + w;
                Xplus(:,col) = x;
                X_blocks(:,t+1,msi) = x; 
            end

            U_blocks(:,:,msi) = U_pert;
        end
    end

    % ------------------ PE diagnostics ------------------
    % Use FIRST nominal block for Hankel stats (clearer than concat).
    U1 = U_nom_all{1};                     % (n_u × n_k)
    %[hrank, hfull, rfrac] = hankel_stats(U1, n_u, L);
    Lraw = max(1, round(getfielddef(pe,'order',2)));
    Leff = min(Lraw, max(1, n_k));    % <= cap by n_k
    [hrank, hfull, rfrac] = hankel_stats(U1, n_u, Leff);

    Z = [Xminus; Uminus];
    Zinfo.rankZ      = rank(Z);
    Zinfo.condZ      = cond(Z*Z');
    Zinfo.hankel_rank= hrank;
    Zinfo.hankel_full= hfull;
    Zinfo.rank_frac  = rfrac;

    % ------------------ NEW: DATASET ------------------
    % If sys has outputs (C,D), generate y; else y = x for convenience.
    Cmat = [];
    Dmat = [];
    try
        % Struct path
        if isstruct(sys)
            if isfield(sys,'C'), Cmat = sys.C; end
            if isfield(sys,'D'), Dmat = sys.D; end
        else
            % Object path (CORA classes, ss/dss, etc.)
            Cmat = sys.C;  
            Dmat = sys.D;  
        end
    catch
        Cmat = []; Dmat = [];
    end
    
    if ~isempty(Cmat)
        ny = size(Cmat,1);
        Y_blocks = zeros(ny, n_k, Mtot);
        for i = 1:Mtot
            for t = 1:n_k
                x_t = X_blocks(:,t+1,i);
                u_t = U_blocks(:,t,i);
                if isempty(Dmat)
                    Y_blocks(:,t,i) = Cmat*x_t;
                else
                    Y_blocks(:,t,i) = Cmat*x_t + Dmat*u_t;
                end
            end
        end
    else
        ny = dim_x;
        Y_blocks = X_blocks(:,2:end,:);   
    end


    % Shapes chosen to map 1:1 into validateReach expectations:
    %  - per testcase: u (n_k x n_u), y (n_k x ny x 1), initialState (dim_x x 1)
    DATASET = struct();
    DATASET.dim_x    = dim_x;
    DATASET.n_u      = n_u;
    DATASET.n_y      = ny;
    DATASET.n_k      = n_k;
    DATASET.n_blocks = Mtot;
    DATASET.x0_list  = x0_list;                 % (dim_x x Mtot)
    DATASET.U_blocks = U_blocks;                % (n_u x n_k x Mtot)
    DATASET.W_blocks = W_blocks;                % (dim_x x n_k x Mtot)
    DATASET.X_blocks = X_blocks;                % (dim_x x (n_k+1) x Mtot)
    DATASET.Y_blocks = Y_blocks;                % (ny    x n_k     x Mtot)
end

% =================== helpers ===================
function v = getfielddef(S, f, d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end

function [hrank, hfull, rfrac] = hankel_stats(Uminus, m, L)
% Uminus: (m × N) sequence; builds block Hankel with window L.
    N = size(Uminus,2);
    if N < L
        hrank = 0; hfull = false; rfrac = 0; return;
    end
    H = block_hankel(Uminus, L);  % (mL × (N-L+1))
    hrank = rank(H);
    rfrac = hrank / min(size(H));
    hfull = (hrank >= m*L) && (hrank >= min(size(H)));
end

function H = block_hankel(U, L)
% U: (m × N) -> H: (mL × (N-L+1))
    [m, N] = size(U);
    cols = N - L + 1;
    if cols <= 0, H = zeros(m*L,0); return; end
    H = zeros(m*L, cols);
    for i = 1:L
        H((i-1)*m+1:i*m, :) = U(:, i:i+cols-1);
    end
end
