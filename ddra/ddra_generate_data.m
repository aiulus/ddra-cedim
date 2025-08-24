function [Xminus, Uminus, Xplus, W, Zinfo] = ddra_generate_data(sys, R0, U, C, pe)
% Simulate blocks to build (X^-, U^-, X^+) and define W.
% Adds PE-driven nominal inputs (randn or sin) with optional randomness
% around the nominal (via U), and reports PE diagnostics in Zinfo.
%
% Inputs:
%   sys:   discrete LTI with fields A,B (as before)
%   R0,U:  zonotopes
%   C:     config with C.shared.{n_m,n_s,n_k,p_extr}, C.ddra.{alpha_w,eta_w}
%   pe:    struct with fields (all optional)
%            .mode  : "randn" | "sinWave"
%            .order : desired PE order (integer >=1)
%            .strength : amplitude scale (default 1)
%            .deterministic : true/false (default true)
%

    % ------------------ unpack ------------------
    n_m = C.shared.n_m; 
    n_s = C.shared.n_s; 
    n_k = C.shared.n_k; 
    dim_x = size(sys.A,1);
    n_u = size(sys.B,2);

    % Process noise W (unchanged)
    eta_w = getfielddef(C.ddra,'eta_w',1);
    W = zonotope(zeros(dim_x,1), C.ddra.alpha_w * ones(dim_x, eta_w));

    % ------------------ PE nominal inputs ------------------
    % Defaults (preserve old behavior if pe not provided)
    if nargin < 5 || isempty(pe)
        pe = struct('mode','randn','order',2,'strength',1,'deterministic',true);
    end
    mode  = lower(char(getfielddef(pe,'mode','randn')));
    L     = max(1, round(getfielddef(pe,'order',2)));
    Aamp  = getfielddef(pe,'strength',1);
    det   = getfielddef(pe,'deterministic',true);
    p_extr= getfielddef(C.shared,'p_extr',0.0);

    % Build n_m nominal input matrices U_nom_all{m} (n_u x n_k)
    U_nom_all = cell(1, n_m);
    for m = 1:n_m
        if det, rng(1000+m,'twister'); end
        % --- Build a Leff-rank temporal basis and map to input space ---
        Leff_time = max(1, min(L, n_k-1));
        Leff      = Leff_time;
        
        % Temporal basis S (Leff x n_k)
        switch mode
            case 'sinwave'
                t = 0:(n_k-1);
                f = linspace(0.08, 0.45, Leff);      % spread frequencies
                phi = 2*pi*rand(1,Leff);
                S = zeros(Leff, n_k);
                for l = 1:Leff
                    S(l,:) = sin(2*pi*f(l)*t + phi(l));
                end
            otherwise  % 'randn'
                S = randn(Leff, n_k);                
        end
        
        % Map temporal bases to input space with Leff directions
        if isa(U,'zonotope') && ~isempty(U.generators)
            % use first Leff generator directions; fall back to random if fewer
            GU = U.generators;                       % dim_u x eta_u
            if size(GU,2) >= Leff
                V = GU(:,1:Leff);
            else
                V = [GU, randn(size(GU,1), Leff-size(GU,2))];
            end
        else
            % no explicit input generators: synthesize orthonormal-ish directions
            V = randn(n_u, Leff); V = V ./ max(1e-12, vecnorm(V));
        end
        
        if isa(U,'zonotope'), cU = center(U); else, cU = zeros(n_u,1); end
        U_nom = cU + Aamp * (V * S);
        
        % Optional: match CORA’s “continuous input” feel
        U_nom = cumsum(U_nom, 2);

        U_nom_all{m} = U_nom;
    end

    % ------------------ assemble blocks ------------------
    % For each (m,s), perturb nominal with a sample from U (like createTestSuite)
    %n_blocks = n_m * n_s;
    %Xminus = []; Uminus = []; Xplus = [];
    N = n_m * n_s * n_k;
    Xminus = zeros(dim_x, N);
    Uminus = zeros(n_u,  N);
    Xplus  = zeros(dim_x, N);

    for m = 1:n_m
        for s = 1:n_s
            U_nom = U_nom_all{m}; % (n_u x n_k)
            if rand < p_extr
                U_pert = U_nom + randPoint(U, n_k, 'extreme'); % (n_u x n_k)
            else
                U_pert = U_nom + randPoint(U, n_k);            % (n_u x n_k)
            end

            % Simulate this block
            x = randPoint(R0);
            for t = 1:n_k
                Xminus(:,end+1) = x; 
                u = U_pert(:,t); 
                Uminus(:,end+1) = u;
                w = randPoint(W);
                x = sys.A*x + sys.B*u + w;
                Xplus(:,end+1) = x;
            end
        end
    end

    % ------------------ rank/conditioning + PE diagnostics ------------------
    % --- PE diagnostics on the FIRST block only (clearer than on concatenation)
    U1 = U_nom_all{1};                % dim_u x n_k (or use the perturbed one if you prefer)
    [hrank, hfull, rfrac] = hankel_stats(U1, n_u, L);

    Z = [Xminus; Uminus];
    Zinfo.rankZ = rank(Z);
    Zinfo.condZ = cond(Z*Z');

    Zinfo.hankel_rank = hrank;
    Zinfo.hankel_full = hfull;
    Zinfo.rank_frac   = rfrac;
end

% =================== helpers ===================
function v = getfielddef(S, f, d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end

function [hrank, hfull, rfrac] = hankel_stats(Uminus, m, L)
% Uminus: (m x N) concatenated inputs
    N = size(Uminus,2);
    if N < L
        hrank = 0; hfull = false; rfrac = 0; return;
    end
    H = block_hankel(Uminus, L);  % (mL x (N-L+1))
    hrank = rank(H);
    rfrac = hrank / min(size(H));
    hfull = (hrank >= m*L) && (hrank >= min(size(H)));
end

function H = block_hankel(U, L)
% U: (m x N) -> H: (mL x (N-L+1))
    [m, N] = size(U);
    cols = N - L + 1;
    if cols <= 0
        H = zeros(m*L,0); return;
    end
    H = zeros(m*L, cols);
    for i = 1:L
        H((i-1)*m+1:i*m, :) = U(:, i:i+cols-1);
    end
end
