function [Xminus, Uminus, Xplus, W, Zinfo] = ddra_generate_data(sys, R0, U, C, pe)
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
    eta_w = getfielddef(C.ddra,'eta_w',1);
    alpha_w = getfielddef(C.ddra,'alpha_w',0);
    W = zonotope(zeros(dim_x,1), alpha_w * ones(dim_x, max(1,eta_w)));

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
    % Build U_nom_all{m} as (n_u × n_k)
    U_nom_all = cell(1, n_m);
    Leff_time = max(1, min(L, max(1, n_k-1)));   % ensure at least 1, at most n_k-1
    Leff      = Leff_time;

    for m = 1:n_m
        if det, rng(1000+m,'twister'); end

        % Temporal basis S (Leff × n_k)
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
            GU = U.generators;                          % (n_u × eta_u)
            if size(GU,2) >= Leff
                V = GU(:,1:Leff);
            else
                V = [GU, randn(size(GU,1), Leff-size(GU,2))];
            end
        else
            V = randn(n_u, Leff);
        end
        % Normalize columns (avoid huge scales if GU small)
        vn = vecnorm(V); vn(vn==0) = 1;
        V = V ./ vn;

        % Center of U (if available)
        if isa(U,'zonotope'), cU = center(U); else, cU = zeros(n_u,1); end

        U_nom = cU + Aamp * (V * S);

        % Optional: accumulate for smoother, "continuous-ish" inputs
        if n_k >= 2
            U_nom = cumsum(U_nom, 2);
        end

        U_nom_all{m} = U_nom;
    end

    % ------------------ assemble (X^-,U^-,X^+) ------------------
    % Preallocate exactly N columns and fill with an explicit counter.
    N = n_m * n_s * n_k;
    Xminus = zeros(dim_x, N);
    Uminus = zeros(n_u,  N);
    Xplus  = zeros(dim_x, N);

    col = 0;
    for m = 1<n_m
        U_nom = U_nom_all{m}; % (n_u × n_k), reuse across samples s

        for s = 1:n_s
            % Perturb nominal with a sample from U (per step)
            if rand < p_extr
                U_pert = U_nom + randPoint(U, n_k, 'extreme'); % (n_u × n_k)
            else
                U_pert = U_nom + randPoint(U, n_k);            % (n_u × n_k)
            end

            % Simulate this block
            x = randPoint(R0);
            for t = 1:n_k
                col = col + 1;
                Xminus(:,col) = x;

                u = U_pert(:,t);
                Uminus(:,col) = u;

                w = randPoint(W);
                x = sys.A*x + sys.B*u + w;

                Xplus(:,col) = x;
            end
        end
    end

    % ------------------ PE diagnostics ------------------
    % Use FIRST nominal block for Hankel stats (clearer than concat).
    U1 = U_nom_all{1};                     % (n_u × n_k)
    [hrank, hfull, rfrac] = hankel_stats(U1, n_u, L);

    Z = [Xminus; Uminus];
    Zinfo.rankZ      = rank(Z);
    Zinfo.condZ      = cond(Z*Z');
    Zinfo.hankel_rank= hrank;
    Zinfo.hankel_full= hfull;
    Zinfo.rank_frac  = rfrac;
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
