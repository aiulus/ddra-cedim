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
    W = zonotope(zeros(dim_x,1), C.ddra.alpha_w*ones(dim_x, C.ddra.eta_w));

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
        switch mode
            case 'randn'
                % One active random segment long enough to help PE of order L
                delay = 0;
                tZero = min(n_k, max(ceil(0.3*n_k), ceil(0.5*n_k + L)));
                U_nom = zeros(n_u, n_k);
                idx1 = max(1, delay+1);
                idx2 = min(n_k, tZero);
                if idx2 >= idx1
                    U_nom(:, idx1:idx2) = Aamp * randn(n_u, idx2-idx1+1);
                end
                % Match CORA default contInput = true
                U_nom = cumsum(U_nom, 2);

            case 'sinwave'
                % Period (in samples) ~ n_k / max(1,L/2); ensure >=8 for shape
                T = max(8, round(n_k / max(1, L/2)));
                Ts = 0; delay = 0;
                t = 0:(n_k-1);
                U_nom = Aamp * cos((2*pi/T) * (t - delay));
                U_nom = repmat(U_nom, n_u, 1); % same phase per channel for now
                % Match CORA default contInput = true
                U_nom = cumsum(U_nom, 2);

            otherwise
                % Fallback: pure random-in-U 
                U_nom = zeros(n_u, n_k);
        end
        U_nom_all{m} = U_nom;
    end

    % ------------------ assemble blocks ------------------
    % For each (m,s), perturb nominal with a sample from U (like createTestSuite)
    n_blocks = n_m * n_s;
    Xminus = []; Uminus = []; Xplus = [];

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
    Z = [Xminus; Uminus];
    Zinfo.rankZ = rank(Z);
    Zinfo.condZ = cond(Z*Z');

    % PE stats based on block-Hankel of the whole U sequence
    [h_rank, h_full, rank_frac] = hankel_stats(Uminus, n_u, L);
    Zinfo.hankel_rank = h_rank;
    Zinfo.hankel_full = h_full;
    Zinfo.rank_frac   = rank_frac; % r / min(size(H))

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
