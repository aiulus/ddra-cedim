function M_AB = ddra_linear_learn(sys, lookup)
    %% Expects:
    %   * sys = c2d(.)
    %   * lookup = struct(...) 
    %
    % Data layout:
    %   X_minus, U_minus, X_plus : dim_x × (n_m*n_s*n_k), ordered as
    %       (m=1,s=1,k=0..n_k-1), (m=1,s=2,k=0..n_k-1), ...,
    %       (m=2,s=1,k=0..n_k-1), ..., (m=n_m,s=n_s,k=0..n_k-1).
    %
    % Original Author: Amr Alanwar (28-Oct-2020)
    % Refactor: Aybüke Ulusarslan (21-Aug-2025)
    % Notes: generalized to n_s >= 1 samples per test case; n_s=1 reproduces original

    %------------- BEGIN CODE --------------

    rng(1, 'twister');  

    %% system dimensions
    dim_x = lookup.dim_x; 
    dim_u = lookup.dim_u; 
    dim_y = lookup.dim_y; 

    n_s = lookup.n_s;   % samples per test case
    n_m = lookup.n_m;   % number of test cases
    n_k = lookup.n_k;   % trajectory length per sample
    totalsamples = n_m * n_s * n_k;  % total time-samples across all (m,s)

    %% initial-state uncertainty X0
    eta_x = lookup.eta_x; 
    alpha_x = lookup.alpha_x; 
    G_x0 = diag(ones(dim_x, eta_x));
    c_x = lookup.c_x;
    c_delta_x = lookup.c_delta_x;
    c_x = c_x + c_delta_x;
    X0 = zonotope(c_x, alpha_x * G_x0);

    %% input uncertainty U
    eta_u = lookup.eta_u; 
    alpha_u = lookup.alpha_u; 
    G_u = diag(ones(dim_u, eta_u));
    c_u = lookup.c_u;
    c_delta_u = lookup.c_delta_u;
    c_u = c_u + c_delta_u;
    U = zonotope(c_u, alpha_u * G_u); 

    %% process noise W and its matrix-zonotope 
    c_w = lookup.c_w;
    alpha_w = lookup.alpha_w;
    eta_w = lookup.eta_w;
    G_w = ones(dim_x, eta_w);
    W = zonotope(c_w, alpha_w * G_w);
    W_mz = aux_buildMatrixZonotope(W, dim_x, totalsamples);

    %% sample all inputs (flattened over (m,s,t))
    u_all = zeros(dim_u, totalsamples);
    for idx = 1:totalsamples
        u_all(:, idx) = randPoint(U);
    end

    %% simulate: block per (m,s), each with n_k steps
    % Layout: block index b = (m-1)*n_s + s,  m=1..n_m, s=1..n_s
    n_blocks = n_m * n_s;

    % State trajectories stacked by blocks in rows, columns are time 1..(n_k+1)
    x_all = zeros(dim_x * n_blocks, n_k + 1);
    % Input trajectories stacked by blocks in rows, columns are time 1..n_k
    u_traj = zeros(dim_u * n_blocks, n_k);

    % initialize each block with a (possibly different) x0 sample from X0
    for b = 1:n_blocks
        row_x = (b-1)*dim_x + 1;
        x_all(row_x:row_x+dim_x-1, 1) = randPoint(X0);
    end

    % advance dynamics, consuming u_all and W in linear index order
    idx = 1;
    for b = 1:n_blocks
        row_x = (b-1)*dim_x + 1;
        row_u = (b-1)*dim_u + 1;
        for t = 1:n_k
            u_bt = u_all(:, idx);
            u_traj(row_u:row_u+dim_u-1, t) = u_bt;
            x_curr = x_all(row_x:row_x+dim_x-1, t);
            w_bt = randPoint(W);
            x_all(row_x:row_x+dim_x-1, t+1) = sys.A * x_curr + sys.B * u_bt + w_bt;
            idx = idx + 1;
        end
    end

    %% concatenate data: build X_minus, U_minus, X_plus in the same order as u_all
    X_minus = zeros(dim_x, n_blocks * n_k);
    U_minus = zeros(dim_u, n_blocks * n_k);
    X_plus = zeros(dim_x, n_blocks * n_k);

    col = 1;
    for b = 1:n_blocks
        row_x = (b-1)*dim_x + 1;
        row_u = (b-1)*dim_u + 1;
        % X_minus / U_minus at times 1..n_k, X_plus at times 2..n_k+1
        X_minus(:, col:col+n_k-1) = x_all(row_x:row_x+dim_x-1, 1:n_k);
        U_minus(:, col:col+n_k-1) = u_traj(row_u:row_u+dim_u-1, 1:n_k);
        X_plus(:, col:col+n_k-1) = x_all(row_x:row_x+dim_x-1, 2:n_k+1);
        col = col + n_k;
    end

    % truncate (already exact), kept for parity with original code
    U_full = U_minus(:, 1:totalsamples);
    X_0T = X_minus(:, 1:totalsamples);
    X_1T = X_plus(:, 1:totalsamples);

    % Compute M_AB 
    X1W_cen = X_1T - W_mz.center;
    X1W = matZonotope(X1W_cen, W_mz.generator);
    M_AB = X1W * pinv([X_0T; U_full]);

    % validate [A B] within interval matrix of M_AB
    box_Mab = intervalMatrix(M_AB);
    bounds = box_Mab.int;
    AB = [sys.A, sys.B];
    isContained = (bounds.sup >= AB) & (bounds.inf <= AB);
    if ~all(isContained(:))
        error('System matrices [A, B] are not contained within the interval matrix M_AB.');
    else
        disp('Validation passed: [A, B] is inside M_AB.');
    end

    %------------- END CODE --------------
end