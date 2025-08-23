function [X_model, X_data] = ddra_linear(sys, lookup)
    %% DDRA_LINEAR  Algorithm 1 (LTI Reachability) 
    %   [X_model, X_data] = ddra_linear(sys, lookup)
    %
    % Expects:
    %   * sys  : discrete-time LTI (from c2d or equivalent) with fields A,B,C
    %   * lookup: struct of options (see below)
    %
    % Changes vs. original:
    %   (1) RNG via rng(1,'twister') for reproducibility parity
    %   (2) Matrix zonotope Mw center = [c_w ... c_w]; generators circulated
    %   (3) Subtraction X^+ - Mw uses NEGATED generators (Alg. 1, eq. 23)
    %   (4) Rank check on [X^-; U^-] before pseudoinverse; optional ridge
    %
    % Refactor: 23-Aug-2025
    
    %------------- BEGIN CODE --------------
    
    rng(lookup.seed_if_present_or(1), 'twister'); % reproducible
    
    %% system dimensions
    if isfield(lookup,'dim_x'), dim_x = lookup.dim_x; else, dim_x = size(sys.A,1); end
    if isfield(lookup,'dim_u'), dim_u = lookup.dim_u; else, dim_u = size(sys.B,2); end
    if isfield(lookup,'dim_y'), dim_y = lookup.dim_y; else, dim_y = size(sys.C,1); end
    
    n_s = lookup.n_s;   % samples per test case
    n_m = lookup.n_m;   % number of test cases
    n_k = lookup.n_k;   % trajectory length per sample
    n_k_val = lookup.n_k_val; % evaluation trajectory length
    T = n_m * n_s * n_k;  % total time-samples across all (m,s)
    
    %% initial-state uncertainty X0 (same paramization as original)
    eta_x   = lookup.eta_x; 
    alpha_x = lookup.alpha_x; 
    G_x0    = diag(ones(dim_x, eta_x));
    c_x         = lookup.c_x;
    c_delta_x   = lookup.c_delta_x;
    X0 = zonotope(c_x + c_delta_x, alpha_x * G_x0);
    
    %% input uncertainty U
    eta_u   = lookup.eta_u; 
    alpha_u = lookup.alpha_u; 
    G_u     = diag(ones(dim_u, eta_u));
    c_u       = lookup.c_u;
    c_delta_u = lookup.c_delta_u;
    U = zonotope(c_u + c_delta_u, alpha_u * G_u);
    
    %% process noise W and its matrix-zonotope (paper-faithful)
    c_w     = lookup.c_w;
    alpha_w = lookup.alpha_w;
    eta_w   = lookup.eta_w;
    G_w     = ones(dim_x, eta_w);
    W       = zonotope(c_w, alpha_w * G_w);
    Mw      = build_Mw_matrix_zonotope(W, dim_x, T);  % Mw = <CMw, G̃Mw>
    
    %% sample all inputs (flattened over (m,s,t))
    u_all = zeros(dim_u, T);
    for idx = 1:T
        u_all(:, idx) = randPoint(U);
    end
    
    %% simulate: block per (m,s), each with n_k steps
    n_blocks = n_m * n_s;
    x_all    = zeros(dim_x * n_blocks, n_k + 1); % states per block
    u_traj   = zeros(dim_u * n_blocks, n_k);     % inputs per block
    
    % init each block with a sample from X0
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
            w_bt   = randPoint(W);
            x_all(row_x:row_x+dim_x-1, t+1) = sys.A * x_curr + sys.B * u_bt + w_bt;
            idx = idx + 1;
        end
    end
    
    %% concatenate data: build X_minus, U_minus, X_plus in same order as u_all
    X_minus = zeros(dim_x, n_blocks * n_k);
    U_minus = zeros(dim_u, n_blocks * n_k);
    X_plus  = zeros(dim_x, n_blocks * n_k);
    
    col = 1;
    for b = 1:n_blocks
        row_x = (b-1)*dim_x + 1;
        row_u = (b-1)*dim_u + 1;
        X_minus(:, col:col+n_k-1) = x_all(row_x:row_x+dim_x-1, 1:n_k);
        U_minus(:, col:col+n_k-1) = u_traj(row_u:row_u+dim_u-1, 1:n_k);
        X_plus(:,  col:col+n_k-1) = x_all(row_x:row_x+dim_x-1, 2:n_k+1);
        col = col + n_k;
    end
    
    %% Identify M_AB as in Algorithm 1
    % Form (X^+ - Mw) = <X_plus - C_Mw,  -G̃_Mw>
    X1W_center = X_plus - Mw.center;                 % center subtraction
    negG = cellfun(@(G) -G, Mw.generator, 'UniformOutput', false); % generator negation
    X1W = matZonotope(X1W_center, negG);
    
    % Pseudoinverse multiplier
    Z = [X_minus; U_minus];
    % Rank check (full row rank recommended)
    full_row = size(Z,1);
    rankZ = rank(Z);
    if rankZ < min(full_row, size(Z,2))
        warning('ddra_linear:rankDeficient', ...
            '[X^-; U^-] is rank-deficient (rank=%d < %d). Adding Tikhonov ridge.', rankZ, min(full_row,size(Z,2)));
        lambda = 1e-8;
        M_AB = X1W * ((Z'/(Z*Z' + lambda*eye(full_row)))'); % (Z)^\dagger with ridge
    else
        M_AB = X1W * pinv(Z);
    end
    
    % Validate [A B] inside interval matrix of M_AB
    box_Mab   = intervalMatrix(M_AB);
    bounds    = box_Mab.int;
    AB        = [sys.A, sys.B];
    isContained = (bounds.sup >= AB) & (bounds.inf <= AB);
    if ~all(isContained(:))
        error('System matrices [A, B] are not contained within the interval matrix M_AB.');
    else
        disp('Validation passed: [A, B] is inside M_AB.');
    end
    
    %% propagate sets (unchanged)
    X_model = cell(n_k_val+1,1);
    X_data  = cell(n_k_val+1,1);
    X_model{1} = X0; 
    X_data{1}  = X0;
    
    for i = 1:n_k_val
        X_model{i,1}   = reduce(X_model{i,1}, 'girard', 400);
        X_model{i+1,1} = sys.A * X_model{i} + sys.B * U + W;
    
        X_data{i,1}    = reduce(X_data{i,1}, 'girard', 400);
        X_data{i+1,1}  = M_AB * (cartProd(X_data{i}, U)) + W;
    end
    
    %% visualization
    projectedDims = lookup.plot_settings.projectedDims;
    aux_visualize_original(X0, X_model, X_data, projectedDims);
end

%--------------------------- END OF CODE ----------------------------

function MZ = build_Mw_matrix_zonotope(W, dim_x, T)
    % BUILD_MW_MATRIX_ZONOTOPE  Mw = <CMw, G̃Mw> with center repetition and
    % generator "circulation" across T columns (paper-faithful)
    CMw = repmat(W.center, 1, T);
    Z = [W.center, W.generators];
    nGen = size(W.generators, 2);
    GW = cell(1, nGen * T);
    idx = 1;
    for i = 1:nGen
        g = Z(:, i+1);
        G0 = [g, zeros(dim_x, T-1)];
        GW{idx} = G0;
        for j = 1:(T-1)
            GW{idx + j} = [GW{idx + j - 1}(:, 2:end), GW{idx + j - 1}(:, 1)];
        end
        idx = idx + T;
    end
    MZ = matZonotope(CMw, GW);
end

function seed = seed_if_present_or(default)
    % helper to keep backward compatibility with older lookups that did not
    % carry a seed; prefer cfg.shared.seed in caller
    seed = default;
end

%% ---------------- Visualization copied from original ------------------
function aux_visualize_original(X0, X_model, X_data, projectedDims)
    axx{1} = [0.75,1.5,0.5,4]; 
    axx{2} = [0.75,3,0.8,2.2];
    axx{3} = [0.75,2.3,0.75,2.8];
    index=1;
    numberofplots = 5; % length(X_model)
    for plotRun=1:length(projectedDims)
        figure('Renderer', 'painters', 'Position', [10 10 700 900])
        index=index+1;
        % plot initial set
        handleX0 = plot(X0,projectedDims{plotRun},'k-','LineWidth',2);
        hold on;
        % plot reachable sets from model
        for iSet=2:numberofplots
            handleModel=  plot(X_model{iSet},projectedDims{plotRun},'b','Filled',true,'FaceColor',[.8 .8 .8],'EdgeColor','b'); 
        end
        % plot reachable sets from data
        for iSet=2:numberofplots
            handleData=   plot(X_data{iSet},projectedDims{plotRun},'r');
        end
        xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
        ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
        warOrig = warning; warning('off','all');
        legend([handleX0,handleModel,handleData],...
            'Initial Set','Set from Model','Set from Data','Location','northwest');
        warning(warOrig);
        ax = gca; ax.FontSize = 22;
        outerpos = ax.OuterPosition; ti = ax.TightInset;
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
    end
end
