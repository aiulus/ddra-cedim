function [X_model, X_data] = ddra_linear_inference(sys, lookup, M_AB)
    %% Expects:
    %   * sys = c2d(.)
    %   * lookup = struct(...) 
    %   * M_AB : learned matrix zonotope of [A B]
    %
    % Returns:
    %   * X_model, X_data (exactly as the original function)

    %------------- BEGIN CODE --------------

    %% system dimensions
    dim_x = lookup.dim_x; 
    dim_u = lookup.dim_u; 
    dim_y = lookup.dim_y; 

    n_k_val = lookup.n_k_val; % length of the evaluation trajectory

    %% initial-state uncertainty X0
    eta_x = lookup.eta_x; 
    alpha_x = lookup.alpha_x; 
    G_x0  = diag(ones(dim_x, eta_x));
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

    %% process noise W
    c_w = lookup.c_w;
    alpha_w = lookup.alpha_w;
    eta_w = lookup.eta_w; 
    G_w = ones(dim_x, eta_w);
    W = zonotope(c_w, alpha_w * G_w);

    %% propagate sets (unchanged)
    X_model = cell(n_k_val+1,1);
    X_data = cell(n_k_val+1,1);
    X_model{1} = X0; 
    X_data{1} = X0;

    for i = 1:n_k_val
        X_model{i,1} = reduce(X_model{i,1}, 'girard', 400);
        X_model{i+1,1}= sys.A * X_model{i} + sys.B * U + W;

        X_data{i,1} = reduce(X_data{i,1}, 'girard', 400);
        X_data{i+1,1} = M_AB * (cartProd(X_data{i}, U)) + W;
    end

    %% visualization (unchanged)
    projectedDims = {[1 2],[3 4],[4 5]};
    aux_visualize_original(X0, X_model, X_data, projectedDims);

    %------------- END CODE --------------
end