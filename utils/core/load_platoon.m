function [sys, R0, U] = load_platoon(N_v,N_k)
    % load the tank dynamics with the specified dimension--
    % called with n_n: #vehicles
    %             max(...
    %                 n_k: #time steps for identification, 
    %             n_k_val: #time test cases for idetification
    %                 )
    
    dt = 0.5; % Time step for discretization
    N_u = N_v; % #vehicles
    N_n = N_v*3; % 3*#vehicles
    sys = platoonN(dt,N_v,N_k);
    
    c_R0 = randn(N_n,1); 
    for i=0:N_v-1
        c_R0(i*N_v + 2) = 3*abs(c_R0(i*N_v + 2));
    end
    alpha_R0 = 2*rand(N_n,1);
    c_U = randn(N_u,1);
    alpha_U = rand(N_u,1);
    R0 = zonotope([c_R0,diag(alpha_R0)]);
    U = zonotope([c_U,diag(alpha_U)]);
    
end