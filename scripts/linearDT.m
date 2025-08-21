%% 0. Define shared parameters
% 0.1 Define the system

dim_x = 5;
dim_u = 1;
dim_y = 1;
A = [-1 -4 0 0 0; 
    4 -1 0 0 0; 
    0 0 -3 1 0; 
    0 0 -1 -3 0; 
    0 0 0 0 -2];
B_ss = ones(dim_x,dim_u);
C = [1,0,0,0,0];
D = 0;
sys_c = ss(A,B_ss,C,D); % define continuous time system
dt = 0.05;
sys_d = c2d(sys_c,dt); % convert to discrete system


% 0.2. Define hyperparameters
lookup_ddra = struct( ...
    'n_s', 1, ... % Number of samples per unique trajectory
    'n_m', 1, ...  % Number of unique test cases (x0, u_{1:n_k})
    'n_k', 120, ... % Trajectory length per test case
    'eta_x', 1, ... % Number of generator vectors for \calX_0
    'alpha_x', 0.1, ... % Scaling factor for G_x0
    'c_x', ones(dim_x,1), ... % Initial guess for x0
    'c_delta_x', zeros(dim_x,1), ... % uncertainty about x0
    'eta_u', 1, ... % Number of generator vectors for \calU_i
    'alpha_u', 0.25, ... % Scaling factor for G_u
    'c_u', ones(dim_u, 1), ... % Initial guess for the control input(s)
    'c_delta_u', zeros(dim_u, 1), ... % Uncertainty about the initial guess c_u
    'c_w', zeros(dim_x, 1), ... % Center vector for process noise
    'alpha_w', 0.005, ... % Scaling factor for the generator matrix of process noise
    'eta_w', 1, ... % Number of generators for process noise
    'c_v', zeros(dim_x, 1), ... % Center vector for process noise
    'alpha_v', 0.0, ... % Scaling factor for the generator matrix of process noise
    'eta_v', 1, ... % Number of generators for process noise
    'dim_x', dim_x, ...
    'dim_y', dim_y, ...
    'dim_u', dim_u ...
);

[X_true, X_est] = a_linear(sys_d, lookup_ddra);

    