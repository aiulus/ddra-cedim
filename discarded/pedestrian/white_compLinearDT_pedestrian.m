%% 0. Define shared parameters
% 0.1 Define the system
% ---- DDRA Std. ----
dyn = "pedestrian";
type = "standard";
[sys_cora, ~, ~] = custom_loadDynamics(dyn, type); % define continuous time system
sys_ddra = cora_to_matlab_ss(sys_cora);
dim_x = size(sys_ddra.A, 1);
dim_u = size(sys_ddra.B, 2);
dim_y = size(sys_ddra.C, 1);

projectedDims = {[1 2],[3 4]};
plot_settings = struct();
plot_settings.projectedDims = projectedDims;

% 0.2. Define hyperparameters
lookup_ddra = struct( ...
    'n_s', 1, ... % Number of samples per unique trajectory
    'n_m', 1, ...  % Number of unique test cases (x0, u_{1:n_k})
    'n_k', 120, ... % Trajectory length per test case
    'n_k_val', 20, ... % Validation suite / length of the evaluation trajectory
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
    'dim_u', dim_u, ...
    'plot_settings', plot_settings ...
);

[X_true, X_est] = ddra_linear(sys_ddra, lookup_ddra);

%% Rechset-Conformant System Identification
dyn = "pedestrian"; % dynamics (choose from "platoon", "pedestrian")

lookup_rcsi = struct( ...
    'sys', struct('dyn', dyn), ...
    'n_s', 1, ... % Number of samples per unique trajectory
    'n_m', 20, ...  % Number of unique test cases (x0, u_{1:n_k})
    'n_m_val', 100, ... % Validation suite
    'n_k', 8, ... % Trajectory length per test case
    'n_k_val', 20, ... % Validation suite / length of the evaluation trajectory
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

conf_opts_white = struct( ...
    'options_reach', struct('zonotopeOrder', 100), ...
    'cs', struct('robustnessMargin', 1e-9, 'cost', "interval", 'verbose', false), ...
    'testS', struct('p_extr', 0.2, 'inputCurve', "rand") ... % (Pr. extreme points, traj. type for input)
);

if dyn == "platoon"
    [sys, params_true.R0, params_true.U] = load_platoon(lookup_rcsi.sys.n_n,...
        max(lookup_rcsi.n_k, lookup_rcsi.n_k_val));
else
    [sys, params_true.R0, params_true.U] = custom_loadDynamics(dyn);
end

R = conformance_white(lookup_rcsi, conf_opts_white);