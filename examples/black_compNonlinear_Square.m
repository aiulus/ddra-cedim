% black_compDT_square
%   DDRA linear demo (optional) + black-box conformance ID (GP/CGP).
%
% Prereqs:
%   - custom_loadDynamics.m (if you want to run the DDRA demo)
%   - cora_to_matlab_ss.m      (if you want to run the DDRA demo)
%   - ddra_linear.m            (if you want to run the DDRA demo)
%   - conformance_black.m      (this is the function we just built)

%% 0) Shared setup ---------------------------------------------------------
rng(1,'twister');  % reproducible outer script; conformance_black() uses its own seed

% -- System (CORA) -> MATLAB ss for DDRA demo (optional)
dyn_ddra  = "lorenz";
type_ddra = "standard";
% [sys_cora, ~, ~] = custom_loadDynamics(dyn_ddra, type_ddra);
% %% TODO conversion method from loadDynamics()
% %% sys_ddra = cora_to_matlab_ss(sys_cora);

% Example dims if you convert sys_cora:
% dim_x = sys_cora.nrOfStates;
% dim_u = sys_cora.nrOfInputs;
% dim_y = sys_cora.nrOfOutputs;

% DDRA plotting (optional)
plot_settings = struct();
plot_settings.projectedDims = {[1 2],[3 4]};

%% 1) DDRA Linear (optional, for side-by-side intuition) -------------------
% lookup_ddra = struct( ...
%     'n_s', 1, ...
%     'n_m', 1, ...
%     'n_k', 120, ...
%     'n_k_val', 20, ...
%     'eta_x', 1, ...
%     'alpha_x', 0.1, ...
%     'c_x', ones(dim_x,1), ...
%     'c_delta_x', zeros(dim_x,1), ...
%     'eta_u', 1, ...
%     'alpha_u', 0.25, ...
%     'c_u', ones(dim_u, 1), ...
%     'c_delta_u', zeros(dim_u, 1), ...
%     'c_w', zeros(dim_x, 1), ...
%     'alpha_w', 0.005, ...
%     'eta_w', 1, ...
%     'c_v', zeros(dim_x, 1), ...
%     'alpha_v', 0.0, ...
%     'eta_v', 1, ...
%     'dim_x', dim_x, ...
%     'dim_y', dim_y, ...
%     'dim_u', dim_u, ...
%     'plot_settings', plot_settings ...
% );
% % [X_true, X_est] = ddra_linear(sys_ddra, lookup_ddra);

%% 2) Black-box Reachset-Conformance ID -----------------------------------
% EXACTLY the same settings as the example_nonlinearARX_conform_03_black

dynamics     = "Square";
cost_norm    = "interval";  % size metric for reachable sets
constraints  = "half";      % constraint type in conformance fitting
methodsBlack = ["blackGP","blackCGP"];

% Sizes (ID/train/validation) – matches the example 1:1
sizes.id  = struct('n_m',  2, 'n_s', 50, 'n_k', 4);
sizes.tr  = struct('n_m',100, 'n_s', 10, 'n_k', 4);
sizes.val = struct('n_m',  5, 'n_s', 10, 'n_k', 4);

% Reachability options – same as example
options_reach = struct();
options_reach.zonotopeOrder     = 100;
options_reach.tensorOrder       = 2;
options_reach.errorOrder        = 1;
options_reach.tensorOrderOutput = 2;
options_reach.verbose           = false;

% testSuite settings – same as example
options_testS = struct('p_extr', 0.3);

% Plot settings – same default pattern
plot_settings_bb = struct('plot_Yp', false, 'dims', [1 2]);

% Black-box approximation options – identical to example
approx_options = struct();
approx_options.gp_parallel       = true;
approx_options.gp_pop_size       = 50;
approx_options.gp_num_gen        = 30;
approx_options.gp_func_names     = {'times','plus','square'};
approx_options.gp_max_genes      = 2;
approx_options.gp_max_depth      = 2;
approx_options.gp_parallel       = false;   % final value in the example
approx_options.cgp_num_gen       = 5;
approx_options.cgp_pop_size_base = 5;
approx_options.save_res          = false;

seed = 2;  % same RNG as in the example

% --- Run the black-box wrapper (prints stats and produces validation plots)
[configs, reports] = conformance_black( ...
    dynamics, cost_norm, constraints, methodsBlack, ...
    sizes, options_reach, options_testS, plot_settings_bb, approx_options, seed);

% 'configs' contains entries for "true", "blackGP", "blackCGP".
% 'reports' contains identification & validation percentages and counts.

% Optional: quick peek at validation containment percentages
disp('Validation containment (%):');
disp(table(reports.val.methods.', reports.val.p_contained, ...
    'VariableNames', {'Method','PercentContained'}));
