function  conformance_linear(lookup, conf_opts)
    
    % example_linearSysDT_conform_04_LTV - example for conformance 
    %   identification of nonlinear discrete-time systems to analyze the 
    %   accuracy using simulated data
    % 
    %
    % Syntax:
    %    example_linearSysDT_conform_04_LTV
    %
    % Inputs:
    %    -
    %
    % Outputs:
    %    completed - true/false
    %
    % References:
    %    [1]
    
    % Authors:       Laura Luetzow
    % Written:       31-January-2024
    % Last update:   ---
    % Last revision: ---
    
    % ------------------------------ BEGIN CODE -------------------------------
    
    %% User Specifications ----------------------------------------------------
    % Set Random Number Stream
    rng(1, 'twister');
    
    dyn = lookup.sys.dyn; n_n = lookup.sys.n_n;

    n_s = lookup.n_s; % number of samples per test case
    n_m = lookup.n_m; % number of test cases for identification
    n_m_val = lookup.n_m_val;
    n_k = lookup.n_k; % number of time steps for identification
    n_k_val = lookup.n_k_val;

    constraints = {'gen','half',}; % type of containment constraints
    options = struct();
    
    % Conformance Settings
    options_reach = conf_opts.options_reach;
    options.options_reach = options_reach;
    options.cs = conf_opts.cs;
    options_testS = conf_opts.testS;
        
    % Evaluation settings
    check_contain = false;
    plot_settings.dims = [1 2];
    plot_settings.name = sprintf("Conformance Synthesis: %s", dyn);
    
    % Parameters and System Dynamics
    if dyn == "platoon"
        [sys, params_true.R0, params_true.U] = load_platoon(n_n,...
            max(n_k,n_k_val));
    else
        [sys, params_true.R0, params_true.U] = loadDynamics(dyn);
    end
    params_true.tFinal = sys.dt * n_k - sys.dt;
    
    % Simulation
    params_true.testSuite = createTestSuite(sys, ...
        params_true, n_k, n_m, n_s, options_testS);
    
    % Initial Estimates of the Disturbance Sets
    c_R0 = zeros(size(center(params_true.R0)));
    c_U = zeros(size(center(params_true.U)));
    params_id_init = params_true;
    params_id_init.R0 = zonotope([c_R0 eye(sys.nrOfDims)]);
    params_id_init.U = zonotope([c_U eye(sys.nrOfInputs)]);
    
    %% Conformance Identification ---------------------------------------------
    num_id = length(constraints);
    name_id = cell(num_id,2);
    
    % Struct for saving the identification results for each system
    configs = cell(num_id+1,1);
    configs{1}.sys = sys;
    configs{1}.params = params_true;
    configs{1}.options = options_reach;
    configs{1}.name = "true";
    
    for i_id = 1:num_id
        % run the identification
        options.cs.constraints = constraints{i_id};
        fprintf("Identification with %s-constraints, n_m=%d, " + ...
            "n_k=%d, n_x=%d\n",options.cs.constraints, n_m, n_k, 3*n_n);
        timerVal = tic;
        [configs{i_id+1}.params, ~] = conform(sys,params_id_init,options);
        configs{i_id+1}.sys = sys;
        configs{i_id+1}.options = options_reach;
        configs{i_id+1}.name = options.cs.constraints;
        Ts=toc(timerVal);
        fprintf("Identification time: %.4f\n", Ts);
    end
    
    %% Validation and Visualization -------------------------------------------
    
    % Create Validation Data
    if n_m_val ~= 0
        params_true.tFinal = sys.dt * n_k_val - sys.dt;
        testSuite_val = createTestSuite(sys, params_true, n_k_val, n_m_val, ...
            n_s, options);
    end
    % combine validation and trainings test cases
    testSuite{1} = combineTestCases(params_true.testSuite{1}, testSuite_val{1});
    plot_settings.s_val = size(params_true.testSuite{1}.y,3) + 1; 
        % (setting s_val leads to different color for validation test cases)
    
    % run validation and plotting
    validateReach(testSuite{1}, configs, check_contain, plot_settings);
    
    % example completed
    completed = true;
    
end
  

% ------------------------------ END OF CODE ------------------------------

