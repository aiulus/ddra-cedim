function R = conformance_gray(lookup, conf_opts)
    % CONFORMANCE_GRAY  Gray-box reachset-conformance wrapper.
    %   R = conformance_gray(lookup, conf_opts)
    %
    % Inputs 
    %   lookup: struct with fields
    %       .sys.dcheckyn   (string) dynamics key for custom_loadDynamics()
    %       .sys.n_n   (int)    optional (e.g., platoon vehicle count)
    %       .sys.type  (string) optional CORA uncertainty preset: "rand"|"diag"|"standard" (default "rand")
    %       .methodsGray (string array) optional, default ["graySeq","grayLS","graySim"]
    %       .n_s, .n_m, .n_k    (id data config)
    %       .n_m_val, .n_k_val, .n_s_val (validation config; n_s_val defaults to n_s)
    %
    %   conf_opts: struct with fields
    %       .options_reach (struct) reachability options
    %       .cs            (struct) base conformance options: robustnessMargin, cost, verbose
    %       .testS         (struct) test-suite options: p_extr, inputCurve (or for NL: just p_extr)
    %
    % Output:
    %   R: results from validateReach (containment stats + plots)
    %
    % Notes:
    %   - Falls back to white-box if p_true is empty and you only want to size U's center.
    %
    % Adapted from: example_nonlinearSysDT_conform_02_gray.m (CORA Toolbox)
    % Authors:       Laura Luetzow
    % Written:       14-September-2023
    % Author:  
    % Refactor: Aybüke Ulusarslan (22-Aug-2025)
    
    % ------------------------------ BEGIN CODE -------------------------------
    rng(1,'twister'); 
    
    % Unpack lookup
    sys   = lookup.sys;
    dyn = lookup.dyn;
    n_n   = getfieldwithdefault(lookup.sys, 'n_n', []);
    ldtyp = getfieldwithdefault(lookup.sys, 'type', "rand");
    
    n_s     = lookup.n_s;
    n_m     = lookup.n_m;
    n_k     = lookup.n_k;
    n_m_val = lookup.n_m_val;
    n_k_val = lookup.n_k_val;
    n_s_val = getfieldwithdefault(lookup, 'n_s_val', n_s);
    
    methodsGray = getfieldwithdefault(lookup, 'methodsGray', ["graySeq","grayLS","graySim"]);
    %constraints = getfieldwithdefault(lookup, 'constraints', "half"); % "half" or "gen"
    constraints = "half";
    
    % Optional toggles -- currently not in use
    %chk_id  = getfieldwithdefault(conf_opts, 'check_contain_id', 0); % default fast
    %chk_val = getfieldwithdefault(conf_opts, 'check_contain_val', 1); % default as in paper
    %do_plot = getfieldwithdefault(conf_opts, 'plot', false);
    plot_settings = getfieldwithdefault(conf_opts, 'plot_settings', struct());
    
    % Base options
    options_reach = conf_opts.options_reach;
    options_testS  = conf_opts.testS;
    
    if ~exist('p_true','var'); p_true = []; end
    params_true.tFinal = sys.dt * n_k - sys.dt;
    
    if sys.nrOfOutputs < 2
            constraints = {'gen'};
    end
    
    % Identification options (gray-box)
    options = options_reach;
    options.cs = conf_opts.cs;
    options.cs.constraints = constraints;
    
    if isfield(lookup, 'R0') 
        params_true.R0 = lookup.R0;
    end
    if isfield(lookup, 'U')
        params_true.U = lookup.U;
    end
    
    % Build identification test suite
    params_true.testSuite = createTestSuite(sys, params_true, n_k, n_m, n_s, options_testS);
    
    % derive dim_p from the model definition
    [~, ~, ~, p_true] = custom_loadDynamics(dyn, "diag");
    if isnumeric(p_true) && ~isempty(p_true)
        dim_p = numel(p_true);
    else
        dim_p = 0;
    end
    
    % p = [ model-parameters (if any) ; U-center ]
    options.cs.p0 = 0.01*randn(dim_p + sys.nrOfInputs, 1);
    options.cs.set_p = @(p,params) aux_set_p_gray(p, params, dyn, sys);
    
    % "Configs" container as in CORA examples
    configs = cell(numel(methodsGray) + 1,1);
    configs{1}.sys     = sys;
    configs{1}.params  = rmfield(params_true,'testSuite');
    configs{1}.options = options_reach;
    configs{1}.name    = "true";
    
    % Initial estimates for R0 and U (with slack generators like the example)
    c_R0 = center(params_true.R0);
    c_U  = center(params_true.U);
    params_id_init = params_true;
    params_id_init.R0 = zonotope([c_R0, eye(sys.nrOfDims),  ones(sys.nrOfDims,1)]);
    params_id_init.U  = zonotope([c_U,  eye(sys.nrOfInputs), ones(sys.nrOfInputs,1)]);
    
    % ---- Run gray-box identification (each method) ----
    for i = 1:numel(methodsGray)
        type = string(methodsGray(i));
        fprintf("Identification with method %s\n", type);
        t0 = tic;
        [configs{i+1}.params, results] = conform(sys, params_id_init, options, type);
        Ts = toc(t0);
        configs{i+1}.sys     = results.sys;
        configs{i+1}.options = options_reach;
        configs{i+1}.name    = type;
        fprintf("Identification time: %.4f s\n", Ts);
    end

    % Unify noise for validateReach: W := U  (since E=B, E*W = B*U)
    if sys.nrOfDisturbances > 0
        for i = 1:numel(configs)
            configs{i}.params.W = lookup.U;
        end
    end

    
    % ---- Validation & Visualization ----
    % Identification containment sanity-check on training data
    % ---- Validation & Visualization ----
    % Identification containment sanity-check on training data
    check_contain = getfieldwithdefault(conf_opts, 'check_contain', true);
    methodsList = ["true", methodsGray];
    
    isLinear = isa(sys,'linearSysDT') || isa(sys,'linearARX');
    
    if isLinear
        % Combine the whole training suite into one testCase
        tc_id = params_true.testSuite{1};
        for m = 2:numel(params_true.testSuite)
            tc_id = tc_id.combineTestCases(params_true.testSuite{m});
        end
        % Headless by default for ID check (no popups)
        [~, eval_id] = validateReach(tc_id, configs, check_contain);
   
        num_out = eval_id.num_out;                          % 1×(#configs)
        num_all = n_k * size(tc_id.y, 3);                   % (#steps)×(#traj)
    else
        % Nonlinear fallback: keep per-case, but headless to avoid popups
        num_out = zeros(1, numel(configs));
        num_all = 0;
        for m = 1:length(params_true.testSuite)
            [~, eval_id] = validateReach(params_true.testSuite{m}, configs, check_contain);
            num_out = num_out + eval_id.num_out;
            num_all = num_all + n_k * size(params_true.testSuite{m}.y,3);
        end
    end

    fprintf("IDENTIFICATION DATA:\n");
    for i = 1:length(configs)
        p_contained = 100 - (num_out(i)/num_all)*100;
        fprintf("%s: %.2f%% samples contained in reachable set (should be ~100%%).\n", methodsList(i), p_contained);
    end
    
    % Validation suite
    params_true.tFinal = sys.dt * n_k_val - sys.dt;
    testSuite_val = createTestSuite(sys, params_true, n_k_val, n_m_val, n_s_val, options_testS);
    
    plot_settings.plot_Yp = false;
    plot_settings.dims    = [1 2];
    plot_settings.name    = sprintf("Gray-Box Conformance: %s", dyn);

    do_plot = any(lower(string(getfielddef(cfg.io,'plot_mode','offline'))) == ["online","both"]);

    if do_plot
        PS = struct();
        PS.k_plot   = 1:C.shared.n_k_val;                     % <- plot all steps
        PS.dims     = getfielddef(cfg.io,'plot_dims',[1 2]);  % <- 2D slice
        PS.plot_Yp  = false;                                  % <- like the example
        ps_args = {'plot_settings', PS};
    else
        ps_args = {'plot_settings', []};  % truly suppress plotting
    end

    if ~isfield(plot_settings,'k_plot'), plot_settings.k_plot = 1:n_k_val; end
    if ~isfield(plot_settings,'dims'),   plot_settings.dims    = [1 2];    end
    if ~isfield(plot_settings,'plot_Yp'),plot_settings.plot_Yp = false;    end


    check_contain = getfieldwithdefault(conf_opts, 'check_contain', true);
    
    if isLinear
        % Single validateReach call for the whole validation suite
        tc_val = testSuite_val{1};
        for i = 2:numel(testSuite_val)
            tc_val = tc_val.combineTestCases(testSuite_val{i});
        end
        [R, eval_val] = validateReach(tc_val, configs, check_contain, plot_settings);
        num_out = eval_val.num_out;                      % 1×(#configs)
        num_all = n_k_val * size(tc_val.y, 3);
    else
        % Nonlinear fallback: per-case, but still just one figure per case if plotting
        num_out = zeros(1, numel(configs));
        num_all = 0;
        for i = 1:length(testSuite_val)
            [R, eval_val] = validateReach(testSuite_val{i}, configs, check_contain, plot_settings);
            num_out = num_out + eval_val.num_out;
            num_all = num_all + n_k_val * size(testSuite_val{i}.y,3);
        end
    end

    
    fprintf("VALIDATION DATA:\n");
    for i = 1:length(configs)
        p_contained = 100 - (num_out(i)/num_all)*100;
        fprintf("%s: %.2f%% samples contained in reachable set.\n", methodsList(i), p_contained);
    end
end

% ----- helpers -----
function val = getfieldwithdefault(S, fname, defaultVal)
    if isfield(S, fname); val = S.(fname); else; val = defaultVal; end
end

function [sys_out, params] = aux_set_p_gray(p, params, dyn, sys_base)
    % Probe if this dynamics has a numeric parameter vector
    [~, ~, ~, p_t] = custom_loadDynamics(dyn, "diag");

    if isnumeric(p_t) && ~isempty(p_t)
        % Numeric model params exist -> split [p_model; cU] and rebuild sys
        n_model = numel(p_t);
        p_model = p(1:n_model);
        cU      = p(n_model+1:end);
        [sys_out, ~, ~] = custom_loadDynamics(dyn, "diag", p_model);
    else
        % No numeric params (e.g., k-MSD): keep the current system unchanged
        sys_out = sys_base;
        cU      = p;  % all of p is the U-center
    end

    % Shift U center, keep generators
    params.U = zonotope(cU, params.U.generators);
end


function [sys_out, params] = aux_set_p_gray_functional2(p, params, dyn)
    [sys_def, ~, ~, p_t] = custom_loadDynamics(dyn, "diag"); 
    if isnumeric(p_t) && ~isempty(p_t)
        n_model = numel(p_t);
        assert(numel(p) >= n_model, ...
            'set_p: expected at least %d parameter entries, got %d.', n_model, numel(p));
        p_model = p(1:n_model);
        cU      = p(n_model+1:end);
        [sys_out, ~, ~] = custom_loadDynamics(dyn, "diag", p_model);
    else
        sys_out = sys_def;     % only U-center is optimized
        cU      = p;
    end
    params.U = zonotope(cU, params.U.generators);
end


function [sys_out, params] = aux_set_p_gray_functional1(p, params, dyn)
    % Query default system + the form of p_true without forcing parameters
    [sys_def, ~, ~, p_t] = custom_loadDynamics(dyn, "diag");  

    if isnumeric(p_t) && ~isempty(p_t)
        % We have a numeric parameter vector to identify (e.g., pedestrian, lorenz, NARX)
        n_model = numel(p_t);
        p_model = p(1:n_model);
        cU      = p(n_model+1:end);

        % Rebuild system with updated numeric parameters
        [sys_out, ~, ~] = custom_loadDynamics(dyn, "diag", p_model);
    else
        % No numeric parameter vector (e.g., k-Mass-SD uses a struct) — only shift U-center
        sys_out = sys_def;
        cU      = p;  % entire p is just U-center
    end

    % Shift the center of U, keep its generators
    params.U = zonotope(cU, params.U.generators);
end

function [sys_out, params] = aux_set_p_gray_original(p, params, dyn)
    % Rebuild the system with the first part of p as model parameters
    % (CORA's example uses 'diag' here as well)
    [sys_out, ~, ~, p_t] = custom_loadDynamics(dyn, "diag", p);

    % Remaining entries of p set the center of U
    cU = p(length(p_t) + 1 : end);

    % Keep U's generators, shift only its center
    params.U = zonotope(cU, params.U.generators);
end

% ------------------------------- END CODE --------------------------------