function OUT = kMSD_runner(cfg)
    % SCALABILITY_KMSD  Side-by-side DDRA (Alg.1) vs Gray-Conformance (Algs.2&3)
    %   OUT = scalability_kMSD(cfg)
    %
    % Inputs (struct cfg)
    %   .shared : struct of shared hyperparameters (identical across pipelines)
    %   .ddra   : struct of DDRA-specific knobs
    %   .gray   : struct of gray-conformance-specific knobs
    %
    % This stub builds both pipelines under *equivalent* conditions, runs them,
    % collects results/timings, and saves plots + .mat artifacts under ./results
    % and ./plots. 
    %
    % Return struct OUT with fields:
    %   .cfg, .timing, .ddra, .gray
    %
    % Notes
    %   - Designed for k-Mass-SD (kMSD) use, but dynamics are configurable
    %   - Uses the robust aux_set_p_gray() shape-aware split 
    %
    
    %% ---------- 0) Defaults & I/O ----------
    if nargin==0 || isempty(cfg)
        cfg = default_config();
    end
    
    if ~isfield(cfg,'io') || ~isfield(cfg.io,'save_tag')
        cfg.io.save_tag = datestr(now,'yyyymmdd_HHMMSS');
    end
    plots_dir   = fullfile('plots', cfg.io.save_tag);
    results_dir = fullfile('results', cfg.io.save_tag);
    if ~exist(plots_dir, 'dir'),   mkdir(plots_dir);   end
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
    
    rng(cfg.shared.seed, 'twister');
    
    %% ---------- 1) Build systems (CORA + MATLAB ss) ----------
    % CORA system + sets (for gray pipeline)
    [sys_cora, R0_cora, U_cora] = custom_loadDynamics(cfg.shared.dyn, cfg.shared.type);  % true system; do NOT pass p_init here
    
    % Carry optional W/V for parity if present in cfg.shared
    if isfield(cfg.shared,'W'), W = cfg.shared.W; else, W = []; end
    if isfield(cfg.shared,'V'), V = cfg.shared.V; else, V = []; end
    
    % MATLAB ss for DDRA (same dynamics)
    sys_ddra = cora_to_matlab_ss(sys_cora);  
    
    % Unify dimensions from systems
    dim_x = size(sys_ddra.A,1);
    dim_u = size(sys_ddra.B,2);
    dim_y = size(sys_ddra.C,1);
    
    %% ---------- 2) Shared hyperparameters -> structs ----------
    shared = struct();
    shared.seed    = cfg.shared.seed;
    shared.n_m     = cfg.shared.n_m;
    shared.n_s     = cfg.shared.n_s;
    shared.n_k     = cfg.shared.n_k;
    shared.n_m_val = cfg.shared.n_m_val;
    shared.n_s_val = cfg.shared.n_s_val;
    shared.n_k_val = cfg.shared.n_k_val;
    shared.dim_x   = dim_x; shared.dim_u = dim_u; shared.dim_y = dim_y;
    
    % Use *exact same* R0, U as the CORA loader returned
    shared.R0 = R0_cora;
    shared.U  = U_cora;
    if ~isempty(W), shared.W = W; end
    if ~isempty(V), shared.V = V; end
    
    % Reachability opts for CORA
    shared.options_reach = cfg.shared.options_reach;
    
    % Test-suite shaping to mirror DDRA’s i.i.d. per-step input sampling
    % Test-suite options: choose between CORA-default behavior and DDRA-like i.i.d. per-step inputs
    if isfield(cfg.shared,'testSuite_mode') && cfg.shared.testSuite_mode=="ddra_like"
        shared.options_testS = struct('p_extr', cfg.shared.p_extr, ...
                                      'inputCurve', "randn", ...
                                      'contInput', false, ...
                                      'stateSet', R0_cora);
    else
        % CORA defaults (Algorithm 2/3 examples): randn curve with cumulative integrator
        shared.options_testS = struct('p_extr', cfg.shared.p_extr);
    end
    
    % Base conformance options
    shared.cs_base = cfg.shared.cs_base;  % robustnessMargin, cost, constraints, verbose
    
    %% ---------- 3) Run DDRA pipeline ----------
    lookup_ddra = struct( ...
        'n_s', shared.n_s, 'n_m', shared.n_m, 'n_k', shared.n_k, ...
        'n_k_val', shared.n_k_val, ...
        'eta_x', cfg.ddra.eta_x, 'alpha_x', cfg.ddra.alpha_x, ...
        'c_x', center(shared.R0), 'c_delta_x', zeros(dim_x,1), ...
        'eta_u', cfg.ddra.eta_u, 'alpha_u', cfg.ddra.alpha_u, ...
        'c_u', center(shared.U), 'c_delta_u', zeros(dim_u,1), ...
        'c_w', get_center_if(shared,'W', zeros(dim_x,1)), ...
        'alpha_w', cfg.ddra.alpha_w, 'eta_w', cfg.ddra.eta_w, ...
        'c_v', zeros(dim_x,1), 'alpha_v', 0.0, 'eta_v', 1, ...
        'dim_x', dim_x, 'dim_y', dim_y, 'dim_u', dim_u, ...
        'plot_settings', cfg.ddra.plot_settings);
    
    % Run & time
    T0 = tic;
    [X_model, X_data] = ddra_linear(sys_ddra, lookup_ddra);
    Ts_ddra = toc(T0);
    
    % Save DDRA figures
    save_all_open_figs(plots_dir, 'ddra');
    
    %% ---------- 4) Run Gray-Conformance pipeline ----------
    % Build lookup for gray run (use same R0/U and sampling settings; Algorithms 2 & 3)
    lookup_gray = struct( ...
        'sys', sys_cora, 'dyn', cfg.shared.dyn, ...
        'R0', shared.R0, 'U', shared.U, ...
        'n_m', shared.n_m, 'n_s', shared.n_s, 'n_k', shared.n_k, ...
        'n_m_val', shared.n_m_val, 'n_s_val', shared.n_s_val, 'n_k_val', shared.n_k_val, ...
        'methodsGray', cfg.gray.methodsGray, ...
        'constraints', shared.cs_base.constraints);
    
    % Pass W/V for parity (CORA will augment U internally)
    if isfield(shared,'W'), lookup_gray.W = shared.W; end
    if isfield(shared,'V'), lookup_gray.V = shared.V; end
    
    conf_opts_gray = struct('options_reach', shared.options_reach, ...
                            'cs', shared.cs_base, ...
                            'testS', shared.options_testS);
    
    T0 = tic;
    R_gray = conformance_gray(lookup_gray, conf_opts_gray); 
    Ts_gray = toc(T0);
    
    % Save gray figures
    save_all_open_figs(plots_dir, 'gray');
    
    %% ---------- 5) Store artifacts ----------
    OUT = struct();
    OUT.cfg    = cfg;
    OUT.timing = struct('ddra_s', Ts_ddra, 'gray_s', Ts_gray);
    OUT.ddra   = struct('X_model', {X_model}, 'X_data', {X_data});
    OUT.gray   = struct('note', 'See plots + conformance_gray outputs');
    
    save(fullfile(results_dir, 'scalability_OUT.mat'), 'OUT');
    
    fprintf('DDRA time: %.3f s | Gray time: %.3f s\n', Ts_ddra, Ts_gray);
    fprintf('Saved results -> %s\nSaved plots   -> %s\n', results_dir, plots_dir);
end

%% ===================== Helpers & Runners =====================
function cfg = default_config()
    cfg = struct();
    
    % --- shared ---
    cfg.shared = struct();
    cfg.shared.seed = 1;
    cfg.shared.dyn  = "k-Mass-SD";  
    cfg.shared.type = "diag";        
    cfg.shared.p_extr = 0.3;
    
    % horizons
    cfg.shared.n_m = 5;  cfg.shared.n_s = 30; cfg.shared.n_k = 4;
    cfg.shared.n_m_val = 2; cfg.shared.n_s_val = 30; cfg.shared.n_k_val = 4;
    
    % CORA reachability options
    cfg.shared.options_reach = struct('zonotopeOrder',100,'tensorOrder',2, ...
        'errorOrder',1,'tensorOrderOutput',2,'verbose',false);
    
    % White/Gray LP base
    cfg.shared.cs_base = struct('robustnessMargin',1e-9,'cost',"interval", ...
        'verbose', false, 'constraints', "half");
    
    % Optional W/V parity with DDRA (unset by default)
    % cfg.shared.W = zonotope(zeros(?));
    % cfg.shared.V = zonotope(zeros(?));
    
    % --- ddra ---
    cfg.ddra = struct();
    cfg.ddra.eta_x = 1; cfg.ddra.alpha_x = 0.1;
    cfg.ddra.eta_u = 1; cfg.ddra.alpha_u = 0.25;
    cfg.ddra.eta_w = 1; cfg.ddra.alpha_w = 0.005;
    % plotting
    cfg.ddra.plot_settings = struct('projectedDims', {{[1 2],[3 4]}});
    
    % --- gray ---
    cfg.gray = struct();
    %cfg.gray.methodsGray = ["graySeq","grayLS","graySim"];
    cfg.gray.methodsGray = ["grayLS"];

    % initial numeric parameter guess (only used by gray pipeline during optimization)
    cfg.gray.p_init = [];  % for k-Mass-SD (struct params) leave empty
    
    % --- io ---
    cfg.io = struct('save_tag','kMSD_default');
end

function c = get_center_if(shared, field, fallback)
    if isfield(shared, field)
        c = center(shared.(field));
    else
        c = fallback;
    end
end

function save_all_open_figs(dir_out, prefix)
    figs = findall(0,'Type','figure');
    for i=1:numel(figs)
        f = figs(i);
        nm = sprintf('%s_%02d', prefix, i);
        exportgraphics(f, fullfile(dir_out, nm + ".pdf"), 'ContentType','vector');
    end
end
