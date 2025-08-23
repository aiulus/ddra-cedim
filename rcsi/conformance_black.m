function R = conformance_black(lookup, conf_opts)
% CONFORMANCE_BLACK  Black-box reachset-conformance wrapper (GP/CGP).
%   R = conformance_black(lookup, conf_opts)
%
% Inputs:
%   lookup.sys.dyn, lookup.sys.type   : dynamics key + CORA loader preset ("rand"/"diag"/"standard")
%   Identification data config:        lookup.n_m, lookup.n_s, lookup.n_k
%   Black-box training suites:         lookup.n_m_train, lookup.n_s_train, lookup.n_k_train
%                                      lookup.n_m_val,   lookup.n_s_val,   lookup.n_k_val
%   lookup.methodsBlack                : string array; default ["blackGP","blackCGP"]
%   lookup.constraints                 : "half" | "gen" (default "half")
%   lookup.plot_settings               : struct with fields .plot_Yp, .dims, .name (optional)
%
%   conf_opts.options_reach            : reachability options
%   conf_opts.cs                       : conformance options (robustnessMargin, cost, verbose)
%   conf_opts.testS                    : test-suite options (e.g., .p_extr)
%   conf_opts.approx                   : GP/CGP options (gp_* and cgp_*, plus .p)
%
% Output:
%   R: output of last validateReach call (plots + containment stats)

% ------------------------------ BEGIN CODE -------------------------------
rng(1,'twister');

% ----- Unpack lookup + defaults -----
dyn    = lookup.sys.dyn;
ldtype = iff(isfield(lookup.sys,'type'), lookup.sys.type, "rand");

n_m     = lookup.n_m;     n_s     = lookup.n_s;     n_k     = lookup.n_k;
n_m_val = lookup.n_m_val; n_s_val = lookup.n_s_val; n_k_val = lookup.n_k_val;

n_m_train = iff(isfield(lookup,'n_m_train'), lookup.n_m_train, 100);
n_s_train = iff(isfield(lookup,'n_s_train'), lookup.n_s_train, 10);
n_k_train = iff(isfield(lookup,'n_k_train'), lookup.n_k_train, 4);

methodsBlack = string(iff(isfield(lookup,'methodsBlack'), lookup.methodsBlack, ["blackGP","blackCGP"]));
constraints  = string(iff(isfield(lookup,'constraints'),  lookup.constraints,  "half"));

options_reach = conf_opts.options_reach;
options_testS = conf_opts.testS;

% ----- GP/CGP options (with sensible defaults) -----
approx = struct();
approx.gp_parallel       = false;
approx.gp_pop_size       = 50;
approx.gp_num_gen        = 30;
approx.gp_func_names     = {'times','plus','square'};
approx.gp_max_genes      = 2;
approx.gp_max_depth      = 2;
approx.cgp_num_gen       = 5;
approx.cgp_pop_size_base = 5;
approx.save_res          = false;
approx.p                 = 0;   % set below from sys.n_p if ARX-like
if isfield(conf_opts,'approx')
    fn = fieldnames(conf_opts.approx);
    for i=1:numel(fn), approx.(fn{i}) = conf_opts.approx.(fn{i}); end
end

% ----- Load system + default sets -----
try
    [sys, params_true.R0, params_true.U, ~] = custom_loadDynamics(dyn, ldtype);
catch
    [sys, params_true.R0, params_true.U, ~] = loadDynamics(dyn, ldtype);
end
params_true.tFinal = sys.dt * n_k - sys.dt;

% ----- Build identification / train / validation suites -----
params_true.testSuite        = createTestSuite(sys, params_true, n_k,       n_m,      n_s,      options_testS);
params_true.testSuite_train  = createTestSuite(sys, params_true, n_k_train, n_m_train, n_s_train);
params_true.testSuite_val    = createTestSuite(sys, params_true, n_k_val,   n_m_val,   n_s_val, options_testS);

% ----- If ARX-like (nrOfDims==0), create a shell SS model with dim = dim_y * p -----
is_arx_like = (sys.nrOfDims == 0);
if is_arx_like
    dim_y  = sys.nrOfOutputs;
    dim_u  = sys.nrOfInputs;
    p_dim  = max(1, getfielddef(sys, 'n_p', 1));   % ARX order (default to 1)
    dim_x  = dim_y * p_dim;

    A_sh  = eye(dim_x);
    B_sh  = zeros(dim_x, dim_u);
    C_sh  = eye(dim_y, dim_x);
    D_sh  = zeros(dim_y, dim_u);
    dt_sh = sys.dt;

    sys_eff = linearSysDT(A_sh, B_sh, [], C_sh, D_sh, dt_sh);
    R0_eff  = zonotope(zeros(dim_x,1));

    % If user didn't set approx.p, make learned model stateful with same p
    if ~isfield(approx,'p') || isempty(approx.p) || approx.p==0
        approx.p = p_dim;
    end
else
    sys_eff = sys;                 % keep original stateful system
    R0_eff  = params_true.R0;
end

% ----- Identification options (conformance + approximator) -----
options            = options_reach;
options.cs         = conf_opts.cs;
options.cs.constraints = constraints;
options.approx     = approx;

% ----- Config container (like CORA examples) -----
configs = cell(numel(methodsBlack) + 1, 1);
configs{1}.sys     = sys;
configs{1}.params  = rmfield(params_true,'testSuite');   % keep R0, U, tFinal
configs{1}.options = options_reach;
configs{1}.name    = "true";

% ----- Initial params for conform(): keep all test suites; fix R0/U dims to sys_eff -----
params_id_init               = params_true;  % KEEP testSuite, *_train, *_val
params_id_init.R0            = R0_eff;

dim_u_eff = sys_eff.nrOfInputs;
cU        = center(params_true.U);
if numel(cU) ~= dim_u_eff, cU = zeros(dim_u_eff,1); end
params_id_init.U             = zonotope([cU(:), eye(dim_u_eff), ones(dim_u_eff,1)]);

% Optional: ensure GPTIPS 'extract.m' is preferred over MATLAB's built-in
if exist('prefer_gptips_extract','file') == 2
    prefer_gptips_extract();  %#ok<*UNRCH>
end

% ----- Run black-box identification (GP → CGP) -----
for i = 1:numel(methodsBlack)
    type = string(methodsBlack(i));
    fprintf("Identification with method %s \n", type);

    % Defensive re-harmonization before each call
    params_id_init.R0 = R0_eff;
    cU = center(params_true.U);
    if numel(cU) ~= dim_u_eff, cU = zeros(dim_u_eff,1); end
    params_id_init.U  = zonotope([cU(:), eye(dim_u_eff), ones(dim_u_eff,1)]);

    % -------------------- PATCH --------------------
    % ----- If ARX-like (nrOfDims==0), use the original system and 0-dim R0
    is_arx_like = (sys.nrOfDims == 0);
    
    sys_eff = sys;             % always use the real system for blackCGP
    if is_arx_like
        params_id_init.R0 = zonotope(zeros(0,1));     % key change
        if ~isfield(options,'approx') || ~isfield(options.approx,'p') || options.approx.p==0
            options.approx.p = max(1, getfielddef(sys,'n_p',1));
        end
    else
        params_id_init.R0 = params_true.R0;
    end
    
    % make sure U matches the system’s input dimension
    dim_u_eff = sys.nrOfInputs;
    cU = center(params_true.U);
    if numel(cU) ~= dim_u_eff, cU = zeros(dim_u_eff,1); end
    params_id_init.U = zonotope([cU(:), eye(dim_u_eff), ones(dim_u_eff,1)]);
    % -------------------- PATCH --------------------

    % Preflight assert
    assert(sys_eff.nrOfDims == dim(params_id_init.R0), 'R0/system dimension mismatch before conform().');

    t0 = tic;
    [configs{i+1}.params, results] = conform(sys_eff, params_id_init, options, type);
    Ts = toc(t0);

    configs{i+1}.sys     = results.sys;     % learned black-box model
    configs{i+1}.options = options_reach;
    configs{i+1}.name    = type;

    fprintf("Identification time: %.4f\n", Ts);
end

% ----- Validation & Visualization -----
% 1) Identification data sanity check
num_out = 0; num_in = 0; check_contain = 1;
methodsList = ["true", methodsBlack];
for m=1:length(params_true.testSuite)
    [~, eval_id] = validateReach(params_true.testSuite{m}, configs, check_contain);
    num_out = num_out + eval_id.num_out;
    num_in  = num_in  + eval_id.num_in;
end
num_all = length(params_true.testSuite)*n_k*size(params_true.testSuite{1}.y,3);
fprintf("IDENTIFICATION DATA:\n");
for i = 1:length(configs)
    p_contained = 100 - (num_out(i)/(num_out(i)+num_in(i)))*100;
    p_invalid   = 100*(num_all - (num_out(i)+num_in(i)))/num_all;
    fprintf("%s: %.2f%% contained. %.2f%% invalid.\n", methodsList(i), p_contained, p_invalid);
end

% 2) Validation suite (+ plots)
plot_settings = iff(isfield(lookup,'plot_settings'), lookup.plot_settings, struct('plot_Yp',false,'dims',[1 2]));
if ~isfield(plot_settings,'name')
    plot_settings.name = sprintf("Black-Box Conformance: %s", dyn);
end

num_out = 0; num_in = 0; check_contain = 1;
for m=1:length(params_true.testSuite_val)
    [R, eval_val] = validateReach(params_true.testSuite_val{m}, configs, check_contain, plot_settings);
    num_out = num_out + eval_val.num_out;
    num_in  = num_in  + eval_val.num_in;
end
num_all = length(params_true.testSuite_val)*n_k_val*size(params_true.testSuite_val{1}.y,3);
fprintf("VALIDATION DATA:\n");
for i = 1:length(configs)
    p_contained = 100 - (num_out(i)/(num_out(i)+num_in(i)))*100;
    p_invalid   = 100*(num_all - (num_out(i)+num_in(i)))/num_all;
    fprintf("%s: %.2f%% contained. %.2f%% invalid.\n", methodsList(i), p_contained, p_invalid);
end
end

% ---- tiny helpers ----
function v = iff(cond, a, b); if cond, v=a; else, v=b; end; end
function v = getfielddef(S,f,def); if isobject(S), has=ismethod(S,'properties'); else, has=isstruct(S); end
    if has && (isprop(S,f) || (isstruct(S)&&isfield(S,f))), v = S.(f); else, v = def; end
end
