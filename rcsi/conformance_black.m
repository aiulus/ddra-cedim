function [configs, reports] = conformance_black(dynamics, cost_norm, constraints, methodsBlack, sizes, options_reach, options_testS, plot_settings, approx_options, seed)
% conformance_black - Black-box conformance identification & validation
%   (GP / CGP surrogates), factored from CORA example_nonlinearARX_conform_03_black.
%
% Usage (reproduce the example 1:1):
%   [configs, reports] = conformance_black();
%
% Custom usage (override pieces):
%   sizes.id  = struct('n_m',2,'n_s',50,'n_k',4);
%   sizes.tr  = struct('n_m',100,'n_s',10,'n_k',4);
%   sizes.val = struct('n_m',5,'n_s',10,'n_k',4);
%   [configs, reports] = conformance_black("Square","interval","half",["blackGP","blackCGP"], ...
%                             sizes, [], struct('p_extr',0.3), struct('plot_Yp',false,'dims',[1 2]), [], 2);
%
% Returns
%   configs : cell array with entries {true, blackGP, blackCGP, ...}
%             Each entry has fields: .sys, .params, .options, .name
%   reports : struct with identification/validation summaries (counts & %)
%
% Notes
%   - Defaults match the provided example exactly (RNG, sizes, norms, etc.).
%   - Plots are produced in the validation step, just like the example.

% ------------------------------ Defaults ---------------------------------
if nargin < 1 || isempty(dynamics),    dynamics    = "Square";          end
if nargin < 2 || isempty(cost_norm),   cost_norm   = "interval";        end
if nargin < 3 || isempty(constraints), constraints = "half";            end
if nargin < 4 || isempty(methodsBlack),methodsBlack= ["blackGP","blackCGP"]; end
if nargin < 5 || isempty(sizes)
    sizes.id  = struct('n_m',2,   'n_s',50, 'n_k',4);
    sizes.tr  = struct('n_m',100, 'n_s',10, 'n_k',4);
    sizes.val = struct('n_m',5,   'n_s',10, 'n_k',4);
end
if nargin < 6 || isempty(options_reach)
    options_reach.zonotopeOrder      = 100;
    options_reach.tensorOrder        = 2;
    options_reach.errorOrder         = 1;
    options_reach.tensorOrderOutput  = 2;
    options_reach.verbose            = false;
end
if nargin < 7 || isempty(options_testS), options_testS = struct('p_extr',0.3); end
if nargin < 8 || isempty(plot_settings), plot_settings = struct('plot_Yp',false,'dims',[1 2]); end
if nargin < 9 || isempty(approx_options)
    % Defaults mirror the example (including gp_parallel toggle to false at end)
    approx_options.gp_parallel      = true;
    approx_options.gp_pop_size      = 50;
    approx_options.gp_num_gen       = 30;
    approx_options.gp_func_names    = {'times','plus','square'};
    approx_options.gp_max_genes     = 2;
    approx_options.gp_max_depth     = 2;
    approx_options.gp_parallel      = false; % final value in the example
    approx_options.cgp_num_gen      = 5;
    approx_options.cgp_pop_size_base= 5;
    approx_options.save_res         = false;
end
if nargin < 10 || isempty(seed), seed = 2; end

% ---------------------------- Begin Code ----------------------------------
rng(seed);

% Load system & base uncertainty sets (exactly as in example)
[sys, R0, U, p_true] = loadDynamics(dynamics,"rand");
params_true.R0 = R0;
params_true.U  = U;
params_true.tFinal = sys.dt * sizes.id.n_k - sys.dt;

% Create identification/train/validation suites (as in example)
params_true.testSuite       = createTestSuite(sys, params_true, sizes.id.n_k,  sizes.id.n_m,  sizes.id.n_s,  options_testS);
params_true.testSuite_train = createTestSuite(sys, params_true, sizes.tr.n_k,  sizes.tr.n_m,  sizes.tr.n_s);
params_true.testSuite_val   = createTestSuite(sys, params_true, sizes.val.n_k, sizes.val.n_m, sizes.val.n_s);

% Identification options (copy of example behavior)
options = options_reach;
options.cs.robustnessMargin = 1e-9;
options.cs.verbose          = false;
options.cs.cost             = cost_norm;
options.cs.constraints      = constraints;

% Black-box approximation options (copy into options.approx and set p)
options.approx = approx_options;
options.approx.p = sys.n_p;

% Set up configs cell with "true" baseline first
configs = cell(length(methodsBlack) + 1,1);
configs{1}.sys     = sys;
configs{1}.params  = rmfield(params_true,'testSuite');
configs{1}.options = options_reach;
configs{1}.name    = "true";

% Initial disturbance set estimates (exactly as in example)
c_R0 = center(params_true.R0);
c_U  = center(params_true.U);
params_id_init = params_true;
params_id_init.R0 = zonotope([c_R0]);
params_id_init.U  = zonotope([c_U eye(size(c_U,1)) ones(size(c_U))]);

% ------------------------ Identification loop ----------------------------
for i = 1:length(methodsBlack)
    type = methodsBlack(i);
    fprintf("Identification with method %s \n", type);
    timerVal = tic;
    [configs{i+1}.params, results] = conform(sys, params_id_init, options, type);
    Ts = toc(timerVal);
    configs{i+1}.sys     = results.sys;
    configs{i+1}.options = options_reach;
    configs{i+1}.name    = type;
    fprintf("Identification time: %.4f\n", Ts);
end

% ---------------------- Identification data check -------------------------
num_out = zeros(length(configs),1);
num_in  = zeros(length(configs),1);
check_contain = 1;
methods = ["true" methodsBlack];

for m = 1:length(params_true.testSuite)
    [~, eval] = validateReach(params_true.testSuite{m}, configs, check_contain);
    num_out = num_out + eval.num_out;
    num_in  = num_in  + eval.num_in;
end

% NOTE: We keep the example's reporting formula for num_all (uses n_k_val).
num_all_ident = length(params_true.testSuite) * sizes.val.n_k * size(params_true.testSuite{1}.y,3);

fprintf("IDENTIFICATION DATA: \n");
p_contained_ident = zeros(length(configs),1);
p_invalid_ident   = zeros(length(configs),1);
for i = 1:length(configs)
    p_contained_ident(i) = 100 - (num_out(i) / (num_out(i) + num_in(i))) * 100;
    p_invalid_ident(i)   = (num_all_ident - (num_out(i) + num_in(i))) / num_all_ident * 100;
    fprintf("%s: %.2f%% of the samples are contained in the reachable set (must be 100%%!). \n", ...
        methods(i), p_contained_ident(i));
    fprintf("%s: %.2f%% of the samples were not valid (measurement was nan or reachable set could not be computed). \n", ...
        methods(i), p_invalid_ident(i));
end

% ---------------------- Validation data & plots ---------------------------
params_true.tFinal = sys.dt * sizes.val.n_k - sys.dt;
testSuite_val = createTestSuite(sys, params_true, sizes.val.n_k, sizes.val.n_m, sizes.val.n_s, options_testS);

num_out_val = zeros(length(configs),1);
num_in_val  = zeros(length(configs),1);
check_contain = 1;

for m = 1:length(testSuite_val)
    [~, eval] = validateReach(testSuite_val{m}, configs, check_contain, plot_settings);
    num_out_val = num_out_val + eval.num_out;
    num_in_val  = num_in_val  + eval.num_in;
end

num_all_val = length(testSuite_val) * sizes.val.n_k * size(testSuite_val{1}.y,3);

fprintf("VALIDATION DATA: \n");
p_contained_val = zeros(length(configs),1);
p_invalid_val   = zeros(length(configs),1);
for i = 1:length(configs)
    p_contained_val(i) = 100 - (num_out_val(i) / (num_out_val(i) + num_in_val(i))) * 100;
    p_invalid_val(i)   = (num_all_val - (num_out_val(i) + num_in_val(i))) / num_all_val * 100;
    fprintf("%s: %.2f%% of the samples are contained in the reachable set. \n", ...
        methods(i), p_contained_val(i));
    fprintf("%s: %.2f%% of the samples were not valid. \n", ...
        methods(i), p_invalid_val(i));
end

% ----------------------------- Reports ------------------------------------
reports.ident.methods          = methods;
reports.ident.num_out          = num_out;
reports.ident.num_in           = num_in;
reports.ident.total            = num_all_ident;
reports.ident.p_contained      = p_contained_ident;
reports.ident.p_invalid        = p_invalid_ident;

reports.val.methods            = methods;
reports.val.num_out            = num_out_val;
reports.val.num_in             = num_in_val;
reports.val.total              = num_all_val;
reports.val.p_contained        = p_contained_val;
reports.val.p_invalid          = p_invalid_val;

reports.testSuite_val          = testSuite_val;   % for convenience
reports.params_true            = params_true;     % includes testSuite (ident)

end

% ----------------------------- Aux (optional) -----------------------------
function [sys, params] = aux_set_p(p, params, dyn)
% Kept verbatim from the example (handy for param sweeps).
    [sys, ~, ~, p_t] = loadDynamics(dyn, 'diag', p);
    c = p(length(p_t)+1:end);
    params.U = zonotope(c, params.U.generators);
end
