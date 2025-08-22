function R = conformance_gray(lookup, conf_opts)
% CONFORMANCE_GRAY  Gray-box reachset-conformance wrapper.
%   R = conformance_gray(lookup, conf_opts)
%
% Inputs 
%   lookup: struct with fields
%       .sys.dyn   (string) dynamics key for custom_loadDynamics()
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
%   - Implements GraySim & GraySeq/GrayLS as in the paper’s Sec. V (Alg. 2 & 3). :contentReference[oaicite:1]{index=1}
%   - Falls back to white-box if p_true is empty and you only want to size U’s center.
%
% Author:  Aybüke Ulusarslan (22-Aug-2025)

% ------------------------------ BEGIN CODE -------------------------------
rng(1,'twister'); % keep seed aligned with your pipeline

% Unpack lookup
dyn   = lookup.sys.dyn;
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

% Base options
options_reach = conf_opts.options_reach;
options_testS  = conf_opts.testS;

% Load dynamics + default sets
if dyn == "platoon"
    [sys, params_true.R0, params_true.U] = load_platoon(n_n, max(n_k, n_k_val));
else
    [sys, params_true.R0, params_true.U, p_true] = custom_loadDynamics(dyn, ldtyp);
end
if ~exist('p_true','var'); p_true = []; end
params_true.tFinal = sys.dt * n_k - sys.dt;

if sys.nrOfOutputs < 2
        constraints = {'gen'};
end

% Build identification test suite
params_true.testSuite = createTestSuite(sys, params_true, n_k, n_m, n_s, options_testS);

% Identification options (gray-box)
options = options_reach;
options.cs = conf_opts.cs;
options.cs.constraints = constraints;

% p0 stacks [model params; U-center] like the example
dim_p = numel(p_true);
p0    = 0.01*randn(dim_p + sys.nrOfInputs, 1);
options.cs.p0 = p0;
options.cs.set_p = @(p,params) aux_set_p_gray(p, params, dyn);

% “Configs” container as in CORA examples
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

% ---- Validation & Visualization ----
% Identification containment sanity-check on training data
num_out = 0; 
check_contain = 0; % Takes forever if set to 1
methodsList = ["true", methodsGray];
for m=1:length(params_true.testSuite)
    [~, eval_id] = validateReach(params_true.testSuite{m}, configs, check_contain);
    num_out = num_out + eval_id.num_out;
end
num_all = length(params_true.testSuite)*n_k*size(params_true.testSuite{1}.y,3);
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

num_out = 0; check_contain = 1;
for m=1:length(testSuite_val)
    [R, eval_val] = validateReach(testSuite_val{m}, configs, check_contain, plot_settings);
    num_out = num_out + eval_val.num_out;
end
num_all = length(testSuite_val)*n_k_val*size(testSuite_val{1}.y,3);
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

function [sys_out, params] = aux_set_p_gray(p, params, dyn)
    % Rebuild the system with the first part of p as model parameters
    % (CORA's example uses 'diag' here as well)
    [sys_out, ~, ~, p_t] = custom_loadDynamics(dyn, "diag", p);

    % Remaining entries of p set the center of U
    cU = p(length(p_t) + 1 : end);

    % Keep U's generators, shift only its center
    params.U = zonotope(cU, params.U.generators);
end

% ------------------------------- END CODE --------------------------------
