function configs = gray_identify(sys_cora, R0, U, C, pe, varargin)
% GRAY_IDENTIFY  Build configs {true, gray} with optional noise override W.
%
% Usage:
%   configs = gray_identify(sys_cora, R0, U, C, pe)
%   configs = gray_identify(..., 'overrideW', W)         % process-noise zonotope
%   configs = gray_identify(..., 'options_testS', S)     % testSuite opts
%   configs = gray_identify(..., 'externalTS_train', TS) % CORA testCase objects

    % ---- parse name/value ----
    p = inputParser;
    addParameter(p,'overrideW',[]);
    addParameter(p,'options_testS',[]);
    addParameter(p,'externalTS_train',[]);
    addParameter(p,'externalTS_val',[]);  
    parse(p,varargin{:});
    TS_train_ext = p.Results.externalTS_train;
    W_override   = p.Results.overrideW;
    optTS_in     = p.Results.options_testS;

    % ---- build/choose test-suite options (use same PE as DDRA) ----
    optTS_auto = ts_options_from_pe(C, pe, sys_cora);
    optTS_auto.stateSet = R0;
    if isempty(optTS_in)
        optTS = optTS_auto;
    else
        optTS = optTS_in;
        if ~isfield(optTS,'stateSet'), optTS.stateSet = R0; end
    end

    % ---- optional wrapper hook (kept as-is) ----
    try
        LM      = C.lowmem;
        ccflag  = true; 
        if isfield(LM,'gray_check_contain'), ccflag = LM.gray_check_contain; end
        conf_opts = struct('options_reach', C.shared.options_reach, ...
                           'cs',            C.shared.cs_base, ...
                           'testS',         optTS, ...
                           'check_contain', ccflag);
        if ~isempty(W_override), conf_opts.W = W_override; end
        lookup = struct('sys',sys_cora,'dyn',C.shared.dyn,'R0',R0,'U',U, ...
                        'n_m',C.shared.n_m,'n_s',C.shared.n_s,'n_k',C.shared.n_k, ...
                        'n_m_val',C.shared.n_m_val,'n_s_val',C.shared.n_s_val,'n_k_val',C.shared.n_k_val, ...
                        'methodsGray', C.gray.methodsGray, ...
                        'constraints', C.shared.cs_base.constraints);
        conformance_gray(lookup, conf_opts);
    catch
        % wrapper absent / different API -> ignore
    end

    % ---- base params (true) ----
    params_true = struct('R0',R0,'U',U,'tFinal', sys_cora.dt*(C.shared.n_k-1));

    % ---- CORA-style initial estimates (identity generators) ----
    c_R0 = zeros(sys_cora.nrOfDims,1);
    c_U  = zeros(sys_cora.nrOfInputs,1);
    params_id_init = struct();
    params_id_init.tFinal = params_true.tFinal;
    params_id_init.R0     = zonotope([c_R0, eye(sys_cora.nrOfDims)]);
    params_id_init.U      = zonotope([c_U,  eye(sys_cora.nrOfInputs)]);

    % ---- NEW: set the identification testSuite explicitly ----
    if ~isempty(TS_train_ext)
        % Accept cell-of-testCase or single testCase
        if isa(TS_train_ext,'testCase')
            params_id_init.testSuite = {TS_train_ext};
        else
            assert(iscell(TS_train_ext), ...
                'externalTS_train must be a testCase or cell array of testCase.');
            % quick type check (first element)
            assert(~isempty(TS_train_ext) && isa(TS_train_ext{1},'testCase'), ...
                'externalTS_train must contain CORA testCase objects.');
            params_id_init.testSuite = TS_train_ext;
        end
        % For plotting/sanity “true” config doesn’t need the suite:
        params_true = rmfield_safe(params_true,'testSuite');
    else
        % Fallback to CORA-generated suite (same PE settings)
        params_id_init.testSuite = createTestSuite(sys_cora, ...
            struct('R0',R0,'U',U,'tFinal',params_true.tFinal), ...
            C.shared.n_k, C.shared.n_m, C.shared.n_s, optTS);
    end

    % ---- OPTIONAL: fold W_override into inputs for identification ----
    if ~isempty(W_override) && hasGenerators(W_override)
        % Augment inputs with disturbance channels for identification
        U_aug = cartProd(U, W_override);
        % true params (only for potential plotting)
        params_true.U   = U_aug;
        % identification params
        params_id_init.U = zonotope([center(U_aug), eye(dim(U_aug))]);
    end

    % ---- conformance options ----
    options      = C.shared.options_reach;
    options.cs   = C.shared.cs_base;
    if exist('aux_set_p_gray','file') == 2
        options.cs.set_p = @(p,params) aux_set_p_gray(p, params, C.shared.dyn, sys_cora);
        % minimal p0 (only U-center if no numeric plant params)
        options.cs.p0    = 0.01*randn(sys_cora.nrOfInputs,1);
    end

    % ---- identify (first gray method) ----
    type = C.gray.methodsGray(1);
    [params_hat, results] = conform(sys_cora, params_id_init, options, type);

    % ---- pack configs ----
    configs = cell(2,1);
    configs{1} = struct('sys',sys_cora, ...
                        'params',params_true, ...
                        'options',rmfield_safe(C.shared.options_reach,'cs'), ... % pure reach opts
                        'name','true');

    configs{2} = struct('sys',results.sys, ...
                        'params',params_hat, ...
                        'options',rmfield_safe(C.shared.options_reach,'cs'), ...
                        'name',type);
end

% --------- tiny helpers ----------
function S = rmfield_safe(S, f)
    if isstruct(S) && isfield(S,f), S = rmfield(S,f); end
end
function tf = hasGenerators(Z)
    try
        G = generators(Z);
        tf = ~isempty(G) && any(abs(G(:))>0);
    catch
        tf = false;
    end
end
