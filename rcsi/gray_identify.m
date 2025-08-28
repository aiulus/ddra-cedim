function configs = gray_identify(sys_cora, R0, U, C, pe, varargin)
% GRAY_IDENTIFY  Build configs {true, gray} with optional noise override W.
%
% Usage:
%   configs = gray_identify(sys_cora, R0, U, C, pe)
%   configs = gray_identify(..., 'overrideW', W)       % pass process-noise zonotope
%   configs = gray_identify(..., 'options_testS', S)   % test-suite opts (fallback = ts_options_from_pe)

    % ---- parse name/value ----
    p = inputParser;
    addParameter(p,'overrideW',[]);
    addParameter(p,'options_testS',[]);
    addParameter(p,'externalTS_train',[]);
    addParameter(p,'externalTS_val',[]);
    parse(p,varargin{:});
    TS_train_ext = p.Results.externalTS_train;
    TS_val_ext   = p.Results.externalTS_val;
    W_override = p.Results.overrideW;
    optTS_in   = p.Results.options_testS;

    % ---- build/choose test-suite options (use same PE as DDRA) ----
    optTS_auto = ts_options_from_pe(C, pe, sys_cora);    
    optTS_auto.stateSet = R0;
    if isempty(optTS_in)
        optTS = optTS_auto;
    else
        optTS = optTS_in;
        if ~isfield(optTS,'stateSet'), optTS.stateSet = R0; end
    end

    % ---- optional wrapper hook ----
    try
        LM = C.lowmem; 
        ccflag = true; 
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

    % ---- params_true (kept for configs{1}) + identification suite ----
    params_true = struct('R0',R0,'U',U,'tFinal', sys_cora.dt*(C.shared.n_k-1));

    if ~isempty(TS_train_ext)
        % Use *your* sequences (locked with DDRA)
        params_true.testSuite = TS_train_ext;
    else
        % Fallback to CORA-generated suite
        params_true.testSuite = createTestSuite(sys_cora, params_true, ...
            C.shared.n_k, C.shared.n_m, C.shared.n_s, optTS);
    end
    

    % ---- CORA-style initial estimates (identity generators) ----
    c_R0 = zeros(sys_cora.nrOfDims,1);
    c_U  = zeros(sys_cora.nrOfInputs,1);
    params_id_init = params_true;
    params_id_init.R0 = zonotope([c_R0, eye(sys_cora.nrOfDims)]);
    params_id_init.U  = zonotope([c_U,  eye(sys_cora.nrOfInputs)]);

    if ~isempty(W_override) && hasGenerators(W_override)
        % augment the system to expose w as extra input channels
        % fold W into U so CORA sees one combined input set 
        U_aug = cartProd(U, W_override);
        params_true.U   = U_aug;
        params_id_init.U = zonotope([center(U_aug), eye(dim(U_aug))]);
    end
    
    % To include additional measurement noise V:
    % sys_cora = augment_u_with_v(sys_cora); U_aug = cartProd(U_aug_or_U, V);


    % ---- conformance options ----
    options = C.shared.options_reach;
    options.cs = C.shared.cs_base;

    if exist('aux_set_p_gray','file') == 2
        options.cs.set_p = @(p,params) aux_set_p_gray(p, params, C.shared.dyn, sys_cora);
        options.cs.p0    = 0.01*randn(sys_cora.nrOfInputs,1);
    end

    % ---- identify (first gray method) ----
    type = C.gray.methodsGray(1);
    [params_hat, results] = conform(sys_cora, params_id_init, options, type);

    % ---- pack configs ----
    configs = cell(2,1);
    configs{1} = struct('sys',sys_cora, ...
                        'params',rmfield(params_true,'testSuite'), ...
                        'options',C.shared.options_reach, ...
                        'name','true');

    configs{2} = struct('sys',results.sys, ...
                        'params',params_hat, ...
                        'options',C.shared.options_reach, ...
                        'name',type);
end
