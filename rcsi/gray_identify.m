function configs = gray_identify(sys_cora, R0, U, C, pe)
    lookup = struct('sys',sys_cora,'dyn',C.shared.dyn,'R0',R0,'U',U, ...
        'n_m',C.shared.n_m,'n_s',C.shared.n_s,'n_k',C.shared.n_k, ...
        'n_m_val',C.shared.n_m_val,'n_s_val',C.shared.n_s_val,'n_k_val',C.shared.n_k_val, ...
        'methodsGray', C.gray.methodsGray, 'constraints', C.shared.cs_base.constraints);

    % Build testSuite options
    optTS = ts_options_from_pe(C, pe, sys_cora);   
    optTS.stateSet = R0;
    
    conf_opts = struct('options_reach', C.shared.options_reach, ...
                       'cs', C.shared.cs_base, ...
                       'testS', optTS);
    
    % run conformance as before
    R_tmp = conformance_gray(lookup, conf_opts); 
  
    % Build a minimal configs set (true + first gray method) to evaluate
    params_true = struct('R0',R0,'U',U,'tFinal', sys_cora.dt*(C.shared.n_k-1));
    params_true.testSuite = createTestSuite(sys_cora, params_true, ...
        C.shared.n_k, C.shared.n_m, C.shared.n_s, optTS);

    configs = cell(2,1);
    configs{1} = struct('sys',sys_cora,'params',rmfield(params_true,'testSuite'),'options',C.shared.options_reach,'name','true');
    %options = C.shared.options_reach; options.cs = C.shared.cs_base;
    %options.cs.set_p = @(p,params) aux_set_p_gray(p, params, C.shared.dyn);
    options.cs.set_p = @(p,params) aux_set_p_gray(p, params, C.shared.dyn, sys_cora);
    options.cs.p0 = 0.01*randn(sys_cora.nrOfInputs,1);
    [cfg_i.params, results] = conform(sys_cora, params_true, options, C.gray.methodsGray(1)); 
    configs{2} = struct('sys',results.sys,'params',cfg_i.params,'options',C.shared.options_reach,'name',C.gray.methodsGray(1));
end