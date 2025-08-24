function Tvalidate_g = gray_time_validate(configs, sys_cora, R0, U, C, pe)
    % Build validation test suite (same as gray_containment)
    optTS = ts_options_from_pe(C, pe);
    optTS.stateSet = R0;                         
    params_true = struct('R0',R0,'U',U,'tFinal', sys_cora.dt*C.shared.n_k_val - sys_cora.dt);
    testSuite_val = createTestSuite(sys_cora, params_true, ...
                       C.shared.n_k_val, C.shared.n_m_val, C.shared.n_s_val, optTS);

    % Time reachable-set computation without containment / plotting
    check_contain = 0;
    t0 = tic;
    for m=1:length(testSuite_val)
        validateReach(testSuite_val{m}, configs, check_contain);  
    end
    Tvalidate_g = toc(t0);
end
