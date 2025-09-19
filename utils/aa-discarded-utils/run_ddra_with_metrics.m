function [stats, R_ddra] = run_ddra_with_metrics(sys_ss, R0_ss, U_ss, look, val_counts)
% val_counts: struct with n_m_val, n_s_val, n_k_val
% Returns stats struct with containment (%) and size proxy

    t0 = tic;
    [~, R_ddra] = ddra_nonlinear(sys_ss, look);
    t_total = toc(t0);

    % validation suite on the same state-space system
    params_val = struct('R0', R0_ss, 'U', U_ss, ...
                        'tFinal', sys_ss.dt * val_counts.n_k_val - sys_ss.dt);
    testSuite_val = createTestSuite(sys_ss, params_val, ...
                        val_counts.n_k_val, val_counts.n_m_val, val_counts.n_s_val);

    [contain_pct, sizeI] = eval_reach_against_samples(R_ddra, testSuite_val);

    stats = struct('contain_pct', contain_pct, ...
                   'sizeI', sizeI, ...
                   't_total', t_total, ...
                   't_learn', NaN, 't_check', NaN, 't_infer', NaN);
end
