rng(1,'twister');

% --------- Config ----------
dyn_ss  = "lipschitz2D";       % or "poly2D"
dyn_arx = "lipschitz2D_ARX";   % or "poly2D_ARX"
ldtype  = "standard";

budgets = struct('n_s', 20, 'n_k', 6, 'n_k_val', 6, 'n_m_val', 3);
ddra_noise = struct('eta_w', 1, 'alpha_w', 0.01, 'c_w', [0;0]);

sweep_n_m = [2 4 8 16 32];

% Conformance options
conf_opts = struct();
conf_opts.options_reach = struct('zonotopeOrder',100,'tensorOrder',2,'errorOrder',1,'tensorOrderOutput',2,'verbose',false);
conf_opts.cs   = struct('robustnessMargin',1e-9,'verbose',false,'cost',"interval",'constraints',"half");
conf_opts.testS = struct('p_extr',0.3);
conf_opts.approx = struct('gp_parallel',false,'gp_pop_size',40,'gp_num_gen',20, ...
                          'gp_func_names', {{'times','plus','square'}}, ...
                          'gp_max_genes',2,'gp_max_depth',2, ...
                          'cgp_num_gen',4,'cgp_pop_size_base',5,'rng_seed',1);

% --------- Load systems once ----------
S = build_nonlinear_pair(dyn_ss, dyn_arx, ldtype);
val_counts = struct('n_m_val', budgets.n_m_val, 'n_s_val', budgets.n_s, 'n_k_val', budgets.n_k_val);

% --------- Sweep ----------
SUMMARY = table();
for nm = sweep_n_m
    budgets.n_m = nm;
    look = make_ddra_lookup(S.sys_ss, S.R0_ss, S.U_ss, budgets, ddra_noise);

    % DDRA
    [stats_ddra, R_ddra] = run_ddra_with_metrics(S.sys_ss, S.R0_ss, S.U_ss, look, val_counts); 

    % RCSI-black
    lookup_bb = struct('sys', struct('dyn', dyn_arx, 'type', ldtype), ...
                       'n_m', nm, 'n_s', budgets.n_s, 'n_k', budgets.n_k, ...
                       'n_m_train', 100, 'n_s_train', 10, 'n_k_train', budgets.n_k, ...
                       'n_m_val', budgets.n_m_val, 'n_s_val', budgets.n_s, 'n_k_val', budgets.n_k_val, ...
                       'methodsBlack', ["blackGP"], 'constraints', "half");
    [stats_bb, ~] = conformance_black(lookup_bb, conf_opts);

    SUMMARY = [SUMMARY; table(nm, ...
                stats_ddra.contain_pct, stats_bb.contain_pct, ...
                stats_ddra.sizeI,       stats_bb.sizeI, ...
                stats_ddra.t_total,     stats_bb.t_total, ...
                'VariableNames', {'n_m','cval_ddra','cval_black','sizeI_ddra','sizeI_black','t_ddra_total','t_black_total'})];
end

disp(SUMMARY);
