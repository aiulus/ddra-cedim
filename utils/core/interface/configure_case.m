function C = configure_case(cfg, D, alpha_w, n_m, n_s, n_k)
    C = cfg; 
    C.shared.n_m = n_m; C.shared.n_s = n_s; C.shared.n_k = n_k;
    C.shared.n_m_val = max(2, min(n_m, getfielddef(cfg.shared,'n_m_val',2)));
    C.shared.n_s_val = n_s; C.shared.n_k_val = n_k;
    C.ddra.alpha_w = alpha_w;
    C.shared.testSuite_mode = "ddra_like";
    C.shared.seed = 1; rng(C.shared.seed,'twister');
    if C.shared.dyn == "k-Mass-SD", C.shared.dyn_p = D; end
end
