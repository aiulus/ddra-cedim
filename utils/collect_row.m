function row = collect_row(C, D, alpha_w, n_m, n_s, n_k, pe, ddra_res, gray_res, T_ddra, T_gray)
    row.dyn      = string(C.shared.dyn);
    row.D        = D;
    row.n_m      = n_m; row.n_s = n_s; row.n_k = n_k;
    row.alpha_w  = alpha_w;
    row.pe_mode  = string(pe.mode);
    row.pe_order = getfielddef(pe,'order',NaN);
    
    % Fidelity
    row.ctrain_gray = gray_res.contain_train;
    row.cval_gray   = gray_res.contain_val;
    row.cval_ddra   = ddra_res.contain_val;
    
    % Conservatism
    row.sizeI_ddra  = ddra_res.size_interval;
    row.sizeI_gray  = gray_res.size_interval;
    
    % Conditioning
    row.rankZ_ddra  = ddra_res.rankZ;
    row.condZ_ddra  = ddra_res.condZ;
    
    % Runtimes
    row.t_ddra_learn   = T_ddra.learn;
    row.t_ddra_check   = T_ddra.check;
    row.t_ddra_infer   = T_ddra.infer;
    row.t_gray_learn   = T_gray.learn;
    row.t_gray_val     = T_gray.validate;
    row.t_gray_infer   = T_gray.infer;
end