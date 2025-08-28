function row = pack_row(C, D, alpha_w, pe, ctrain_g, cval_g, cval_d, sizeI_d, sizeI_g, rankZ, condZ, tL, tC, tI, tLg, tVg, tIg)
    row = struct();
    row.dyn      = string(C.shared.dyn);
    row.D        = D;
    row.n_m      = C.shared.n_m; row.n_s = C.shared.n_s; row.n_k = C.shared.n_k;
    row.alpha_w  = alpha_w; row.pe_mode = string(pe.mode); row.pe_order = getfielddef(pe,'order',NaN);
    row.ctrain_gray = ctrain_g; row.cval_gray = cval_g; row.cval_ddra = cval_d;
    row.sizeI_ddra = sizeI_d; row.sizeI_gray = sizeI_g;
    row.rankZ_ddra = rankZ; row.condZ_ddra = condZ;
    row.t_ddra_learn=tL; row.t_ddra_check=tC; row.t_ddra_infer=tI;
    row.t_gray_learn=tLg; row.t_gray_val=tVg; row.t_gray_infer=tIg;
end