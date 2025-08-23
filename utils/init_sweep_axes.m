function [axes, baseC] = init_sweep_axes(cfg, grid)
    axes.n_m   = pick(getfielddef(grid,'n_m_list',   cfg.shared.n_m));
    axes.n_s   = pick(getfielddef(grid,'n_s_list',   cfg.shared.n_s));
    axes.n_k   = pick(getfielddef(grid,'n_k_list',   cfg.shared.n_k));
    axes.D     = pick(getfielddef(grid,'D_list',     getfielddef(cfg,'D', 3)));
    axes.alpha_w = pick(getfielddef(grid,'alpha_w_list', cfg.ddra.alpha_w));
    axes.pe    = pick(getfielddef(grid,'pe_list',    {struct('mode','randn')}));
    baseC = cfg; baseC.shared.testSuite_mode = "ddra_like";
end