function [axes, defaults] = init_axes(cfg, grid)
    axes.NMS = pick(getfielddef(grid,'n_m_list',   cfg.shared.n_m));
    axes.NSS = pick(getfielddef(grid,'n_s_list',   cfg.shared.n_s));
    axes.NKS = pick(getfielddef(grid,'n_k_list',   cfg.shared.n_k));
    axes.DS  = pick(getfielddef(grid,'D_list',     getfielddef(cfg,'D', 3)));
    axes.AWS = pick(getfielddef(grid,'alpha_w_list', cfg.ddra.alpha_w));
    axes.PEs = pick(getfielddef(grid,'pe_list',    {struct('mode','randn')}));
    
    defaults = struct('seed',1);
end