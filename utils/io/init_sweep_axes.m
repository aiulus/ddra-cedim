function [axes, baseC] = init_sweep_axes(cfg, grid)
    % Safe access to nested structs (avoid eager field access)
    shared = getfielddef(cfg, 'shared', struct());
    ddra   = getfielddef(cfg, 'ddra',   struct());

    axes.n_m     = pick(getfielddef(grid, 'n_m_list',     getfielddef(shared, 'n_m', 1)));
    axes.n_s     = pick(getfielddef(grid, 'n_s_list',     getfielddef(shared, 'n_s', 1)));
    axes.n_k     = pick(getfielddef(grid, 'n_k_list',     getfielddef(shared, 'n_k', 1)));
    axes.D       = pick(getfielddef(grid, 'D_list',       getfielddef(cfg,    'D',   3)));
    axes.alpha_w = pick(getfielddef(grid, 'alpha_w_list', getfielddef(ddra,   'alpha_w', 0)));
    axes.pe      = pick(getfielddef(grid, 'pe_list',      {struct('mode','randn')}));

    baseC = cfg;
    baseC.shared.testSuite_mode = "ddra_like";
end
