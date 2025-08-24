function look = make_ddra_lookup(sys_ss, R0_ss, U_ss, budgets, ddra_noise)
% budgets: struct with fields n_m, n_s, n_k, n_k_val
% ddra_noise: struct with fields alpha_w, eta_w, c_w (vector)

    dim_x = sys_ss.nrOfDims; dim_u = sys_ss.nrOfInputs;

    look = struct();
    look.dim_x = dim_x; look.dim_u = dim_u;
    look.n_m   = budgets.n_m;  look.n_s = budgets.n_s;
    look.n_k   = budgets.n_k;  look.n_k_val = budgets.n_k_val;

    look.c_x = center(R0_ss);  look.c_delta_x = zeros(dim_x,1);
    look.alpha_x = 1;          look.eta_x = dim_x;

    look.c_u = center(U_ss);   look.c_delta_u = zeros(dim_u,1);
    look.alpha_u = 1;          look.eta_u = dim_u;

    look.c_w = ddra_noise.c_w; look.alpha_w = ddra_noise.alpha_w; look.eta_w = ddra_noise.eta_w;

    look.stepsLip = 1; look.initpointsLip = 50; look.addZeps = true;

    % dynamics handle for ddra_nonlinear
    look.fun = sys_ss.mFile;
end
