function S = build_nonlinear_pair(dyn_ss, dyn_arx, ldtype)
% Return both flavors + their default sets
% Usage: S = build_nonlinear_pair("lipschitz2D","lipschitz2D_ARX","standard");

    if nargin < 3, ldtype = "standard"; end

    [sys_ss, R0_ss, U_ss]   = custom_loadDynamics(dyn_ss,  ldtype);
    [sys_ax, R0_ax, U_ax]   = custom_loadDynamics(dyn_arx, ldtype);

    S = struct('sys_ss', sys_ss, 'R0_ss', R0_ss, 'U_ss', U_ss, ...
               'sys_ax', sys_ax, 'R0_ax', R0_ax, 'U_ax', U_ax);
end
