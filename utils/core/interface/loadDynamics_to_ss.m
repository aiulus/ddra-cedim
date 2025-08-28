function [sysd_ss, R0, U, p_true] = loadDynamics_to_ss(dynamics, type, p)
% Returns discrete ss equivalent when the CORA system is linearSysDT.

    if nargin==1, [sys_cora, R0, U, p_true] = loadDynamics(dynamics);
    elseif nargin==2, [sys_cora, R0, U, p_true] = loadDynamics(dynamics,type);
    else, [sys_cora, R0, U, p_true] = loadDynamics(dynamics,type,p);
    end
    sysd_ss = cora_to_matlab_ss(sys_cora);
end
