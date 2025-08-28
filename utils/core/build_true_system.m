function [sys_cora, sys_ddra, R0, U] = build_true_system(C)
    if isfield(C.shared,'dyn_p')
        [sys_cora, R0, U] = custom_loadDynamics(C.shared.dyn, C.shared.type, C.shared.dyn_p);
    else
        [sys_cora, R0, U] = custom_loadDynamics(C.shared.dyn, C.shared.type);
    end
    sys_ddra = cora_to_matlab_ss(sys_cora);
end