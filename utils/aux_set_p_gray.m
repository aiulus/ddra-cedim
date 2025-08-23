function [sys_out, params] = aux_set_p_gray(p, params, dyn, sys_base)
    % Probe if this dynamics has a numeric parameter vector
    [~, ~, ~, p_t] = custom_loadDynamics(dyn, "diag");

    if isnumeric(p_t) && ~isempty(p_t)
        % Numeric model params exist -> split [p_model; cU] and rebuild sys
        n_model = numel(p_t);
        p_model = p(1:n_model);
        cU      = p(n_model+1:end);
        [sys_out, ~, ~] = custom_loadDynamics(dyn, "diag", p_model);
    else
        % No numeric params (e.g., k-MSD): keep the current system unchanged
        sys_out = sys_base;
        cU      = p;  % all of p is the U-center
    end

    % Shift U center, keep generators
    params.U = zonotope(cU, params.U.generators);
end
