function [sys_out, params] = aux_set_p_gray(p, params, dyn, sys_base)
    % Robust RCSI "set_p" for gray-box: optionally rebuild system from p_model,
    % and ALWAYS align cU length to the effective input dimension in params.U.
    %
    % Handles cases where U was augmented (e.g., cartProd(U, W)).

    % columnize
    p = p(:);

    % Decide if this dynamics exposes numeric model parameters
    [~, ~, ~, p_t] = custom_loadDynamics(dyn, "diag");

    if isnumeric(p_t) && ~isempty(p_t)
        % Split p = [p_model; cU]
        n_model = numel(p_t);
        assert(numel(p) >= n_model, 'aux_set_p_gray: |p| < |p_model|');
        p_model = p(1:n_model);
        cU      = p(n_model+1:end);

        % Rebuild system from model params
        [sys_out, ~, ~] = custom_loadDynamics(dyn, "diag", p_model);
    else
        % No numeric model params (e.g., k-MSD as currently wired): keep system
        sys_out = sys_base;
        cU      = p;   % all of p is the input-center
    end

    % Shift U center, keeping generators; be robust to augmented U
    % Use accessors to avoid deprecation warnings.
    if ~isa(params.U, 'zonotope')
        % Defensive: most code paths already have a zonotope here.
        params.U = zonotope(params.U);
    end
    G = generators(params.U);      % (nu_eff × eta)
    nu_eff = size(G, 1);           % effective input dimension seen by CORA now

    % Pad/truncate cU to match nu_eff
    if numel(cU) < nu_eff
        cU = [cU; zeros(nu_eff - numel(cU), 1)];
    elseif numel(cU) > nu_eff
        cU = cU(1:nu_eff);
    end

    % Rebuild U with same generators, updated center
    params.U = zonotope(cU, G);
end
