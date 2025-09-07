function Wd = normalizeWForGray(sys_, W_in)
    % Return [] for no disturbances, or a zonotope in R^{nrOfDisturbances}
    nw = sys_.nrOfDisturbances;
    if nw == 0 || isempty(W_in), Wd = []; return; end
    if ~isa(W_in,'zonotope')
        warning('overrideW is not a zonotope; dropping to zero-disturbance.');
        Wd = zonotope(zeros(nw,1)); return;
    end

    Ex = getfielddef(sys_, 'E', []);
    if isempty(Ex) || size(center(W_in),1) ~= size(Ex,1)
        % Fallback with correct dimension, but warn
        warning('coerceWToSys: incompatible dims or empty E; using zero Wd.');
        Wd = zonotope(zeros(nw,1)); return;
    end

    % --- conservative preimage mapping: E*Wd ⊇ W_in ---
    c  = center(W_in);                % nx×1
    Gx = generators(W_in);            % nx×g
    % min-norm preimages via lsqminnorm (more stable than pinv for scale)
    rc = lsqminnorm(Ex, c);           % nd×1
    R  = zeros(nw, size(Gx,2));       % nd×g
    res = zeros(1, size(Gx,2));
    for j = 1:size(Gx,2)
        rj = lsqminnorm(Ex, Gx(:,j)); % nd×1
        R(:,j) = rj;
        res(j) = norm(Ex*rj - Gx(:,j), inf);
    end
    % pick a tiny padding so that E*Wd ⊇ W_in even if some gj not exact
    eps_inf = max([res, 0]) + 1e-12;
    if eps_inf > 0
        pad = eps_inf * eye(nw);      % nd×nd
        Wd = zonotope([rc, [R, pad]]);
    else
        Wd = zonotope([rc, R]);
    end
end
