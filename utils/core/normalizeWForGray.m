function Wd = normalizeWForGray(sys_, W_in, eps_inf)
% Map state-noise W_in to a disturbance-space set Wd s.t. E*Wd \supseteq W_in (conservative).
% If eps_inf is omitted, choose it adaptively from reconstruction residuals.

    % No disturbance channel or empty input --> return [] or zero-disturbance set
    if ~isa(W_in,'zonotope') || isempty(W_in)
        warning('normalizeWForGray: W_in must be a zonotope; returning zero set.');
        Wd = zonotope(zeros(max(1, sys_.nrOfDisturbances),1));
        return;
    end

    E = getfielddef(sys_, 'E', []);
    if isempty(E)
        % No disturbance channels in the model
        Wd = zonotope(zeros(max(1, sys_.nrOfDisturbances),1));
        return;
    end

    % Dimension checks
    nxW = size(center(W_in),1);
    [nE, nd] = size(E);
    if nxW ~= nE
        warning('normalizeWForGray: dim(W_in)=%d does not match rows(E)=%d; returning zero set.', nxW, nE);
        Wd = zonotope(zeros(nd,1));
        return;
    end

    % Preimages via min-norm LS
    c  = center(W_in);      % nE x 1
    Gx = generators(W_in);  % nE x g (allow g=0)

    rc = lsqminnorm(E, c);          % nd x 1
    R  = zeros(nd, size(Gx,2));
    res = zeros(1, size(Gx,2));
    for j = 1:size(Gx,2)
        rj     = lsqminnorm(E, Gx(:,j));  % nd x 1
        R(:,j) = rj;
        res(j) = norm(E*rj - Gx(:,j), inf);
    end

    % Pick/compute padding to guarantee E*Wd ⊇ W_in
    if ~exist('eps_inf','var') || isempty(eps_inf)
        % Residual-adaptive; fall back to tiny positive
        eps_inf = max([res, 0]) + 1e-12;
        % Optional: cap by gen magnitude to avoid huge padding on ill-conditioned E
        try
            gnorm = max(sum(abs(Gx),1));
            eps_inf = min(eps_inf, 1e-6 * max(1, gnorm));
        catch
        end
    end

    if eps_inf > 0
        pad = eps_inf * eye(nd);     % uniform inf-norm ball in disturbance space
        Wd  = zonotope([rc, [R, pad]]);
    else
        Wd  = zonotope([rc, R]);
    end
end

function v = getfielddef(S,f,d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end
