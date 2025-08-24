function Ublock = gen_input_block(dim_u, n_k, Uset, pe)
%GEN_INPUT_BLOCK  Build a dim_u x n_k input block honoring PE settings.
% pe.mode: 'randn' | 'sinWave' (default 'randn')
% pe.order: positive integer L (effective L is clamped to feasible range)
% pe.strength: scalar (default 1)
% pe.deterministic: logical (default false) — for reproducible bases

    mode  = getfielddef(pe,'mode','randn');
    Lreq  = max(1, round(getfielddef(pe,'order',1)));
    strength = getfielddef(pe,'strength',1);
    det   = getfielddef(pe,'deterministic',false);

    % Effective order cannot exceed horizon-1 or generator count
    Leff_time = max(1, min(Lreq, n_k-1));
    GU = [];
    if isa(Uset,'zonotope')
        GU = Uset.generators; % dim_u x eta_u (possibly empty)
    end
    eta_u = size(GU,2);
    Leff_gen  = max(1, min(Lreq, max(1,eta_u)));
    Leff = min(Leff_time, Leff_gen);

    if det, rng(1,'twister'); end

    % --- Build temporal bases S (Leff x n_k) ---
    switch lower(mode)
        case 'sinwave'
            t = 0:n_k-1;
            f = linspace(0.08, 0.45, Leff);
            phi = 2*pi*rand(1,Leff);
            S = zeros(Leff, n_k);
            for l=1:Leff
                S(l,:) = sin(2*pi*f(l)*t + phi(l));
            end
        otherwise  % 'randn' (default)
            % Low-rank temporal model: Leff independent rows
            S = randn(Leff, n_k);
    end

    % --- Map to input space ---
    cU = isa(Uset,'zonotope') * center(Uset) + ~isa(Uset,'zonotope') * zeros(dim_u,1);

    if ~isempty(GU)
        % Use only Leff generator channels -> temporal rank Leff reflected in coeffs
        coeff = zeros(eta_u, n_k);
        coeff(1:Leff,:) = strength * tanh(S);          % [-1,1] bounded coeffs
        Ublock = cU + GU * coeff;                       % dim_u x n_k
    else
        % No explicit generators — build directions in R^{dim_u}
        V = randn(dim_u, Leff); V = V ./ max(1e-12, vecnorm(V));  % orth-ish
        Ublock = cU + strength * (V * S);               % dim_u x n_k
    end
end
