function [R_true, R_data, out] = ddra_nonlinear(sys, lookup)
% ddra_nonlinear  Data-Driven Reachability for nonlinear DT systems
% Implements Algorithm 6 (Lipschitz Reachability) from Alanwar et al.
% Returns:
%   R_data : reachSet (CORA) with timePoint.set{k} = R'_k
%   R_true : (optional) model-based reach set; left empty unless requested
%   out    : struct with artifacts (data matrices, bounds, etc.)

% ---------- 0) Unpack  ----------
rand('seed',1);
dim_x = lookup.dim_x; dim_u = lookup.dim_u;

n_s = lookup.n_s; n_m = lookup.n_m; n_k = lookup.n_k; n_k_val = lookup.n_k_val;
totalsamples = n_m * n_s * n_k;

% X0
Gx = diag(ones(dim_x, 1));                    % diagonal template
if isfield(lookup,'eta_x') && lookup.eta_x ~= dim_x
    % optional: keep user-requested number of generators if provided
    Gx = zeros(dim_x, max(1,lookup.eta_x));   % fall back to zeros if eta_x < 1
    for i=1:min(dim_x, size(Gx,2)), Gx(i,i) = 1; end
end
X0 = zonotope(lookup.c_x + lookup.c_delta_x, lookup.alpha_x * Gx);

% U
Gu = diag(ones(dim_u, 1));                    % diagonal template
if isfield(lookup,'eta_u') && lookup.eta_u ~= dim_u
    Gu = zeros(dim_u, max(1,lookup.eta_u));
    for i=1:min(dim_u, size(Gu,2)), Gu(i,i) = 1; end
end
U  = zonotope(lookup.c_u + lookup.c_delta_u, lookup.alpha_u * Gu);

% W and its matrix zonotope Mw
Gw = ones(dim_x, 1);
if isfield(lookup,'eta_w') && lookup.eta_w > 1
    Gw = ones(dim_x, lookup.eta_w);
end
W  = zonotope(lookup.c_w, lookup.alpha_w * Gw);
Mw = aux_buildMatrixZonotope(W, dim_x, totalsamples);

% dynamics handle (x_{k+1} = f(x_k,u_k))
if isfield(lookup,'fun') && ~isempty(lookup.fun)
    fun = lookup.fun;
else
    try, fun = sys.mFile; catch, error('Provide lookup.fun: x_next = fun(x,u).'); end
end

% ---------- 1) Generate dataset D = (U_-, X) ----------
u_all = zeros(dim_u, totalsamples);
for i=1:totalsamples, u_all(:,i) = randPoint(U); end

n_blocks = n_m * n_s;
x_all  = zeros(dim_x*n_blocks, n_k+1);
x_free = zeros(dim_x*n_blocks, n_k+1);
for b=1:n_blocks
    r = (b-1)*dim_x + 1;
    x0b = randPoint(X0);
    x_all(r:r+dim_x-1,1)  = x0b;
    x_free(r:r+dim_x-1,1) = x0b;
end

idx = 1;
for b=1:n_blocks
    rx = (b-1)*dim_x + 1;
    for t=1:n_k
        u_bt   = u_all(:,idx);
        x_curr = x_free(rx:rx+dim_x-1,t);
        x_nf   = fun(x_curr, u_bt);
        x_free(rx:rx+dim_x-1,t+1) = x_nf;
        x_all (rx:rx+dim_x-1,t+1) = x_nf + randPoint(W);
        idx = idx + 1;
    end
end

% Flatten blocks to (dim × totalsamples)
X_minus = zeros(dim_x, totalsamples);
X_plus  = zeros(dim_x, totalsamples);
U_minus = zeros(dim_u, totalsamples);
col = 1;
for b=1:n_blocks
    rx = (b-1)*dim_x + 1;
    for t=1:n_k
        X_minus(:,col) = x_all(rx:rx+dim_x-1, t);
        X_plus(:, col) = x_all(rx:rx+dim_x-1, t+1);
        U_minus(:,col) = u_all(:, col);
        col = col + 1;
    end
end

% ---------- 2) Lipschitz inflation Z_eps (L and gamma via helper) ----------
stepsLip      = getf(lookup,'stepsLip',1);
initpointsLip = getf(lookup,'initpointsLip',50);
normType      = getf(lookup,'normType',2);
addZeps       = getf(lookup,'addZeps',true);

[gamma,L] = aux_compLipConst(fun, U, X0, stepsLip, initpointsLip, dim_x, normType);
eps_vec   = (L(:) .* gamma(:)) / 2;               % per-coord L*delta/2
Z_eps     = zonotope([zeros(dim_x,1), diag(eps_vec)]);

% ---------- 3) Algorithm 6 propagation ----------
R_cells = cell(n_k_val,1);
R_cells{1} = X0;

oneRow = ones(1, totalsamples);
C_Mw   = center(Mw);   % safer than Mw.center

for k = 1:n_k_val-1
    % linearization point
    x_star = center(R_cells{k});
    u_star = center(U);

    % local LS model around (x*,u*): M' * [1; x-x*; u-u*]
    Phi = [ oneRow;
            X_minus - x_star .* ones(dim_x, totalsamples);
            U_minus - u_star .* ones(dim_u, totalsamples) ];
    Mprime = (X_plus - C_Mw) * pinv(Phi);   % dim_x × (1+dim_x+dim_u)

    % residuals across data columns → interval enclosure (NO Minkowski difference)
    Res = X_plus - Mprime * Phi;            % each col is r_j
    lb = min(Res, [], 2);
    ub = max(Res, [], 2);
    V_one = zonotope(interval(lb, ub));     % residual enclosure (conservative)

    % Affine propagation: c + A*(R_k - x*) + B*(U - u*) + W + V_one (+ Z_eps)
    c = Mprime(:,1);
    A = Mprime(:, 2:1+dim_x);
    B = Mprime(:, 2+dim_x:end);

    RK_shift = R_cells{k} + (-x_star);
    U_shift  = U + (-u_star);

    R_next = c + A*RK_shift + B*U_shift + W + V_one;
    if addZeps
        R_next = R_next + Z_eps;
    end

    % Light reduction to control generator growth
    R_cells{k+1} = reduce(R_next, 'girard', 200);
end

% wrap as CORA reachSet
tp.set  = R_cells(2:end);
tp.time = num2cell((sys.dt:sys.dt:sys.dt*(n_k_val-1))');
R_data  = reachSet(tp);

% ---------- 4) Model-based (optional) ----------
R_true = [];
if isfield(lookup,'compute_model_reach') && lookup.compute_model_reach
    try
        params.R0 = X0; params.U = U; params.tFinal = sys.dt * n_k_val - sys.dt;
        opt.dim_x = dim_x;
        opt.zonotopeOrder = getf(lookup,{'reach','zonotopeOrder'},100);
        opt.tensorOrder   = getf(lookup,{'reach','tensorOrder'},2);
        opt.errorOrder    = getf(lookup,{'reach','errorOrder'},5);
        opt.W = W; opt.tStart = 0; opt.uTrans = 0;
        [R_true, ~] = reach_DT(sys, params, opt);
    catch
        % leave R_true empty on failure
    end
end

% ---------- 5) (Optional) plots ----------
if isfield(lookup,'plot_settings') && isfield(lookup.plot_settings,'projectedDims')
    aux_visualize( ...
        X0, R_true, R_data, ...
        lookup.plot_settings.projectedDims, ...
        getf(lookup.plot_settings,'numberofplots',5) );
end

% ---------- 6) Artifacts ----------
out = struct();
out.X_0T = X_minus; out.X_1T = X_plus; out.U_full = U_minus;
out.W = W; out.Mw = Mw; out.eps = eps_vec; out.L = L(:); out.gamma = gamma(:);
out.Mprime_last = exist('Mprime','var')*Mprime;
end
