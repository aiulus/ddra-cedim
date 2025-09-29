function [safe_true, safe_alg_tau] = eval_safety_labels( ...
    sysT, sysG, VAL, R0, M_AB, W_eff, W_pred, H, h, k_set, tau_vec, reach_ord)

% Inputs:
%   sysT, sysG  : linearSysDT (true + Gray)
%   VAL         : struct with x0{b}, u{b}, y{b}
%   R0          : zonotope initial-set (state space)  <-- NEW (explicit)
%   M_AB        : DDRA matZonotope
%   W_eff       : W in true state-space coords
%   W_pred      : W mapped to Gray disturbance space (or [])
%   H,h         : polyhedral spec Hy <= h in output space
%   k_set       : steps to enforce (default all)
%   tau_vec     : tightening sweep (default 0:0.01:0.2 or similar)
%   reach_ord   : Girard order for reach
%
% Outputs:
%   safe_true    : logical (B x 1) for 'per_block' OR (sum_B nkv x 1) for 'per_step'
%   safe_alg_tau : struct with .ddra and .gray logical decisions for each tau

if nargin < 11 || isempty(k_set), k_set = 1:size(VAL.y{1},1); end
if nargin < 12 || isempty(tau_vec), tau_vec = linspace(0,0.2,21); end
if nargin < 13 || isempty(reach_ord), reach_ord = 80; end

% Prepare reaches (once per block)
Bv  = numel(VAL.x0);
nkv = size(VAL.y{1},1);
ny  = size(sysT.C,1);

if size(H,2) ~= ny
    error('Safety spec H must have %d columns (ny). Got H: %dx%d. Fix cfg.metrics.safety.H/h.', ny, size(H,1), size(H,2));
end
if size(h,1) ~= size(H,1)
    error('Safety spec h must have %d rows to match H. Got h: %dx1.', size(H,1), size(h,1));
end

% Direction set for support computations = rows of H (each constraint)
Ddirs = H.';  % (ny x nIneq)

optReach = struct('zonotopeOrder', reach_ord, 'reductionTechnique', 'girard');

% Accumulators for per-sample decisions
dec_true   = [];     % true safety labels (column)
dec_ddra_0 = [];     % DDRA decisions at tau=0 (kept for debugging)
dec_ddra_T = [];     % DDRA decisions across tau (nTau columns)
dec_gray_T = [];     % Gray  decisions across tau

for b = 1:Bv
    % -- True and Gray reaches along block b
    params_true = struct('R0', R0, 'u', VAL.u{b}', 'tFinal', sysT.dt*(nkv-1));
    Rt = reach(sysT, params_true, optReach).timePoint.set;   % cell{nkv} (state sets)

    R0_gray = zonotope(zeros(size(center(R0))), generators(R0)) + VAL.x0{b};
    params_gray = struct('R0', R0_gray, 'u', VAL.u{b}', 'tFinal', sysG.dt*(nkv-1));
    if isprop(sysG,'nrOfDisturbances') && sysG.nrOfDisturbances>0 && ~isempty(W_pred)
        params_gray.W = W_pred;
    end
    Rg = reach(sysG, params_gray, optReach).timePoint.set;

    % -- DDRA PRE-update walk for outputs
    Kred = reach_ord;
    Xk = reduce(R0 + VAL.x0{b}, 'girard', Kred);
    Yd = cell(nkv,1);
    for k = 1:nkv
        Xk = reduce(Xk, 'girard', Kred);
        uk = VAL.u{b}(k,:).';
        Yd{k} = sysT.C*Xk + sysT.D*uk;    % PRE-update output set
        Xk = M_AB * cartProd(Xk, zonotope(uk)) + W_eff;
    end

    % -- Build per-step safety for this block
    true_step = false(nkv,1);
    ddra_step = false(nkv, numel(tau_vec));
    gray_step = false(nkv, numel(tau_vec));

    for k = 1:nkv
        uk = VAL.u{b}(k,:).';

        % True output set at k
        Yt = sysT.C*Rt{k} + sysT.D*uk;   Yt = toZono(Yt);
        % Gray output set at k
        Yg = sysG.C*Rg{k} + sysG.D*uk;   Yg = toZono(Yg);
        % DDRA output set at k (already in Yd)
        Yd_k = toZono(Yd{k});

        % True label (no tightening)
        true_step(k) = is_safe_zono(H, h, Yt, Ddirs);

        % Algorithmic decisions for each tau
        for it = 1:numel(tau_vec)
            ht = h - tau_vec(it);
            ddra_step(k,it) = is_safe_zono(H, ht, Yd_k, Ddirs);
            gray_step(k,it) = is_safe_zono(H, ht, Yg,   Ddirs);
        end
    end

    % Reduce to per-block if requested
    mode = getfielddef(getfielddef(evalin('caller','cfg.metrics','struct()'),'safety','struct()'),'mode','per_block');
    if strcmpi(mode,'per_step')
        use_idx = k_set(:);
        dec_true   = [dec_true;   true_step(use_idx)];
        dec_ddra_T = [dec_ddra_T; ddra_step(use_idx,:)];
        dec_gray_T = [dec_gray_T; gray_step(use_idx,:)];
    else
        use_idx = k_set(:);
        dec_true   = [dec_true;   all(true_step(use_idx))];
        dec_ddra_T = [dec_ddra_T; all(ddra_step(use_idx,:), 1)];
        dec_gray_T = [dec_gray_T; all(gray_step(use_idx,:), 1)];
    end
end

safe_true        = logical(dec_true(:));
safe_alg_tau.ddra = logical(dec_ddra_T);
safe_alg_tau.gray = logical(dec_gray_T);
end

% ---------- helpers ----------
function tf = is_safe_zono(H, h, Z, Ddirs)
% Z is a (toZono) zonotope in output space.
% Safety: Hy <= h for all y in Z  <=>  support_Z(H_i^T) <= h_i for all i.
c = center(zonotope(Z));
G = generators(zonotope(Z));
s = support_zono_vec(Ddirs, c, G);    % 1 x nIneq (directions are columns of Ddirs)
tf = all(s(:) <= h(:) + 1e-9);        % small slack
end

function v = getfielddef(S, name, default)
if ~isstruct(S) || ~isfield(S, name) || isempty(S.(name))
    v = default;
else
    v = S.(name);
end
end
