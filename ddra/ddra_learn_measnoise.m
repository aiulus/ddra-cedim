function [ABc, AV_one, V_meas, info] = ddra_learn_measnoise(Xminus, Uminus, Xplus, W, sys_ddra, C)
%DDRA_LEARN_MEASNOISE  Algorithm 4 (no AV assumption) for measured states.
% Learns a *deterministic* AB center plus a single zonotope AV_one that
% bounds residual effects (A*V + modeling mismatch), robust to state meas noise.
%
% Inputs:
%   Xminus, Uminus, Xplus : (nx/nu × N) TRAIN triples (can be *measured* states)
%   W : zonotope in *state space* (baseline, e.g., from B*U; may be inflated later)
%   sys_ddra : struct with fields A,B (used only for dims)
%   C : row config (unused; kept for signature parity)
%
% Outputs:
%   ABc    : numeric [nx × (nx+nu)] least-squares center
%   AV_one : zonotope (nx-dim) bounding residuals (per Alg. 4) + (-W) + (-V)
%   V_meas : zonotope used for state measurement noise (may be zero-zono)
%   info   : struct with fields N, nx, nu, mins, maxs for sanity/debug.

    nx = size(sys_ddra.A,1);
    nu = size(sys_ddra.B,2);

    N  = size(Xminus,2);
    Z  = [Xminus; Uminus];                 % (nx+nu) × N
    X1 = Xplus;                            % (nx × N)

    % --- try to recover V from config; otherwise default to zero ---
    V_meas = zonotope(zeros(nx,1));
    try
        if isfield(C,'data') && isfield(C.data,'train') && isfield(C.data.train,'meas') ...
                && isfield(C.data.train.meas,'enable') && C.data.train.meas.enable
            Vcfg = C.data.train.meas;
            if isfield(Vcfg,'V') && ~isempty(Vcfg.V)
                V_meas = Vcfg.V;   % already a zonotope
            elseif isfield(Vcfg,'gen_scale')
                V_meas = zonotope(zeros(nx,1), Vcfg.gen_scale*eye(nx));
            end
        end
    catch
    end

    % --- centers of W and V (for the center-only regressor) ---
    cW = center(W);  if isempty(cW), cW = zeros(nx,1); end
    cV = center(V_meas); if isempty(cV), cV = zeros(nx,1); end

    % --- (1) Center AB via LS on center-subtracted X1 ---
    %     ABc ≈ (X1 - cW*1^T - cV*1^T) * pinv(Z)
    C1 = X1 - cW * ones(1,N) - cV * ones(1,N);
    ABc = C1 * pinv(Z);     % [nx × (nx+nu)]

    % --- (2) Residuals: R = X1 - ABc*Z   (shape nx×N)
    R = X1 - ABc * Z;

    % --- (3) Per-dimension intervals for residuals
    %     AV_oneterm bounds the *unknown* AV term; include -W and -V to get soundness
    mins = min(R,[],2);
    maxs = max(R,[],2);
    Rbox = zonotope(interval(mins, maxs));

    AV_one = Rbox + (-1)*W + (-1)*V_meas;

    info = struct('N',N,'nx',nx,'nu',nu,'mins',mins,'maxs',maxs);
end
