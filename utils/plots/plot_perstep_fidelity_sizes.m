function plot_perstep_fidelity_sizes(artifact_mat, varargin)
% Plots (top) per-k containment on VAL for Gray & DDRA, averaged over blocks,
% and (bottom) per-k set size (interval width sum) for True, DDRA, Gray.
% Uses only data inside an artifact saved by run_sweeps.

% Args (name/value):
%   'Dims'      : output dims to show in legend (used only for title), default [1 2]
%   'Reduce'    : zonotope order cap for internal reductions (Girard), default 60
%   'SaveBase'  : if nonempty, save PNG/PDF with this base path
%   'Show'      : true/false (default true)

p = inputParser;
addParameter(p,'Dims',[1 2]);
addParameter(p,'Reduce',60);
addParameter(p,'SaveBase',"");
addParameter(p,'Show',true);
parse(p,varargin{:});
prm = p.Results;

S = load(artifact_mat);
sysT = S.sys_ddra;          % "true" linearSysDT
sysG = S.sys_gray;          % identified Gray linearSysDT
VAL  = S.VAL;               % VAL.x0{b}, VAL.u{b}, VAL.y{b}
W_eff= S.W_eff;             % zonotope in state space
M_AB = S.M_AB;              % DDRA matrix zonotope (A,B)

% disturbance mapping (true & gray)
W_true = coerceWToSys(sysT, W_eff);
W_pred = [];
if isprop(sysG,'nrOfDisturbances') && sysG.nrOfDisturbances>0
    W_pred = normalizeWForGray(sysG, W_eff);
end

Kred = prm.Reduce;
nkv  = size(VAL.y{1},1);
Bv   = numel(VAL.x0);
ny   = size(sysT.C,1);
m    = size(sysT.D,2);

% accumulators
cov_gray_k = zeros(nkv,1); cov_ddra_k = zeros(nkv,1);
w_true_k   = zeros(nkv,1); w_gray_k   = zeros(nkv,1); w_ddra_k = zeros(nkv,1);
cnt_k      = zeros(nkv,1);

for b = 1:Bv
    % --- True & Gray reach along block b
    params_true = struct('R0', VAL.R0, 'u', VAL.u{b}', 'tFinal', sysT.dt*(nkv-1));
    if ~isempty(W_true), params_true.W = W_true; end
    Rt = reach(sysT, params_true, struct('zonotopeOrder',Kred,'reductionTechnique','girard')).timePoint.set;

    R0_gray = zonotope(zeros(size(center(VAL.R0))), generators(VAL.R0)) + VAL.x0{b};
    params_gray = struct('R0', R0_gray, 'u', VAL.u{b}', 'tFinal', sysG.dt*(nkv-1));
    if ~isempty(W_pred), params_gray.W = W_pred; end
    Rg = reach(sysG, params_gray, struct('zonotopeOrder',Kred,'reductionTechnique','girard')).timePoint.set;

    % --- DDRA pre-update walk (mirrors streaming)
    Xk = reduce(VAL.R0 + VAL.x0{b}, 'girard', Kred);
    for k = 1:nkv
        Xk   = reduce(Xk, 'girard', Kred);
        uk   = VAL.u{b}(k,:).';
        U_pt = zonotope(uk);

        Yt = sysT.C*Rt{k} + sysT.D*uk;
        Yg = sysG.C*Rg{k} + sysG.D*uk;
        Yd = sysT.C*Xk    + sysT.D*uk;                    % DDRA output at PRE-update k

        % widths (sum of interval widths across outputs)
        w_true_k(k) = w_true_k(k) + interval_width_sum(Yt);
        w_gray_k(k) = w_gray_k(k) + interval_width_sum(Yg);
        w_ddra_k(k) = w_ddra_k(k) + interval_width_sum(Yd);

        % containment - uses interval hulls
        yk = VAL.y{b}(k,:).';
        if contains_interval(yk, Yg, 1e-6),   cov_gray_k(k) = cov_gray_k(k) + 1; end
        if contains_interval(yk, Yd, 1e-6),   cov_ddra_k(k) = cov_ddra_k(k) + 1; end
        cnt_k(k) = cnt_k(k) + 1;

        % propagate DDRA
        Xk = M_AB * cartProd(Xk, U_pt) + W_eff;
    end
end

% average over blocks, convert coverage to %
w_true_k = w_true_k ./ max(1,cnt_k);
w_gray_k = w_gray_k ./ max(1,cnt_k);
w_ddra_k = w_ddra_k ./ max(1,cnt_k);
cov_gray_k = 100*(cov_gray_k ./ max(1,cnt_k));
cov_ddra_k = 100*(cov_ddra_k ./ max(1,cnt_k));

% ------------- plotting -------------
if ~prm.Show, return; end
f = figure('Color','w','Name','Per-step fidelity & sizes');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% Fidelity
nexttile; hold on; grid on;
plot(1:nkv, cov_gray_k,'-s','LineWidth',1.6,'DisplayName','Gray coverage (%)');
plot(1:nkv, cov_ddra_k,'-o','LineWidth',1.6,'DisplayName','DDRA coverage (%)');
xlabel('time step k'); ylabel('coverage (%)'); title('Per-step containment on VAL'); legend('Location','best');

% Sizes
nexttile; hold on; grid on;
plot(1:nkv, w_true_k,'-','LineWidth',1.6,'DisplayName','True set size');
plot(1:nkv, w_ddra_k,'-o','LineWidth',1.6,'DisplayName','DDRA set size');
plot(1:nkv, w_gray_k,'-s','LineWidth',1.6,'DisplayName','Gray set size');
xlabel('time step k'); ylabel('\Sigma interval widths across outputs');
ttl = sprintf('Per-step interval-size proxy (dims ~ [%s])', num2str(prm.Dims));
title(ttl); legend('Location','best');

% optional save
if strlength(prm.SaveBase)>0
    saveas(f, prm.SaveBase + "_perstep_sizes.png");
    try, exportgraphics(f, prm.SaveBase + "_perstep_sizes.pdf",'ContentType','vector'); catch, end
end
end

% ---------- helpers ----------
function w = interval_width_sum(S)
    Iv = interval(S);
    try, lo = infimum(Iv); hi = supremum(Iv);
    catch, lo = Iv.inf;    hi = Iv.sup;
    end
    w = sum(max(double(hi - lo), 0), 'all');
end

function tf = contains_interval(y, S, tol)
    Iv = interval(S);
    try, lo = infimum(Iv); hi = supremum(Iv);
    catch, lo = Iv.inf;    hi = Iv.sup;
    end
    tf = all(y <= hi + tol) && all(y >= lo - tol);
end
