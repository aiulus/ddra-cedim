function plot_reach_all_onepanel(sys_true, sys_gray, VAL, varargin)
% PLOT_REACH_ALL_ONEPANEL
%   Single figure with: nominal traj, TRUE reach (R0/U), DDRA reach (M_AB,W),
%   and RCSI/Gray reach. All in output space, same axes, same inputs/horizon.
%
% Required:
%   sys_true : CORA linearSysDT (the “true” model used for data)
%   sys_gray : CORA linearSysDT (identified gray/RCSI model)
%   VAL      : struct with fields x0{b}, u{b}, optional y{b}, R0 (zonotope), U (zonotope)
%
% Name–Value (optional):
%   'Block'     : which block b in VAL to plot (default 1)
%   'Dims'      : [i j] output dims (default [1 2])
%   'MAB'       : matZonotope for DDRA recursion (default [])
%   'W'         : process-noise zonotope in state space (default [])
%   'Order'     : zonotope reduction order (default 100)
%   'Reduce'    : reduction technique (default 'girard')
%   'ShowSamples': true/false to overlay VAL.y (default true)
%   'Save'      : filename (no save if empty; default '')
%
% Notes:
%   • TRUE/GRAY reach use exact R0/U if provided (VAL.R0, VAL.U); otherwise nominal-only.
%   • DDRA applies X_{k+1} \in M_AB * [X_k; (u_k \minksum U)] \minksum W, then maps to output.
%   • All plotting is in output space (C,D).

    % ----- parse options -----
    opt = struct('Block',1,'Dims',[1 2],'MAB',[],'W',[], ...
                 'Order',100,'Reduce','girard','ShowSamples',true,'Save','');
    for i=1:2:numel(varargin), if isfield(opt,varargin{i}), opt.(varargin{i}) = varargin{i+1}; end, end
    b = opt.Block; dims = opt.Dims(:)'; assert(numel(dims)==2,'Dims must be 2-vector.');

    % ----- sanity checks & shapes -----
    assert(isa(sys_true,'linearSysDT') && isa(sys_gray,'linearSysDT'), 'Use CORA linearSysDT systems.');
    assert(abs(sys_true.dt - sys_gray.dt) < 1e-12, 'dt mismatch between true and gray systems.');
    nx = size(sys_true.A,1); ny = size(sys_true.C,1); mu = size(sys_true.B,2);
    assert(all(dims>=1 & dims<=ny), 'Dims out of range for outputs.');

    x0 = VAL.x0{b};                  % (nx × 1)
    u  = VAL.u{b};                   % (m × n_k) or (n_k × m)
    if size(u,2) < size(u,1), u = u.'; end     % make (n_k × m) → then transpose for CORA
    if size(u,1) ~= mu && size(u,2) == mu, u = u.'; end
    assert(size(u,1)==mu, 'VAL.u{%d} must be m×n_k.', b);
    n_k = size(u,2);

    % Uncertainty sets
    hasR0 = isfield(VAL,'R0') && isa(VAL.R0,'zonotope') && ~isempty(VAL.R0);
    hasU  = isfield(VAL,'U')  && isa(VAL.U, 'zonotope') && ~isempty(VAL.U);
    if ~(hasR0 && hasU)
        warning('plot_reach_all_onepanel:NoSets', ...
                'No VAL.R0/VAL.U provided — plotting nominal-only reach (no uncertainty).');
    end
    R0_set = iff(hasR0, @() VAL.R0, @() zonotope(x0));      % will offset by x0 below anyway
    U_set  = iff(hasU,  @() VAL.U,  @() zonotope(zeros(mu,1)));

    % Options (same for all)
    options = struct('zonotopeOrder', opt.Order, 'verbose', false);

    % ----- TRUE REACH -----
    paramsT = struct();
    paramsT.tStart = 0; paramsT.tFinal = sys_true.dt * (n_k-1);
    paramsT.R0     = R0_set + x0;     % exact initial uncertainty + nominal x0
    paramsT.U      = U_set;           % exact input uncertainty
    paramsT.u      = u.';             % CORA wants (n_k × m)
    Rtrue = reach(sys_true, paramsT, options);
    Xtrue = Rtrue.timePoint.set;
    % map to outputs: y = Cx + D u_k
    Ytrue = map_X_to_Y(Xtrue, sys_true.C, sys_true.D, u);

    % ----- GRAY/RCSI REACH (same R0/U/u) -----
    paramsG = paramsT;
    Rgray = reach(sys_gray, paramsG, options);
    Xgray = Rgray.timePoint.set;
    Ygray = map_X_to_Y(Xgray, sys_gray.C, sys_gray.D, u);

    % ----- DDRA REACH via M_AB (optional) -----
    Yddra = cell(1,n_k);
    if ~isempty(opt.MAB)
        Xdd = cell(1,n_k); Xdd{1} = R0_set + x0;
        for k = 1:n_k-1
            if hasU
                Uk = (U_set + u(:,k));
            else 
                Uk = zonotope(u(:,k));                         % u_k \minksum U  (zonotope with center u_k)
            end
            Zk = cartProd(Xdd{k}, Uk);                         % [x; u]
            Xk1 = opt.MAB * Zk;                                % matrix zonotope multiplication
            if ~isempty(opt.W), Xk1 = Xk1 + opt.W; end
            Xdd{k+1} = reduce(Xk1, opt.Reduce, opt.Order);
        end
        Yddra = map_X_to_Y(Xdd, sys_true.C, sys_true.D, u);    % map with true C,D (output space)
    end

    % ----- nominal trajectory (true model) -----
    xt = zeros(nx, n_k); xt(:,1) = x0;
    for k = 1:n_k-1, xt(:,k+1) = sys_true.A*xt(:,k) + sys_true.B*u(:,k); end
    yt = sys_true.C*xt + sys_true.D*u;

    % (optional) Samples
    haveY = isfield(VAL,'y') && numel(VAL.y)>=b && ~isempty(VAL.y{b});
    if haveY
        Ymeas = VAL.y{b};              % (n_k × ny) or (ny × n_k)
        if size(Ymeas,2)==ny && size(Ymeas,1)==n_k
            % ok
        elseif size(Ymeas,1)==ny && size(Ymeas,2)==n_k
            Ymeas = Ymeas.';           % make (n_k × ny)
        else
            haveY = false;              % ignore weird shapes
        end
    end

    % ----- Plot -----
    Cc = struct('true',[0.20 0.20 0.80], 'ddra',[0.85 0.20 0.20], 'gray',[0.20 0.55 0.30], ...
                'nom',[0 0 0], 'meas',[0.10 0.10 0.10]);
    figure('Color','w'); hold on; box on; grid on;
    title('Reachable output sets (true vs DDRA vs RCSI)'); xlabel(sprintf('y_{%d}',dims(1))); ylabel(sprintf('y_{%d}',dims(2)));

    % True
    [hTrue,~] = plot_series(Ytrue, dims, Cc.true, 0.10);
    set(hTrue(1),'DisplayName','True reach');

    % DDRA (if available)
    if ~isempty(opt.MAB)
        [hD,~] = plot_series(Yddra, dims, Cc.ddra, 0.12);
        set(hD(1),'DisplayName','DDRA reach');
    end

    % Gray
    [hG,~] = plot_series(Ygray, dims, Cc.gray, 0.10);
    set(hG(1),'DisplayName','RCSI/Gray reach');

    % Nominal
    hNom = plot(yt(dims(1),:), yt(dims(2),:), 'k-', 'Color', Cc.nom, 'LineWidth', 1.4, 'DisplayName','Nominal (true)');

    % Samples
    if haveY && opt.ShowSamples
        plot(Ymeas(:,dims(1)), Ymeas(:,dims(2)), 'o', 'Color', Cc.meas, 'MarkerSize', 3, ...
            'MarkerFaceColor', [0.6 0.6 0.6], 'HandleVisibility','off');
    end

    % shared limits
    lims = auto_limits({Ytrue, Ygray, iff(~isempty(opt.MAB), @() Yddra, @() {})}, dims, iff(haveY,@()Ymeas, @()[]));
    axis(lims);

    % legend
    L = {'Nominal (true)','True reach','DDRA reach','RCSI/Gray reach'};
    H = [hNom, hTrue(1), iff(~isempty(opt.MAB), @() hD(1), @() gobjects(0)), hG(1)];
    legend(H, L(1:numel(H)), 'Location','best');

    set(gca,'FontSize',12);

    if ~isempty(opt.Save)
        [p,~,~] = fileparts(opt.Save); if ~isempty(p) && ~exist(p,'dir'), mkdir(p); end
        exportgraphics(gcf, [opt.Save '.png'], 'Resolution', 200);
        exportgraphics(gcf, [opt.Save '.pdf'], 'ContentType','vector');
    end
end

% ---------- helpers ----------
function Y = map_X_to_Y(X, C, D, u)
    if isempty(D), D = zeros(size(C,1), size(u,1)); end
    K = numel(X); Y = cell(1,K);
    for k=1:K
        Xk = zonotope(X{k});
        Y{k} = zonotope(C*center(Xk) + D*u(:,k), C*generators(Xk));
    end
end

function [hs, lastIv] = plot_series(Ycells, dims, color, alpha)
    hs = gobjects(0); lastIv = [];
    for k=1:numel(Ycells)
        S = Ycells{k}; if isempty(S) || representsa(S,'emptySet'), continue; end
        S = reduce(S,'girard',min(1000,size(generators(S),2)));  %# keep plottable
        Iv = interval(project(S,dims)); lastIv = Iv;
        w = Iv.sup - Iv.inf; tol = 1e-12;
        if all(w <= tol)                 % point
            h = plot(S, dims, 'Color', color, 'Marker','.', 'LineStyle','none');
        elseif any(w <= tol)             % segment / degenerate
            h = plot(S, dims, 'Color', color, 'LineWidth', 0.9);
        else                              % 2D filled
            h = plot(S, dims, 'FaceColor', color, 'EdgeColor', color, 'FaceAlpha', alpha, 'LineWidth', 0.6);
        end
        if numel(hs)==0, hs = h; else, set(h,'HandleVisibility','off'); end
    end
end

function lims = auto_limits(groups, dims, Ymeas)
    minx = inf; maxx = -inf; miny = inf; maxy = -inf;
    for gi = 1:numel(groups)
        G = groups{gi}; if isempty(G), continue; end
        for k=1:numel(G)
            S = G{k}; if isempty(S), continue; end
            Iv = interval(project(S,dims));
            minx = min(minx, Iv.inf( dims(1) )); maxx = max(maxx, Iv.sup( dims(1) ));
            miny = min(miny, Iv.inf( dims(2) )); maxy = max(maxy, Iv.sup( dims(2) ));
        end
    end
    if ~isempty(Ymeas)
        minx = min(minx, min(Ymeas(:,dims(1)))); maxx = max(maxx, max(Ymeas(:,dims(1))));
        miny = min(miny, min(Ymeas(:,dims(2)))); maxy = max(maxy, max(Ymeas(:,dims(2))));
    end
    if ~isfinite(minx) || ~isfinite(miny), minx=-1; maxx=1; miny=-1; maxy=1; end
    pad = 0.08;
    sx = maxx-minx; sy = maxy-miny;
    lims = [minx-pad*(sx+eps), maxx+pad*(sx+eps), miny-pad*(sy+eps), maxy+pad*(sy+eps)];
end

function v = iff(cond,a,b), if cond, v=a(); else, v=b(); end, end
