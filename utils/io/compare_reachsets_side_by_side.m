function compare_reachsets_side_by_side(sys_gray, sys_ddra, VAL, W_eff, b, dims, varargin)
%COMPARE_REACHSETS_SIDE_BY_SIDE
%  Two-panel figure: (Left) DDRA/true reach; (Right) Gray/RCSI reach.
%  Uses same nominal input, identical axes, optional nominal curves & samples.
%
%  Required:
%    sys_gray  : CORA system returned by Gray/RCSI
%    sys_ddra  : "true" system used to generate data (DDRA nominal)
%    VAL       : struct with fields x0{b}, u{b} (and optionally y{b})
%    W_eff     : [] or a zonotope (only used if disturbance dim matches sys)
%    b         : block index into VAL
%    dims      : [i j] output dims to plot
%
%  Options (Name,Value):
%    'Colors'       : struct from colorscheme('tum') (default auto)
%    'SaveDir'      : folder to save figures (default '')
%    'TikZ'         : true/false (default false)
%    'Name'         : base filename (default 'reach_side_by_side')
%    'ShowNominal'  : true/false (default true)
%    'VisualInflate': scalar eps to add in output dims for visible area (default 1e-2)
%    'Every'        : plot every k-th step to declutter (default 1 = all)
%    'M_AB'         : matZonotope for DDRA recursion (default [])
%    'OverlayTrue'  : overlay true sets on left when DDRA is shown (default false)

opt = struct('Colors', [], 'SaveDir', '', 'TikZ', false, 'Name', 'reach_side_by_side', ...
             'ShowNominal', true, 'VisualInflate', 1e-2, 'Every', 1, ...
             'M_AB', [], 'OverlayTrue', false);   % %% NEW

if ~isempty(varargin), opt = parseOpts(opt, varargin{:}); end
if isempty(opt.Colors), C = colorscheme('tum'); else, C = opt.Colors; end
darken = @(rgb,f) max(min(rgb*f,1),0);
if ~isfield(C,'ddra_dark'), C.ddra_dark = darken(C.ddra,0.7); end
if ~isfield(C,'rcsi_dark'), C.rcsi_dark = darken(C.rcsi,0.7); end
if ~isfield(C,'nominal'),   C.nominal   = [0 0 0];         end

% sanity
ny = sys_gray.nrOfOutputs;
assert(numel(dims)==2 && all(dims>=1) && all(dims<=ny), 'dims must be 2-vector within 1..ny');

% build a "true" CORA object with gray’s dt to keep horizons aligned
sys_true = linearSysDT(linearSys(sys_ddra.A, sys_ddra.B, [], sys_ddra.C, sys_ddra.D), sys_gray.dt);

% pull block data, ensure u shape
x0 = VAL.x0{b};
u  = VAL.u{b};              % (m×n_k) or (n_k×m)
m  = sys_gray.nrOfInputs;
if size(u,1)~=m && size(u,2)==m, u = u.'; end
assert(size(u,1)==m, 'VAL.u{%d} must be m×n_k; got %dx%d', b, size(u,1), size(u,2));
n_k = size(u,2); dt = sys_true.dt;

% cora params (deterministic R0/U; we inflate later for visibility)
params = struct('tStart',0, 'tFinal',dt*(n_k-1), ...
                'R0', zonotope(x0), ...
                'U',  zonotope(zeros(m,1)), ...
                'u',  u.');                     % CORA expects n_k×m
% include W only if dims match disturbances
if ~isempty(W_eff)
    try nW = size(center(W_eff),1); catch, nW = NaN; end
    if isequal(sys_gray.nrOfDisturbances, nW), params.W = W_eff; end
end
options = struct('zonotopeOrder', inf, 'verbose', false);
params.u = params.u.';

% reach states
Rg = reach(sys_gray, params, options);
Rt = reach(sys_true, params, options);
Xg = Rg.timePoint.set;
Xt = Rt.timePoint.set;

% map states -> outputs: y = Cx + D u_k
Cmat = sys_ddra.C; Dmat = sys_ddra.D; if isempty(Dmat), Dmat = zeros(size(Cmat,1), m); end
[Yg, Yt] = map_to_output_traces(Xg, Xt, Cmat, Dmat, u);
% --- DDRA (data-driven) output sets if M_AB is available ----------------  
Ydd = [];  % will remain empty if M_AB not provided or fails
if ~isempty(opt.M_AB)
    try
        Xdd = cell(1,n_k); Xdd{1} = zonotope(x0);      % start from point initial state
        for k = 1:n_k-1
            Uk  = zonotope(u(:,k));                   % singleton input {u_k}
            Zk  = cartProd(Xdd{k}, Uk);               % [x;u] product
            Xk1 = opt.M_AB * Zk;                      % apply learned matrix zonotope
            if ~isempty(W_eff), Xk1 = Xk1 + W_eff; end
            % optionally: Xk1 = reduce(Xk1,'girard',200);
            Xdd{k+1} = Xk1;
        end
        Ydd = cell(1,n_k);
        for k = 1:n_k
            off   = Dmat*u(:,k);
            Ydd{k}= map_to_output_zono(Xdd{k}, Cmat, off);
        end
    catch ME
        warning('compare_reachsets_side_by_side:DDRAPlot', ...
                'DDRA plot skipped (M_AB path): %s', ME.message);
        Ydd = [];
    end
end

% visual inflation for filled polygons even if uncertainty is zero
if opt.VisualInflate > 0
    for k=1:n_k
        Yg{k} = inflate_output_dims(Yg{k}, dims, opt.VisualInflate, ny);
        Yt{k} = inflate_output_dims(Yt{k}, dims, opt.VisualInflate, ny);
        if ~isempty(Ydd), Ydd{k} = inflate_output_dims(Ydd{k}, dims, opt.VisualInflate, ny); end  % %% NEW
    end
end


% nominal trajectories (deterministic)
if opt.ShowNominal
    xg = zeros(size(x0,1), n_k); xg(:,1) = x0;
    xt = xg;
    for k=1:n_k-1
        xg(:,k+1) = sys_gray.A*xg(:,k) + sys_gray.B*u(:,k);
        xt(:,k+1) = sys_ddra.A*xt(:,k) + sys_ddra.B*u(:,k);
    end
    yg = sys_ddra.C*xg + sys_ddra.D*u;
    yt = sys_ddra.C*xt + sys_ddra.D*u;
end

% measured outputs if present
hasY = isfield(VAL,'y') && numel(VAL.y)>=b && ~isempty(VAL.y{b});
if hasY
    Ymeas = squeeze(VAL.y{b}); if size(Ymeas,1)~=n_k && size(Ymeas,2)==n_k, Ymeas = Ymeas.'; end
end

% compute shared axes from both sides (+samples)
lims = auto_axes_from_sets([Yt Yg Ydd], dims, hasY, tern(hasY, Ymeas, []));   % %% CHANGED

% ==== figure layout ====
figure('Color','w'); tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% LEFT: DDRA/true
% LEFT: DDRA (data-driven) if available;
nexttile; hold on; box on;
leftTitle = 'True model';
if ~isempty(Ydd), leftTitle = 'DDRA (data-driven)'; end
title(leftTitle);
xlabel(sprintf('y_{%d}',dims(1))); ylabel(sprintf('y_{%d}',dims(2)));

hDD = []; hTrue = []; hNomT = []; hSamp = [];
if ~isempty(Ydd)
    % DDRA reach
    for k = 1:opt.Every:n_k
        h = plot_set2d(Ydd{k}, dims, C.ddra, 0.14);
        if isempty(hDD), set(h,'DisplayName','DDRA reach'); hDD = h; else, set(h,'HandleVisibility','off'); end
    end
    % optional overlay: true model outline
    if opt.OverlayTrue
        for k = 1:opt.Every:n_k
            plot_set2d(Yt{k}, dims, C.rcsi_dark, 0.06);  % faint
        end
    end
else
    % fallback: True model reach (what you had before)
    for k = 1:opt.Every:n_k
        h = plot_set2d(Yt{k}, dims, C.ddra, 0.14);
        if isempty(hTrue), set(h,'DisplayName','True reach'); hTrue = h; else, set(h,'HandleVisibility','off'); end
    end
end

if opt.ShowNominal
    % nominal for the "true" side (left panel baseline)
    h = plot(yt(dims(1),:), yt(dims(2),:), '-', 'Color', C.ddra_dark, 'LineWidth',1.3);
    if isempty(hNomT), set(h,'DisplayName','DDRA/True nominal'); hNomT = h; else, set(h,'HandleVisibility','off'); end
end
if hasY
    h = plot(Ymeas(:,dims(1)), Ymeas(:,dims(2)), 'o-', 'Color', C.nominal, 'LineWidth',1.0, 'MarkerSize',4);
    if isempty(hSamp), set(h,'DisplayName','samples'); hSamp = h; else, set(h,'HandleVisibility','off'); end
end
axis(lims); legend('Location','best');


% RIGHT: Gray/RCSI
nexttile; hold on; box on; title('Gray / RCSI'); xlabel(sprintf('y_{%d}',dims(1))); ylabel(sprintf('y_{%d}',dims(2)));
hGray = []; hNomG = [];
for k = 1:opt.Every:n_k
    h = plot_set2d(Yg{k}, dims, C.rcsi, 0.14);
    if isempty(hGray), set(h,'DisplayName','Gray reach'); hGray = h; else, set(h,'HandleVisibility','off'); end
end
if opt.ShowNominal
    h = plot(yg(dims(1),:), yg(dims(2),:), '--', 'Color', C.rcsi_dark, 'LineWidth',1.3);
    if isempty(hNomG), set(h,'DisplayName','RCSI nominal'); hNomG = h; else, set(h,'HandleVisibility','off'); end
end
if hasY
    h = plot(Ymeas(:,dims(1)), Ymeas(:,dims(2)), 'o-', 'Color', C.nominal, 'LineWidth',1.0, 'MarkerSize',4);
    set(h,'HandleVisibility','off'); % avoid duplicate legend entry
end
axis(lims); legend('Location','best');

% save outputs
if ~isempty(opt.SaveDir)
    if ~exist(opt.SaveDir,'dir'), mkdir(opt.SaveDir); end
    base = fullfile(opt.SaveDir, opt.Name);
    exportgraphics(gcf, [base '.png'], 'Resolution', 200);
    exportgraphics(gcf, [base '.pdf'], 'ContentType','vector');
    if opt.TikZ && exist('matlab2tikz','file')==2
        matlab2tikz([base '.tex'], 'height','\figureheight','width','\figurewidth');
    end
end

end % === main ===

% ------------ helpers ------------
function [Yg, Yt] = map_to_output_traces(Xg, Xt, C, D, u)
n_k = numel(Xg); m = size(u,1);
if isempty(D), D = zeros(size(C,1), m); end
Yg = cell(1,n_k); Yt = cell(1,n_k);
for k=1:n_k
    off = D*u(:,k);
    Yg{k} = map_to_output_zono(Xg{k}, C, off);
    Yt{k} = map_to_output_zono(Xt{k}, C, off);
end
end

function ZY = map_to_output_zono(ZX, C, d)
ZX = zonotope(ZX);
ZY = zonotope(C*center(ZX) + d, C*generators(ZX));
end

function Zout = inflate_output_dims(Zin, dims, eps, ny)
G = zeros(ny, numel(dims));
for i=1:numel(dims), G(dims(i),i) = eps; end
Zout = Zin + zonotope(zeros(ny,1), G);
end

function lims = auto_axes_from_sets(Ycells, dims, hasY, Ymeas)
minx = inf; maxx = -inf; miny = inf; maxy = -inf;
for i=1:numel(Ycells)
    S = Ycells{i}; if isempty(S), continue; end
    Iv = interval(project(S,dims));
    minx = min(minx, infimum(Iv(1))); maxx = max(maxx, supremum(Iv(1)));
    miny = min(miny, infimum(Iv(2))); maxy = max(maxy, supremum(Iv(2)));
end
if hasY
    minx = min(minx, min(Ymeas(:,dims(1)))); maxx = max(maxx, max(Ymeas(:,dims(1))));
    miny = min(miny, min(Ymeas(:,dims(2)))); maxy = max(maxy, max(Ymeas(:,dims(2))));
end
padx = 0.1*(maxx-minx + eps); pady = 0.1*(maxy-miny + eps);
lims = [minx-padx, maxx+padx, miny-pady, maxy+pady];
end

function h = plot_set2d(S, dims, col, alpha)
if representsa(S,'emptySet'), h = plot(nan,nan); return; end
Iv = interval(S); w = Iv.sup - Iv.inf; w2 = w(dims); tol = 1e-12;
if all(w2 <= tol)
    h = plot(S, dims, 'Color', col, 'Marker','.', 'LineStyle','none', 'LineWidth', 1.0);
elseif any(w2 <= tol)
    h = plot(S, dims, 'Color', col, 'LineWidth', 1.0);
else
    h = plot(S, dims, 'FaceColor', col, 'EdgeColor', col, 'FaceAlpha', alpha, 'LineWidth', 0.6);
end
end

function v = tern(c,a,b), if c, v=a; else, v=b; end, end

function o = parseOpts(o, varargin)
for i=1:2:numel(varargin), k=varargin{i}; v=varargin{i+1};
    if isfield(o,k), o.(k) = v; end
end
end
