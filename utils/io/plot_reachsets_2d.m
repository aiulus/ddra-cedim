function plot_reachsets_2d(sys_gray, sys_ddra, VAL, W_eff, b, dims, varargin)
% Plot 2D OUTPUT reachable sets (Gray/RCSI vs true/DDRA) for block b.
% Robust to VAL.u orientation. Adds options for legend clarity, nominal traces,
% and visual inflation to show filled sets when uncertainty is zero.

opt = struct('Colors', [], 'SaveDir', '', 'TikZ', false, 'Name', 'reach_2d', ...
             'ShowNominal', true, 'VisualInflate', 0);   % <-- new knobs
if ~isempty(varargin), opt = parseOpts(opt, varargin{:}); end
if isempty(opt.Colors), C = colorscheme('tum'); else, C = opt.Colors; end

% small helpers for darker lines
darken = @(rgb, f) max(min(rgb * f,1),0);
if ~isfield(C,'ddra_dark'), C.ddra_dark = darken(C.ddra,0.7); end
if ~isfield(C,'rcsi_dark'), C.rcsi_dark = darken(C.rcsi,0.7); end
if ~isfield(C,'nominal'),   C.nominal   = [0 0 0];          end

% sanity on dims
ny = sys_gray.nrOfOutputs;
assert(numel(dims)==2 && all(dims>=1) && all(dims<=ny), ...
    'dims must be a 2-vector within 1..nrOfOutputs');

% build CORA system for the "true"/DDRA side (use sys_gray.dt for consistency)
sys_true = linearSysDT(linearSys(sys_ddra.A, sys_ddra.B, [], sys_ddra.C, sys_ddra.D), sys_gray.dt);

% pull data for block b and coerce shapes
x0 = VAL.x0{b};             % (n_x × 1)
u  = VAL.u{b};              % (m × n_k) or (n_k × m)
m  = sys_gray.nrOfInputs;
if size(u,1) ~= m && size(u,2) == m, u = u.'; end
assert(size(u,1) == m, 'VAL.u{%d} has wrong shape; expected m=%d.', b, m);
n_k = size(u,2);
dt  = sys_true.dt;

% reachability params (no input uncertainty)
params = struct();
params.tStart = 0;
params.tFinal = dt*(n_k-1);
params.R0     = zonotope(x0);
params.U      = zonotope(zeros(m,1));
params.u      = u;                    

% only pass W if it matches; otherwise skip (keeps reach valid)
if ~isempty(W_eff)
    try nW = size(center(W_eff),1); catch, nW = NaN; end
    if isequal(sys_gray.nrOfDisturbances, nW)
        params.W = W_eff;
    end
end

options = struct('zonotopeOrder', inf, 'verbose', false);

% state reach
Rg = reach(sys_gray, params, options);
Rt = reach(sys_true, params, options);
Xg = Rg.timePoint.set;   % Gray/RCSI
Xt = Rt.timePoint.set;   % true/DDRA

% map to output: y = C x + D u_k
Cmat = sys_ddra.C;  Dmat = sys_ddra.D;
if isempty(Dmat), Dmat = zeros(size(Cmat,1), m); end

Yg = cell(1,n_k); Yt = cell(1,n_k);
for k = 1:n_k
    uk  = u(:,k);
    off = Dmat*uk;
    Yg{k} = map_to_output_zono(Xg{k}, Cmat, off);
    Yt{k} = map_to_output_zono(Xt{k}, Cmat, off);

    % optional: visually inflate in OUTPUT space so sets show as polygons
    if opt.VisualInflate > 0
        Yg{k} = inflate_output_dims(Yg{k}, dims, opt.VisualInflate, ny);
        Yt{k} = inflate_output_dims(Yt{k}, dims, opt.VisualInflate, ny);
    end
end

% optional nominal trajectories (deterministic)
if opt.ShowNominal
    xg = zeros(size(x0,1), n_k); xg(:,1) = x0;
    xt = xg;
    for k = 1:n_k-1
        xg(:,k+1) = sys_gray.A*xg(:,k) + sys_gray.B*u(:,k);
        xt(:,k+1) = sys_ddra.A*xt(:,k) + sys_ddra.B*u(:,k);
    end
    yg = sys_ddra.C*xg + sys_ddra.D*u;
    yt = sys_ddra.C*xt + sys_ddra.D*u;
end

% plot
figure('Color','w'); hold on; box on;
title(sprintf('Output sets (block %d), dims [%d %d]', b, dims(1), dims(2)));
xlabel(sprintf('y_{%d}', dims(1))); ylabel(sprintf('y_{%d}', dims(2)));

% pre-create legend prototypes (stable icons)
hTrue = []; hGray = []; hSamp = []; hNomT = []; hNomG = [];

% DDRA/true sets (blue)
for k = 1:n_k
    h = plot_set2d(Yt{k}, dims, C.ddra, 0.12);
    if isempty(hTrue), set(h,'DisplayName','DDRA/true reach'); hTrue = h;
    else, set(h,'HandleVisibility','off'); end
end

% Gray/RCSI sets (orange)
for k = 1:n_k
    h = plot_set2d(Yg{k}, dims, C.rcsi, 0.12);
    if isempty(hGray), set(h,'DisplayName','Gray/RCSI reach'); hGray = h;
    else, set(h,'HandleVisibility','off'); end
end

% measured samples (if any)
if isfield(VAL,'y') && numel(VAL.y) >= b && ~isempty(VAL.y{b})
    Ymeas = squeeze(VAL.y{b});
    if size(Ymeas,1)~=n_k && size(Ymeas,2)==n_k, Ymeas = Ymeas.'; end
    h = plot(Ymeas(:,dims(1)), Ymeas(:,dims(2)), 'o-', 'LineWidth',1.1, ...
             'Color', C.nominal, 'MarkerSize',4);
    set(h,'DisplayName','samples'); hSamp = h;
end

% nominal tracks (optional)
if opt.ShowNominal
    h = plot(yt(dims(1),:), yt(dims(2),:), '-',  'Color', C.ddra_dark, 'LineWidth',1.2);
    if isempty(hNomT), set(h,'DisplayName','DDRA nominal'); hNomT = h;
    else, set(h,'HandleVisibility','off'); end

    h = plot(yg(dims(1),:), yg(dims(2),:), '--', 'Color', C.rcsi_dark, 'LineWidth',1.2);
    if isempty(hNomG), set(h,'DisplayName','RCSI nominal'); hNomG = h;
    else, set(h,'HandleVisibility','off'); end
end

legend('Location','best');

% save?
if ~isempty(opt.SaveDir)
    if ~exist(opt.SaveDir,'dir'), mkdir(opt.SaveDir); end
    base = fullfile(opt.SaveDir, opt.Name);
    exportgraphics(gcf, [base '.png'], 'Resolution', 200);
    exportgraphics(gcf, [base '.pdf'], 'ContentType','vector');
    if opt.TikZ && exist('matlab2tikz','file')==2
        matlab2tikz([base '.tex'], 'height','\figureheight','width','\figurewidth');
    end
end
end

% --- helpers -------------------------------------------------------------
function ZY = map_to_output_zono(ZX, Cmat, d)
ZX = zonotope(ZX);
c  = center(ZX);  G = generators(ZX);
ZY = zonotope(Cmat*c + d, Cmat*G);
end

function Zout = inflate_output_dims(Zin, dims, eps, ny)
% add tiny generators in the requested output dims (visual-only thickness)
Gx = zeros(ny, numel(dims));
for i = 1:numel(dims), Gx(dims(i), i) = eps; end
Zout = Zin + zonotope(zeros(ny,1), Gx);
end

function o = parseOpts(o, varargin)
for i=1:2:numel(varargin), k = varargin{i}; v = varargin{i+1};
    if isfield(o,k), o.(k) = v; end
end
end

function h = plot_set2d(S, dims, col, alpha)
% Robust 2D plot for CORA contSets; returns a handle for legend control.
if representsa(S,'emptySet'); h = plot(nan,nan); return; end
Iv = interval(S); w = Iv.sup - Iv.inf; w2 = w(dims); tol = 1e-12;
if all(w2 <= tol)
    h = plot(S, dims, 'Color', col, 'Marker','.', 'LineStyle','none', 'LineWidth', 1.0);
elseif any(w2 <= tol)
    h = plot(S, dims, 'Color', col, 'LineWidth', 1.0);
else
    h = plot(S, dims, 'FaceColor', col, 'EdgeColor', col, 'FaceAlpha', alpha, 'LineWidth', 0.6);
end
end
