function overlay_reach_sets(artifact_mat, varargin)
% Overlays True, DDRA, Gray output sets at selected time steps.
% Also scatters validation outputs for those steps.

% Args (name/value):
%   'Dims'     : 2-D output indices to plot, default [1 2]
%   'Klist'    : vector of time indices (1-based), default [1, round(nk/2), nk]
%   'Reduce'   : zonotope order for reductions, default 60
%   'SaveBase' : base path to save figs (optional)
%   'Show'     : bool

p = inputParser;
addParameter(p,'Dims',[1 2]);
addParameter(p,'Klist',[]);
addParameter(p,'Reduce',60);
addParameter(p,'SaveBase',"");
addParameter(p,'Show',true);
parse(p,varargin{:});
prm = p.Results;

S = load(artifact_mat);
sysT = S.sys_ddra; sysG = S.sys_gray; VAL = S.VAL;
W_eff = S.W_eff;   M_AB = S.M_AB;

W_true = coerceWToSys(sysT, W_eff);
W_pred = []; if isprop(sysG,'nrOfDisturbances') && sysG.nrOfDisturbances>0
    W_pred = normalizeWForGray(sysG, W_eff);
end

Kred = prm.Reduce;
nkv  = size(VAL.y{1},1);
if isempty(prm.Klist), prm.Klist = unique([1, round(nkv/2), nkv]); end
dims = prm.Dims(:).';

% pick a representative block 
b = 1;
params_true = struct('R0', VAL.R0, 'u', VAL.u{b}', 'tFinal', sysT.dt*(nkv-1));
if ~isempty(W_true), params_true.W = W_true; end
Rt = reach(sysT, params_true, struct('zonotopeOrder',Kred,'reductionTechnique','girard')).timePoint.set;

R0_gray = zonotope(zeros(size(center(VAL.R0))), generators(VAL.R0)) + VAL.x0{b};
params_gray = struct('R0', R0_gray, 'u', VAL.u{b}', 'tFinal', sysG.dt*(nkv-1));
if ~isempty(W_pred), params_gray.W = W_pred; end
Rg = reach(sysG, params_gray, struct('zonotopeOrder',Kred,'reductionTechnique','girard')).timePoint.set;

% DDRA walk
Xk = reduce(VAL.R0 + VAL.x0{b}, 'girard', Kred);
Yd = cell(nkv,1);
for k = 1:nkv
    Xk = reduce(Xk, 'girard', Kred);
    uk = VAL.u{b}(k,:).';
    Yd{k} = sysT.C*Xk + sysT.D*uk;          % pre-update output set
    Xk = M_AB * cartProd(Xk, zonotope(uk)) + W_eff;
end

% plotting
ncols = numel(prm.Klist);
f = figure('Color','w','Name','Overlay of reachable sets');
tiledlayout(1,ncols,'TileSpacing','compact','Padding','compact');

for ii = 1:ncols
    k = prm.Klist(ii); k = max(1, min(nkv, k));
    nexttile; hold on; grid on; axis equal;
    title(sprintf('k = %d', k));

    % True
    plot_set_2d(sysT.C*Rt{k} + sysT.D*VAL.u{b}(k,:).', dims, [0.2 0.6 0.2], 0.25, 'True');

    % DDRA
    plot_set_2d(Yd{k}, dims, [0.2 0.2 0.7], 0.25, 'DDRA');

    % Gray
    plot_set_2d(sysG.C*Rg{k} + sysG.D*VAL.u{b}(k,:).', dims, [0.8 0.3 0.1], 0.25, 'Gray');

    % Samples
    scatter(VAL.y{b}(k,dims(1)), VAL.y{b}(k,dims(2)), 18, 'k', 'filled', 'MarkerFaceAlpha',0.7);

    xlabel(sprintf('y_{%d}',dims(1))); ylabel(sprintf('y_{%d}',dims(2)));
    legend('Location','bestoutside'); box on;
end

if strlength(prm.SaveBase)>0
    saveas(f, prm.SaveBase + "_overlays.png");
    try, exportgraphics(f, prm.SaveBase + "_overlays.pdf",'ContentType','vector'); catch, end
end
end

% ------ helpers ------
function plot_set_2d(S, dims, col, alpha, name)
    % Project via linearMap onto dims
    P = zeros(2, size(center(S),1)); P(1,dims(1))=1; P(2,dims(2))=1;
    Y = linearMap(S, P);
    try
        h = plot(Y, [1 2]); 
        if ~isempty(h) && isgraphics(h(1)), set(h(1),'FaceColor',col,'FaceAlpha',alpha,'EdgeColor','none'); end
    catch
        % fallback: outline via sampled support (rarely needed if CORA handles plot)
        outline_zono(Y, col, alpha);
    end
    if nargin>=5, plot(NaN,NaN,'s','MarkerFaceColor',col,'MarkerEdgeColor',col,'DisplayName',name); end
end

function outline_zono(Z, col, alpha)
    % crude convex hull via support on a circle (used only if CORA plot fails)
    th = linspace(0,2*pi,128);
    D  = [cos(th); sin(th)];
    Zz = zonotope(Z);
    s  = D' * center(Zz) + sum(abs(D' * generators(Zz)),2);
    P  = (D .* s')';
    K  = convhull(P(:,1), P(:,2));
    patch(P(K,1), P(K,2), col, 'FaceAlpha',alpha, 'EdgeColor','none');
end
