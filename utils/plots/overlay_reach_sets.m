function overlay_reach_sets(artifact_mat, varargin)
% Overlays True, DDRA, Gray output sets at selected time steps.
% Also scatters validation outputs for those steps.

% Args (name/value):
%   'Dims'     : 2-D output indices to plot; if [], auto-pick from output energy
%   'Klist'    : vector of time indices (1-based), default [1, round(nk/2), nk]
%   'Reduce'   : zonotope order cap for reductions (Girard), default 60
%   'SaveBase' : base path to save figs (optional)
%   'Show'     : bool

p = inputParser;
addParameter(p,'Dims',[]);
addParameter(p,'Klist',[]);
addParameter(p,'Reduce',60);
addParameter(p,'SaveBase',"");
addParameter(p,'Show',true);
parse(p,varargin{:});
prm = p.Results;

S = load(artifact_mat);
sysT = S.sys_ddra; 
sysG = S.sys_gray; 
VAL  = S.VAL;         % contains R0, x0{b}, u{b}, y{b}, U (input set!)
W_eff= S.W_eff;       
M_AB = S.M_AB;

% Disturbance sets mapped like in run_sweeps
W_true = coerceWToSys(sysT, W_eff);
W_pred = [];
if isprop(sysG,'nrOfDisturbances') && sysG.nrOfDisturbances>0
    W_pred = normalizeWForGray(sysG, W_eff);
end

Kred = prm.Reduce;
nkv  = size(VAL.y{1},1);
if isempty(prm.Klist), prm.Klist = unique([1, round(nkv/2), nkv]); end

% --- pick a representative block 
b = 1;
uk_all = VAL.u{b};         % (n_k x m)

% ---------- Reach: TRUE ----------
% IMPORTANT: include U (input set), like in run_sweeps
params_true = struct('R0', VAL.R0, 'U', VAL.U, 'u', uk_all', ...
                     'tFinal', sysT.dt*(nkv-1));
if ~isempty(W_true), params_true.W = W_true; end
optReach = struct('zonotopeOrder',Kred,'reductionTechnique','girard');
Rt = reach(sysT, params_true, optReach).timePoint.set;   % cell{k} in state space

% ---------- Reach: GRAY ----------
% match run_sweeps: zero-centered R0 + x0 offset for block b
R0_gray = zonotope(zeros(size(center(VAL.R0))), generators(VAL.R0)) + VAL.x0{b};
params_gray = struct('R0', R0_gray, 'u', uk_all', 'tFinal', sysG.dt*(nkv-1));
if ~isempty(W_pred), params_gray.W = W_pred; end
Rg = reach(sysG, params_gray, optReach).timePoint.set;   % cell{k} in state space

% ---------- DDRA PRE-update walk ----------
Xk = reduce(VAL.R0 + VAL.x0{b}, 'girard', Kred);
Yd = cell(nkv,1);

is_meas = isfield(S,'AB_meas_center') && isfield(S,'AV_oneterm');
if is_meas
    ABc   = S.AB_meas_center;
    AVone = S.AV_oneterm;
    Vmeas = getfield(S,'V_meas', zonotope(zeros(size(center(VAL.R0),1),1)));
end

for k = 1:nkv
    Xk = reduce(Xk, 'girard', Kred);
    uk = VAL.u{b}(k,:).';

    if ~is_meas
        % standard
        Yd{k} = sysT.C*Xk + sysT.D*uk;
        Xk    = S.M_AB * cartProd(Xk, zonotope(uk)) + S.W_eff;
    else
        % meas-noise aware (Alg. 4)
        Yd{k} = sysT.C*Xk + sysT.D*uk;
        X_for_AB = Xk;
        if k>1 && ~isempty(generators(Vmeas))
            X_for_AB = minkowskiSum(Xk, Vmeas);
        end
        Xk = (ABc * cartProd(X_for_AB, zonotope(uk))) + AVone + S.W_eff;
    end
end


% ---------- Build output sets for TRUE/GRAY (coerced to zonotopes) ----------
Yt = cell(nkv,1);
Yg = cell(nkv,1);
for k = 1:nkv
    uk = uk_all(k,:).';
    Yt{k} = asZono( linearMap(Rt{k}, sysT.C) + sysT.D*uk );
    Yg{k} = asZono( linearMap(Rg{k}, sysG.C) + sysG.D*uk );
end

% ---------- Choose output dims ----------
ny = size(sysT.C,1);
dims = prm.Dims(:).';
if isempty(dims)
    dims = pick_output_dims(Yt);  % auto from True output richness
end
if numel(dims)~=2
    dims = [1, min(2,ny)];
end

% ---------- plotting ----------
ncols = numel(prm.Klist);
f = figure('Color','w','Name','Overlay of reachable sets');
tiledlayout(1,ncols,'TileSpacing','compact','Padding','compact');

for ii = 1:ncols
    k = prm.Klist(ii); k = max(1, min(nkv, k));
    nexttile; hold on; grid on; axis equal;
    title(sprintf('k = %d', k));

    % True
    plot_set_2d(Yt{k}, dims, [0.2 0.6 0.2], 0.25, 'True');

    % DDRA
    plot_set_2d(Yd{k}, dims, [0.2 0.2 0.7], 0.25, 'DDRA');

    % Gray
    plot_set_2d(Yg{k}, dims, [0.8 0.3 0.1], 0.25, 'Gray');

    % Samples
    scatter(VAL.y{b}(k,dims(1)), VAL.y{b}(k,dims(2)), 18, 'k', 'filled', 'MarkerFaceAlpha',0.7);

    xlabel(sprintf('y_{%d}',dims(1))); ylabel(sprintf('y_{%d}',dims(2)));
    legend('Location','bestoutside'); box on;
end

if strlength(prm.SaveBase)>0
    saveas(f, prm.SaveBase + "_overlays.png");
    try, exportgraphics(f, prm.SaveBase + "_overlays.pdf",'ContentType','vector'); catch, end
end
if ~prm.Show, close(f); end
end

% ====================== helpers ======================

function Z = asZono(S)
% Coerce any set S to a simple zonotope (outer-approx) that plot() understands.
    try
        if isa(S,'zonotope'); Z = reduce(S,'girard',max(1,get_zonotope_order(S))); return; end
    catch, end
    try, Z = toZono(S);  % works for polyZonotope → zonotope
    catch
        % last resort: interval hull → zonotope
        Iv = interval(S);
        try, lo = infimum(Iv); hi = supremum(Iv);
        catch, lo = Iv.inf;    hi = Iv.sup;
        end
        c  = (double(lo)+double(hi))/2;
        g  = (double(hi)-double(lo))/2;
        Z  = zonotope(c, diag(g));
    end
    % be nice on order
    Z = reduce(Z, 'girard', 60);
end

function ord = get_zonotope_order(Z)
    G = generators(Z); if isempty(G), ord = 1; else, ord = size(G,2); end
end

function dims = pick_output_dims(Ycell)
% pick two output indices with the largest generator energy at a representative step
    k = min(numel(Ycell), max(1, round(numel(Ycell)/2)));
    Z = asZono(Ycell{k});
    G = generators(Z);
    if isempty(G)
        ny = size(center(Z),1); dims = [1, min(2,ny)]; return;
    end
    e = sum(G.^2,2);
    [~,I] = sort(e,'descend');
    if numel(I)>=2, dims = sort(I(1:2)); else, dims = [I(1), max(1,I(1)+1)]; end
end

function plot_set_2d(S, dims, col, alpha, name)
    % S is already a zonotope in output space; project & plot robustly
    ny = size(center(S),1);
    dims = max(1,min(ny,dims));
    P = zeros(2, ny); P(1,dims(1))=1; P(2,dims(2))=1;
    Y = asZono(linearMap(S, P));   % ensure 2D zono
    drew = false;
    try
        h = plot(Y, [1 2]); 
        if ~isempty(h) && isgraphics(h(1))
            set(h(1),'FaceColor',col,'FaceAlpha',alpha,'EdgeColor','none');
            drew = true;
        end
    catch, drew = false;
    end
    if ~drew
        outline_zono(Y, col, alpha);
    end
    if nargin>=5
        plot(NaN,NaN,'s','MarkerFaceColor',col,'MarkerEdgeColor',col,'DisplayName',name);
    end
end

function outline_zono(Z, col, alpha)
    % convex hull via sampled support (always works for 2D)
    th = linspace(0,2*pi,128);
    D  = [cos(th); sin(th)];
    c  = center(Z); G = generators(Z);
    s  = D' * c + sum(abs(D' * G),2);
    P  = (D .* s')';
    K  = convhull(P(:,1), P(:,2));
    patch(P(K,1), P(K,2), col, 'FaceAlpha',alpha, 'EdgeColor','none');
end
