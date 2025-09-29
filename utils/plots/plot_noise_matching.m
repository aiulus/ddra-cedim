function plot_noise_matching(artifact_mat, varargin)
% Visual confirmation of noise matching:
%   (i)  overlay W_eff (true state disturbance) vs E_g * W_pred (Gray pushforward)
%   (ii) per-state interval half-width bars
%   (iii) directional support comparison and max deficiency
%
% Inputs
%   artifact_mat : path to a .mat saved by run_sweeps (contains sys_gray, W_eff; W_pred is
%                  recomputed if not present)
% Name/Value
%   'StateDims'  : 2 indices for state projection overlay (auto if []), default []
%   'Nd'         : # of directions for directional support, default 64
%   'Title'      : figure title suffix, default ""
%   'SaveBase'   : base path (no ext) to save png/pdf, default ""
%
% Requires CORA's zonotope ops: center, generators, linearMap, interval
p = inputParser;
addParameter(p,'StateDims',[]);
addParameter(p,'Nd',64);
addParameter(p,'Title',"");
addParameter(p,'SaveBase',"");
parse(p,varargin{:});
prm = p.Results;

S = load(artifact_mat);

% --- Pull what we need from the artifact
sysG  = S.sys_gray;
W_eff = S.W_eff;                      % true/DDRA state disturbance
W_pred = [];                          % Gray disturbance-space zonotope (n_w)
if isfield(S,'W_pred'),  W_pred = S.W_pred; end
if isempty(W_pred)
    % Rebuild Gray disturbance set from state-space W (repo helper)
    try
        W_pred = normalizeWForGray(sysG, W_eff);
    catch
        % If normalizeWForGray not on path, try build_W_pred (also in repo)
        try
            U = S.VAL.U;                      %#ok<NASGU> % some versions need U shape
            W_pred = build_W_pred(sysG, [], W_eff);
        catch
            warning('Could not reconstruct Gray disturbance set (W_pred).'); 
            W_pred = [];
        end
    end
end

% If Gray has no disturbance channels, just show W_eff bars and return
hasDist = isprop(sysG,'nrOfDisturbances') && sysG.nrOfDisturbances>0;
if ~hasDist || isempty(W_pred)
    figure('Color','w','Name','Noise matching');
    t = tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
    title(t, compose_title("Noise matching", prm.Title, " (Gray has no disturbance channels)"));
    bar_state_halfwidths(W_eff, [], 'True W_{eff}', 'Gray E_gZ_d (n/a)');
    if strlength(prm.SaveBase)>0
        saveas(gcf, prm.SaveBase + "_noise_match.png");
        try, exportgraphics(gcf, prm.SaveBase + "_noise_match.pdf",'ContentType','vector'); catch, end
    end
    return;
end

% --- Map Gray disturbance to state space: E_g * W_pred
Eg = []; 
try
    Eg = sysG.E;               % CORA linearSysDT uses E for disturbances
catch
    % Some versions store disturbance matrix as Bd or Ed; attempt fallback
    if isprop(sysG,'Bd'), Eg = sysG.Bd; end
end
if isempty(Eg)
    error('Could not access Gray disturbance injection matrix E_g.');
end
Wg_state = linearMap(W_pred, Eg);   % state-space disturbance induced by Gray

% --- Decide state dims for overlay
nx = size(center(W_eff),1);
dims = prm.StateDims;
if isempty(dims)
    dims = pick_dims_for_plot(W_eff);     % 2 strongest rows in generator energy
end

% --- Compute coordinate half-widths (interval radii)
[hw_true, hw_gray] = coord_halfwidths(W_eff, Wg_state);

% --- Directional support comparison (in state space)
Nd = max(8, round(prm.Nd));
Ddirs = sample_unit_dirs(nx, Nd, 12345);
[su_true, su_gray] = support_on_dirs(W_eff, Wg_state, Ddirs);  % (Nd x 1) each
gap_dir = max( su_true - su_gray, 0 );                         % violation if >0
max_gap = max(gap_dir);

% --- Plotting
nrows = 2 + (nx==2);   % if nx==2 we show a polar panel, else skip it
figure('Color','w','Name','Noise matching');
t = tiledlayout(nrows,1,'Padding','compact','TileSpacing','compact');
title(t, compose_title("Noise matching", prm.Title, sprintf('  (max dir deficiency = %.3g)', max_gap)));

% (A) 2D overlay in state space
nexttile; hold on; grid on; axis equal;
plot_set_2d(W_eff, dims, [0.2 0.2 0.7], 0.25, 'True W_{eff}');
plot_set_2d(Wg_state, dims, [0.8 0.3 0.1], 0.25, 'Gray E_g Z_d');
xlabel(sprintf('x_%d',dims(1))); ylabel(sprintf('x_%d',dims(2)));
legend('Location','bestoutside'); title('State disturbance: overlay (projection)'); box on;

% (B) Per-state interval half-widths
nexttile; bar_state_halfwidths(hw_true, hw_gray, 'True W_{eff}', 'Gray E_gZ_d');
title('Per-state interval half-widths (||G||_1 by row)'); grid on;

% (C) Directional supports (only if nx==2 use polar; else line)
if nx==2
    nexttile; 
    th = atan2(Ddirs(2,:), Ddirs(1,:));     % angles
    th = [th, th(1)];                        % close loop
    rT = [su_true(:); su_true(1)];
    rG = [su_gray(:); su_gray(1)];
    polarplot(th, rT, '-','LineWidth',1.4); hold on;
    polarplot(th, rG, '-','LineWidth',1.4);
    legend('True','Gray','Location','best'); title('Directional support (state space)');
else
    nexttile; hold on; grid on;
    plot(1:Nd, su_true, '-o', 'DisplayName','True','LineWidth',1.2);
    plot(1:Nd, su_gray, '-s', 'DisplayName','Gray','LineWidth',1.2);
    stem(find(gap_dir>0), su_true(gap_dir>0)-su_gray(gap_dir>0), ':', 'DisplayName','deficiency');
    xlabel('direction index'); ylabel('support radius'); legend('Location','best');
    title('Directional support (state space)');
end

% optional save
if strlength(prm.SaveBase)>0
    saveas(gcf, prm.SaveBase + "_noise_match.png");
    try, exportgraphics(gcf, prm.SaveBase + "_noise_match.pdf",'ContentType','vector'); catch, end
end
end

% ================= helpers =================
function dims = pick_dims_for_plot(Z)
    G = generators(zonotope(Z));
    if isempty(G), dims = [1 2]; return; end
    e = sum(G.^2,2);               % energy per state coordinate
    [~,I] = sort(e,'descend');
    dims = sort(I(1:min(2,numel(I))));
    if numel(dims)==1, dims = [dims 1+mod(dims, max(2,numel(e)))]; end
end

function [hw_true, hw_gray] = coord_halfwidths(W_true, W_gray)
    Gt = generators(zonotope(W_true));
    if isempty(Gt), Gt = zeros(size(center(W_true),1),0); end
    hw_true = sum(abs(Gt),2);
    if nargin<2 || isempty(W_gray)
        hw_gray = nan(size(hw_true));
    else
        Gg = generators(zonotope(W_gray));
        if isempty(Gg), Gg = zeros(size(center(W_gray),1),0); end
        hw_gray = sum(abs(Gg),2);
    end
end

function bar_state_halfwidths(a, b, labA, labB)
    if ~isvector(a)
        a = a(:);
    end
    if isempty(b)
        b = nan(size(a));
    else
        b = b(:);
    end
    X = [(1:numel(a))' (1:numel(a))'];
    Y = [a b];
    h = bar(X, Y, 1.0, 'grouped'); %#ok<NASGU>
    xlabel('state index i'); ylabel('half-width (L1)'); 
    legend(labA, labB, 'Location','best');
end

function [su_true, su_gray] = support_on_dirs(W_true, W_gray, Ddirs)
    Zt = zonotope(W_true);
    ct = center(Zt); Gt = generators(Zt);
    su_true = Ddirs' * ct + sum(abs(Ddirs' * Gt), 2);
    Zg = zonotope(W_gray);
    cg = center(Zg); Gg = generators(Zg);
    su_gray = Ddirs' * cg + sum(abs(Ddirs' * Gg), 2);
end

function D = sample_unit_dirs(n, Nd, seed)
    if nargin<3, seed = 12345; end
    rng(seed,'twister');
    X = randn(n, Nd);
    D = X ./ max(sqrt(sum(X.^2,1)), eps);
end

function ttl = compose_title(base, extra, tail)
    s = string(base);
    if strlength(extra)>0, s = s + " — " + string(extra); end
    if nargin>=3 && strlength(tail)>0, s = s + string(tail); end
    ttl = s;
end

function plot_set_2d(S, dims, col, alpha, name)
    % Project set via linearMap onto dims
    P = zeros(2, size(center(zonotope(S)),1)); P(1,dims(1))=1; P(2,dims(2))=1;
    Y = linearMap(zonotope(S), P);
    try
        h = plot(Y, [1 2]); 
        if ~isempty(h) && isgraphics(h(1)), set(h(1),'FaceColor',col,'FaceAlpha',alpha,'EdgeColor','none'); end
    catch
        % fallback: outline via sampled support
        th = linspace(0,2*pi,128);
        D  = [cos(th); sin(th)];
        c  = center(zonotope(Y)); G = generators(zonotope(Y));
        s  = D' * c + sum(abs(D' * G), 2);
        Pts= (D .* s')';
        K  = convhull(Pts(:,1), Pts(:,2));
        patch(Pts(K,1), Pts(K,2), col, 'FaceAlpha',alpha, 'EdgeColor','none');
    end
    if nargin>=5, plot(NaN,NaN,'s','MarkerFaceColor',col,'MarkerEdgeColor',col,'DisplayName',name); end
end
