function plot_reach_all_onepanel(sys_true, sys_gray, VAL, varargin)
% PLOT_REACH_ALL_ONEPANEL

    % ---------- parse options ----------
    %opt = struct('Block',1,'Dims',[1 2],'MAB',[],'W',[], ...
    %             'Order',100,'Reduce','girard','ShowSamples',true,'Save','');
    opt = struct('Block',1,'Dims',[1 2],'MAB',[],'W',[], ...
             'Order',100,'Reduce','girard','ShowSamples',true,'Save','', ...
             'ReachOptions',[]); 

    for i=1:2:numel(varargin)
        if isfield(opt,varargin{i}), opt.(varargin{i}) = varargin{i+1}; end
    end
    b = opt.Block; dims = opt.Dims(:)'; assert(numel(dims)==2,'Dims must be 2-vector.');

    % ---------- shapes ----------
    assert(isa(sys_true,'linearSysDT') && isa(sys_gray,'linearSysDT'), 'Use CORA linearSysDT systems.');
    assert(abs(sys_true.dt - sys_gray.dt) < 1e-12, 'dt mismatch between true and gray systems.');

    nx = size(sys_true.A,1);
    m  = size(sys_true.B,2);

    % C/D safe defaults
    Ctrue = safe_mat(sys_true,'C',eye(nx));
    Dtrue = safe_mat(sys_true,'D',zeros(size(Ctrue,1), m));
    Cgray = safe_mat(sys_gray,'C',Ctrue);
    Dgray = safe_mat(sys_gray,'D',Dtrue);
    ny    = size(Ctrue,1);

    assert(all(dims>=1 & dims<=ny), 'Dims out of range for outputs.');

    % ---------- x0, u, y ----------
    x0  = VAL.x0{b};  if isrow(x0), x0 = x0.'; end
    U   = normalize_u(VAL.u{b}, m);         % m × n_k
    n_k = size(U,2);

    % samples (optional) — make Ymeas always exist
    Ymeas = [];
    haveY = isfield(VAL,'y') && numel(VAL.y)>=b && ~isempty(VAL.y{b});
    if haveY
        Ymeas = VAL.y{b};                   % (n_k×ny) or (ny×n_k)
        if size(Ymeas,2)==ny && size(Ymeas,1)==n_k
            % ok
        elseif size(Ymeas,1)==ny && size(Ymeas,2)==n_k
            Ymeas = Ymeas.';                % (n_k×ny)
        else
            haveY = false;
            Ymeas = [];
        end
    end
    Ymeas_for_limits = Ymeas;               % [] when haveY==false

    % ---------- sets & options ----------
    hasR0  = isfield(VAL,'R0') && isa(VAL.R0,'zonotope') && ~isempty(VAL.R0);
    if hasR0
        R0_set = VAL.R0;
    else
        R0_set = zonotope(zeros(nx,1));
    end
      
      
    % ---------- disturbances (prefer opt.W; else VAL.W; else empty) ----------
    W_in = [];
    if ~isempty(opt.W)
        W_in = opt.W;
    elseif isfield(VAL,'W') && ~isempty(VAL.W)
        W_in = VAL.W;
    end

    % ---------- reach options (CORA) ----------
    % Build a minimal options struct; allow caller to override via opt.ReachOptions.
    options = struct('zonotopeOrder', opt.Order, 'reductionTechnique', opt.Reduce);
    if isfield(opt,'ReachOptions') && ~isempty(opt.ReachOptions)
        fns = fieldnames(opt.ReachOptions);
        for ii = 1:numel(fns)
            options.(fns{ii}) = opt.ReachOptions.(fns{ii});
        end
    end

    % ---------- TRUE reach ----------
    paramsT = struct('R0', R0_set + x0, 'u', U, ...
                     'tStart', 0, 'tFinal', sys_true.dt*(n_k-1));
    if sys_true.nrOfDisturbances > 0 && ~isempty(W_in)
        paramsT.W = map_W_for_sys(sys_true, W_in);   % map to disturbance space
    end
    Rtrue = reach(sys_true, paramsT, options);

    Xtrue = Rtrue.timePoint.set;
    Ytrue = map_sets_to_output(Xtrue, Ctrue, Dtrue, U);

    % ---------- GRAY reach ----------
    paramsG = struct('R0', R0_set + x0, 'u', U, ...
                     'tStart', 0, 'tFinal', sys_gray.dt*(n_k-1));
    if sys_gray.nrOfDisturbances > 0 && ~isempty(W_in)
        try
            % map state-disturbance W into gray's disturbance space (conservative)
            paramsG.W = normalizeWForGray(sys_gray, W_in);
        catch
            paramsG.W = W_in; % safe fallback if helper not available
        end
    end
    Rgray = reach(sys_gray, paramsG, options);
    Xgray = Rgray.timePoint.set;
    Ygray = map_sets_to_output(Xgray, Cgray, Dgray, U);

    % ---------- DDRA (optional) ----------
    Yddra = {};
    if ~isempty(opt.MAB)
        Xdd = cell(1,n_k); Xdd{1} = R0_set + x0;
        for k = 1:n_k-1
            u_k = U(:,k);
            Xk  = reduce(toZono(Xdd{k}), opt.Reduce, opt.Order);
            Zk  = cartProd(Xk, zonotope(u_k));
            Xn  = opt.MAB * Zk + opt.W;
            Xdd{k+1} = reduce(Xn, opt.Reduce, opt.Order);
        end
        Yddra = map_sets_to_output(Xdd, Ctrue, Dtrue, U);
    end

    % ---------- nominal trajectory ----------
    xt = zeros(nx, n_k); xt(:,1) = x0;
    for k = 1:n_k-1, xt(:,k+1) = sys_true.A*xt(:,k) + sys_true.B*U(:,k); end
    yt = Ctrue*xt + Dtrue*U;                            % ny×n_k

    % ---------- figure and colors (define Cc BEFORE using it) ----------
    Cc = struct(...
        'true',[0.93 0.49 0.19], ...   % orange (for outline + hatch)
        'ddra',[0.55 0.55 0.55], ...   % gray
        'gray',[0.00 0.40 0.74], ...   % TUM blue
        'nom',[0 0 0], ...
        'meas',[0.10 0.10 0.10]);

    figure('Color','w'); hold on; box on; grid on;
    title('Reachable output sets: TRUE vs DDRA vs GRAY'); 
    xlabel(sprintf('y_{%d}',dims(1))); 
    ylabel(sprintf('y_{%d}',dims(2)));

    % 1) Filled sets FIRST (so True can sit on top)
    if ~isempty(Yddra)
        [hD,~] = plot_series(Yddra, dims, Cc.ddra, 0.12, 'filled');  % gray
        if ~isempty(hD) && isgraphics(hD(1)), set(hD(1),'DisplayName','DDRA reach'); end
    end
    [hG,~] = plot_series(Ygray, dims, Cc.gray, 0.15, 'filled');      % TUM blue
    if ~isempty(hG) && isgraphics(hG(1)), set(hG(1),'DisplayName','Gray reach'); end

    % 2) TRUE reach as diagonal hatch + outline (do this ON this axes)
    % (i) invisible filled patch to capture polygon vertices
    [~,polyTrueFilled] = plot_series(Ytrue, dims, Cc.true, 0.001, 'filled'); % FaceAlpha ~0
    % (ii) overlay diagonal hatch inside each polygon
    for i = 1:numel(polyTrueFilled)
        P = polyTrueFilled{i};
        if ~isempty(P), apply_diagonal_hatch(P.X, P.Y, lighten(Cc.true, 0.45)); end
    end
    % (iii) bold orange outline on top
    [hTrue,~] = plot_series(Ytrue, dims, Cc.true, 0.00, 'outline');
    if ~isempty(hTrue) && isgraphics(hTrue(1))
        set(hTrue(1),'DisplayName','True reach','LineWidth',1.4);
    end

    % 3) Nominal traj (thinner), crosses
    hNom = plot(yt(dims(1),:), yt(dims(2),:), '-', ...
        'Color', Cc.nom, 'LineWidth', 0.9, 'DisplayName','Nominal (true)');
    plot(yt(dims(1),:), yt(dims(2),:), 'x', ...
        'Color', Cc.nom, 'MarkerSize', 5, 'LineWidth', 0.9, 'HandleVisibility','off');

    % axes limits before drawing arrows (so arrows scale conversions are correct)
    lims = auto_limits({Ytrue, Ygray, Yddra}, dims, Ymeas_for_limits);
    axis(lims);

    % measurements overlay (markers only, no connecting lines)
    if haveY && opt.ShowSamples
        plot(Ymeas(:,dims(1)), Ymeas(:,dims(2)), 'o', ...
            'Color', Cc.meas, 'MarkerSize', 3, ...
            'MarkerFaceColor', [0.6 0.6 0.6], 'LineStyle','none', ...
            'HandleVisibility','off');
    end

    % triangular arrow heads between steps (consistent, figure-anchored)
    draw_annotation_arrows(gca, yt(dims(1),:), yt(dims(2),:), Cc.nom);

    % initial state marker
    hInit = plot(yt(dims(1),1), yt(dims(2),1), 'o', ...
        'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'MarkerSize', 6, ...
        'DisplayName','Initial state');

    % keep nominal on top (in case filled sets overlap)
    if isgraphics(hNom), uistack(hNom,'top'); end

    % legend (robust)
    H = gobjects(0); L = {};
    if isgraphics(hInit), H(end+1)=hInit; L{end+1}='Initial state'; end
    if isgraphics(hNom),  H(end+1)=hNom;  L{end+1}='Nominal (true)'; end
    if ~isempty(hTrue) && isgraphics(hTrue(1)), H(end+1)=hTrue(1); L{end+1}='True reach'; end
    if ~isempty(hG)    && isgraphics(hG(1)),    H(end+1)=hG(1);    L{end+1}='Gray reach'; end
    if exist('hD','var') && ~isempty(hD) && isgraphics(hD(1)), H(end+1)=hD(1); L{end+1}='DDRA reach'; end
    if ~isempty(H), legend(H,L,'Location','best'); else, legend('off'); end
    set(gca,'FontSize',12);

    % ---------- save (optional) ----------
    if ~isempty(opt.Save)
        try
            if exist('save_plot','file') == 2
                save_plot(gcf, fileparts2(opt.Save), basename2(opt.Save), 'Formats',{'png','pdf'}, 'Resolution',200);
            else
                [p,~,b] = fileparts(opt.Save); if ~isempty(p) && ~exist(p,'dir'), mkdir(p); end
                exportgraphics(gcf, fullfile(p,[b '.png']), 'Resolution',200);
                exportgraphics(gcf, fullfile(p,[b '.pdf']), 'ContentType','vector');
            end
        catch ME
            warning('plot_reach_all_onepanel:saveFailed','Could not save: %s', ME.message);
        end
    end
end


% ================= helpers (local) =================

function [h, polys] = plot_series(sets, dims, color, alpha, ~)
% DRAW a list of CORA sets; fill only if the projection has area.
% RETURNS:
%   h     : graphics handles (one per set; patch or line)
%   polys : cell array same length as sets; for area-like sets contains
%           struct('X',Xvals,'Y',Yvals), otherwise [].

if isempty(sets), h = gobjects(0); polys = {}; return; end
h = gobjects(numel(sets),1);
polys = cell(1,numel(sets));

for k = 1:numel(sets)
    Z = sets{k};
    if isempty(Z) || (isnumeric(Z) && isempty(Z)), continue; end

    if is_area_like(Z, dims)
        % 2D area -> ok to fill (FaceAlpha supported)
        hh = plot(Z, dims, 'FaceColor', color, 'EdgeColor', color, 'FaceAlpha', alpha);
        if numel(hh) >= 1, h(k) = hh(1); else, h(k) = hh; end

        % Try to capture polygon vertices for hatching later
        try
            X = get(h(k),'XData'); Y = get(h(k),'YData');
            if isnumeric(X) && isnumeric(Y) && ~isempty(X) && ~isempty(Y)
                polys{k} = struct('X', X(:).', 'Y', Y(:).');
            else
                polys{k} = [];
            end
        catch
            polys{k} = [];
        end
    else
        % point/segment -> NO FaceAlpha/FaceColor (Line object)
        try
            hh = plot(Z, dims, 'Color', color, 'LineWidth', 1.2);
            if numel(hh) >= 1, h(k) = hh(1); else, h(k) = hh; end
        catch
            % last resort: plot center marker
            c = center(Z);
            h(k) = plot(c(dims(1)), c(dims(2)), '.', 'Color', color, 'MarkerSize', 12);
        end
        polys{k} = [];   % nothing to hatch
    end
end
end


function tf = is_area_like(Z, dims)
% True if the projection onto dims has nonzero area.
tf = false;
try
    % Fast check via interval width
    I = interval(Z);
    w = I.sup - I.inf;
    w = w(dims);
    tf = all(isfinite(w)) && all(w > 1e-12);
    if ~tf
        % Fallback: check generator rank in projected dims
        G = generators(Z);        % n x g
        if ~isempty(G)
            G2 = G(dims, :);      % 2 x g
            tf = rank(G2) >= 2;   % spans 2D
        end
    end
catch
    % If anything fails, be conservative: do NOT fill
    tf = false;
end
end

function draw_annotation_arrows(ax, xs, ys, color)
% Consistent triangular arrowheads using figure-normalized coordinates.
    if numel(xs) < 2, return; end
    fig = ancestor(ax,'figure');
    % axis limits → axis→figure normalization
    xl = xlim(ax); yl = ylim(ax);
    axpos = getpixelposition(ax,true);  % in pixels, absolute
    figpos = getpixelposition(fig,true);
    % convert axes pixel rect to normalized figure coords
    axnorm = [axpos(1:2)./figpos(3:4), axpos(3:4)./figpos(3:4)];

    toNorm = @(x,y)[ axnorm(1) + (x - xl(1))/(xl(2)-xl(1))*axnorm(3), ...
                     axnorm(2) + (y - yl(1))/(yl(2)-yl(1))*axnorm(4) ];

    for k = 1:numel(xs)-1
        p = toNorm(xs(k),   ys(k));
        q = toNorm(xs(k+1), ys(k+1));
        % Draw short arrow finishing at q (trim back 10% to show head nicely)
        u = [q(1)-p(1), q(2)-p(2)];
        p2 = [q(1)-0.10*u(1), q(2)-0.10*u(2)];
        h = annotation(fig,'arrow',[p(1) p2(1)],[p(2) p2(2)], ...
            'HeadStyle','vback2','LineWidth',0.8,'Color',color);
        % Optional: uniform head size (pixels)
        try
            set(h,'HeadLength',7,'HeadWidth',7);
        end
    end
end

function apply_diagonal_hatch(X, Y, color)
% Clip diagonal lines to polygon using polyxpoly (Mapping Toolbox).
    if exist('polyxpoly','file') ~= 2, return; end
    X = X(:).'; Y = Y(:).'; mask = ~isnan(X) & ~isnan(Y);
    X = X(mask); Y = Y(mask);
    if numel(X) < 3, return; end
    minx = min(X); maxx = max(X); miny = min(Y); maxy = max(Y);
    W = max(maxx-minx, maxy-miny); if ~isfinite(W) || W<=0, return; end
    step = 0.07 * W;
    cmin = (miny - maxx) - 0.2*W; cmax = (maxy - minx) + 0.2*W;
    xline = linspace(minx-0.2*W, maxx+0.2*W, 200);
    cs   = cmin:step:cmax;
    for c = cs
        yline = xline + c;
        [xi, yi] = polyxpoly(X, Y, xline, yline);
        for j = 1:2:numel(xi)-1
            line([xi(j) xi(j+1)],[yi(j) yi(j+1)], 'Color', color, ...
                 'LineWidth', 0.7, 'HandleVisibility','off');
        end
    end
end

function c = lighten(c, f)
% Mix color with white by fraction f in [0,1] to mimic "semi-transparent" look.
    c = (1-f)*c + f*[1 1 1];
end


function U = normalize_u(Uin, m)
    Uin = squeeze(Uin);
    if ndims(Uin) > 2
        error('VAL.u should be 2D (m×n_k or n_k×m).');
    end
    if size(Uin,1) == m
        U = Uin;
    elseif size(Uin,2) == m
        U = Uin.';          % n_k×m -> m×n_k
    else
        error('VAL.u has incompatible size; expected one dimension to be m=%d.', m);
    end
end

function Y = map_sets_to_output(Xcells, C, D, U)
    if isempty(Xcells), Y = {}; return; end
    Y = cell(1, numel(Xcells));
    ny = size(C,1); nxC = size(C,2);
    for k = 1:numel(Xcells)
        if isempty(Xcells{k}), Y{k} = []; continue; end
        Z = toZono(Xcells{k});
        dimZ = size(center(Z),1);
        if ~isempty(C) && dimZ == nxC
            yc = C * center(Z);  yg = C * generators(Z);
        elseif dimZ == ny || isempty(C)
            yc = center(Z);      yg = generators(Z);
            if isempty(C), ny = dimZ; end
        else
            error('map_sets_to_output:dimMismatch', ...
                 'Set dim %d incompatible with C (%dx%d).', dimZ, size(C,1), size(C,2));
        end
        if ~isempty(D) && ~isempty(U) && size(D,2) == size(U,1)
            uk = U(:, min(k, size(U,2)));
            yc = yc + D * uk;
        end
        Y{k} = zonotope(yc, yg);
    end
end

function Z = toZono(S)
    if isa(S,'zonotope'); Z = S; else; Z = zonotope(S); end
end

function A = safe_mat(sys, field, defaultA)
    try, A = sys.(field); if isempty(A), A = defaultA; end
    catch, A = defaultA; end
end


function lims = auto_limits(groups, dims, Ymeas)
    minx = inf; maxx = -inf; miny = inf; maxy = -inf;
    for gi = 1:numel(groups)
        G = groups{gi}; if isempty(G), continue; end
        for k=1:numel(G)
            S = G{k}; if isempty(S), continue; end
            Iv = interval(project(toZono(S),dims));
            minx = min(minx, Iv.inf(1)); maxx = max(maxx, Iv.sup(1));
            miny = min(miny, Iv.inf(2)); maxy = max(maxy, Iv.sup(2));
        end
    end
    if ~isempty(Ymeas)
        minx = min(minx, min(Ymeas(:,dims(1)))); maxx = max(maxx, max(Ymeas(:,dims(1))));
        miny = min(miny, min(Ymeas(:,dims(2)))); maxy = max(maxy, max(Ymeas(:,dims(2))));
    end
    if ~isfinite(minx) || ~isfinite(miny), minx=-1; maxx=1; miny=-1; maxy=1; end
    pad = 0.08; sx = maxx-minx; sy = maxy-miny;
    lims = [minx-pad*(sx+eps), maxx+pad*(sx+eps), miny-pad*(sy+eps), maxy+pad*(sy+eps)];
end

function v = iff(cond, a, b), if cond, v=a; else, v=b; end, end
function S = rmfield_if_exist(S,f), if isfield(S,f), S = rmfield(S,f); end, end
function [p,b] = split_pathbase(pathstr), [p,~,b] = fileparts(pathstr); end
function p = fileparts2(s), [p,~] = split_pathbase(s); end
function b = basename2(s), [~,b] = split_pathbase(s); end

function Wd = map_W_for_sys(sys, Win)
% Map a state-space or mismatched W into the disturbance space expected by CORA.
% If Win is already q-dim (q = nrOfDisturbances), return it as-is.

    % Dimensions
    try
        nx = size(sys.A,1);
    catch
        nx = size(center(Win),1);
    end
    try
        q = max(0, sys.nrOfDisturbances);
    catch
        q = 0;
    end
    d = size(center(Win),1);

    if q == 0
        % System has no disturbance channels; return empty so CORA doesn't check it
        Wd = [];
        return;
    end

    if d == q
        % Already in disturbance space
        Wd = Win;
        return;
    end

    if d == nx
        % State-dimensional W: get a conservative preimage in disturbance space
        try
            % Preferred if available in your toolbox
            Wd = coerceWToSys(sys, Win);
        catch
            % Fallback: normalizeWForGray builds a conservative disturbance-set preimage
            Wd = normalizeWForGray(sys, Win);
        end
        return;
    end

    % Last resort: build a diagonal zonotope in q-dim using Win's radius
    try
        rad = sum(abs(generators(Win)), 2);
        if isempty(rad), rad = 0; end
        if isscalar(rad), rad = repmat(rad, q, 1); end
        Wd = zonotope(zeros(q,1), diag(rad(1:q)));
    catch
        % If everything fails, drop W to avoid dimension errors
        Wd = [];
    end
end
