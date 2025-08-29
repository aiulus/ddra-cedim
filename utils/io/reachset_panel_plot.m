function reachset_panel_plot(artfile, varargin)
% PLOT_FOURWAY_SINGLEPANEL
% Overlay in ONE figure:
%   (1) nominal validation trajectory (n_m_val must be 1)
%   (2) TRUE reachable output sets (true model + true uncertainty)
%   (3) RCSI/Gray reachable output sets
%   (4) DDRA reachable output sets (from M_AB)
%
% Usage:
%   plot_fourway_singlepanel('.../artifacts/row_0001.mat', ...
%       'Dims',[1 2], 'Every',1, 'Order',50, 'Title','Row 1 overlay');

    p = inputParser;
    addParameter(p,'Dims',[1 2]);
    addParameter(p,'Every',1);
    addParameter(p,'Order',50);
    addParameter(p,'Title','Validation overlay');
    parse(p,varargin{:});
    dims  = p.Results.Dims;
    Every = p.Results.Every;
    Kcap  = p.Results.Order;
    ttl   = p.Results.Title;

    S = load(artfile);   % expects: sys_gray, sys_ddra, VAL, W_eff, M_AB, meta
    sys_true = S.sys_ddra;       % CORA linearSysDT (true/"plant" model)
    sys_gray = S.sys_gray;       % CORA linearSysDT (identified)
    VAL      = S.VAL;            % struct with R0, U, y, (optionally u)
    W_eff    = S.W_eff;          % zonotope (may be zero)
    M_AB     = S.M_AB;           % matrix zonotope

    % Basic dims and caps
    nk = size(VAL.y,1);                  % time steps available on VAL
    if nk == 0, error('VAL.y is empty'); end
    if Every < 1, Every = 1; end

    % ---- Build TRUE and GRAY params (augment inputs with W if needed) ----
    Uset = VAL.U; R0 = VAL.R0;

    % Try to augment the true system with W (if the model supports it)
    U_true = Uset;  sys_true_aug = sys_true;
    try
        if exist('W_eff','var') && isa(W_eff,'zonotope') && ~isempty(generators(W_eff))
            U_true       = cartProd(Uset, W_eff);
            sys_true_aug = augment_u_with_w(sys_true);
        end
    catch
        % fall back: run without explicit augmentation
        U_true = Uset;
        sys_true_aug = sys_true;
    end

    % Match gray model input dimension
    U_gray = Uset;
    try
        nu_gray = sys_gray.nrOfInputs;
        if exist('W_eff','var') && isa(W_eff,'zonotope') && ~isempty(generators(W_eff))
            if nu_gray == dim(Uset) + dim(W_eff)
                U_gray = cartProd(Uset, W_eff);
            end
        end
    catch
        U_gray = Uset;
    end

    % CORA reach options (keep light + consistent)
    options_reach = struct( ...
        'zonotopeOrder',     Kcap, ...
        'tensorOrder',       2, ...
        'errorOrder',        1, ...
        'tensorOrderOutput', 2, ...
        'verbose',           false);

    % ---- TRUE reachable OUTPUT sets via CORA ----
    params_true = struct('R0', R0, 'U', U_true);
    R_true = reach(sys_true_aug, params_true, options_reach);
    Y_true = R_true.timePoint.set;    % cell {k} of output sets

    % ---- GRAY reachable OUTPUT sets via CORA ----
    params_gray = struct('R0', R0, 'U', U_gray);
    R_gray = reach(sys_gray, params_gray, options_reach);
    Y_gray = R_gray.timePoint.set;

    % ---- DDRA OUTPUT sets: propagate with M_AB, then map to y ----
    Xk = R0; Xsets = cell(nk,1);
    for k = 1:nk
        Xk    = reduce(Xk,'girard',Kcap);
        Xnext = M_AB * cartProd(Xk, Uset) + W_eff;
        Xnext = reduce(Xnext,'girard',Kcap);
        Xsets{k} = Xnext;
        Xk = Xnext;
    end
    Y_ddra = cell(nk,1);
    for k = 1:nk
        % y_k = C x_k (+ D u_k). Use full input set at step k for consistency.
        Y_ddra{k} = linearMap(sys_true.C, Xsets{k});
        if isprop(sys_true,'D') || isfield(sys_true,'D')
            try
                if any(sys_true.D(:) ~= 0)
                    Y_ddra{k} = Y_ddra{k} + linearMap(sys_true.D, Uset);
                end
            catch
                % ignore if dimensions don't align
            end
        end
    end

    % ---- Nominal validation trajectory (must be single trajectory) ----
    y_nom = squeeze(VAL.y(:, :, 1));   % [nk x ny]
    if size(y_nom,1) ~= nk
        nk = size(y_nom,1);
    end

    % ---- Plot (single panel) ----
    figure('Color','w','Name','Four-way overlay'); hold on;

    % To avoid legend spam, only label the first poly of each family
    didTrue=false; didGray=false; didDDRA=false;

    for k = 1:Every:nk
        % TRUE
        if iscell(Y_true) && k <= numel(Y_true) && isa_wrap(Y_true{k})
            lbl = iff(~didTrue,'True reachable','');
            pset(Y_true{k}, dims, [0.85 0.85 0.85], 0.6, '-', 0.5, lbl);
            didTrue = true;
        end
        % GRAY
        if iscell(Y_gray) && k <= numel(Y_gray) && isa_wrap(Y_gray{k})
            lbl = iff(~didGray,'RCSI/Gray','');
            pset(Y_gray{k}, dims, [0.90 0.45 0.10], 0.25, '--', 1.0, lbl);
            didGray = true;
        end
        % DDRA
        if iscell(Y_ddra) && k <= numel(Y_ddra) && isa_wrap(Y_ddra{k})
            lbl = iff(~didDDRA,'DDRA','');
            pset(Y_ddra{k}, dims, [0.20 0.45 0.80], 0.25, '-', 1.2, lbl);
            didDDRA = true;
        end
    end

    % Nominal points
    plot(y_nom(:,dims(1)), y_nom(:,dims(2)), 'kx-', 'DisplayName','Nominal y');

    grid on; axis tight; axis equal;
    xlabel(sprintf('y_{%d}', dims(1))); ylabel(sprintf('y_{%d}', dims(2)));
    title(ttl);
    legend('Location','best');

    % --------- helpers ----------
    function tf = isa_wrap(S)
        tf = ~isempty(S) && (isa(S,'contSet') || isfield(S,'generators') || isa(S,'zonotope'));
    end
    function pset(Z, d, face, alpha, ls, lw, name)
        % plot a CORA set (zonotope/contSet) on dims d
        try
            if isa(Z,'contSet')
                plot(Z, d, 'FaceColor', face, 'FaceAlpha', alpha, ...
                    'LineStyle', ls, 'LineWidth', lw, ...
                    iff(~isempty(name),'DisplayName',name, 'HandleVisibility','off'));
            else
                Zz = zonotope(Z); 
                plot(zonotope(Z), d, 'FaceColor', face, 'FaceAlpha', alpha, ...
                    'LineStyle', ls, 'LineWidth', lw, ...
                    iff(~isempty(name),'DisplayName',name, 'HandleVisibility','off'));
            end
        catch
            % best effort: ignore if plotting fails
        end
    end
    function out = iff(cond, a, b)
        if cond, out = a; else, out = b; end
    end
end
