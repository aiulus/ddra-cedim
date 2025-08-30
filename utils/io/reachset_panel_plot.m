function reachset_panel_plot(artfile, varargin)
% REACHSET_PANEL_PLOT
% Overlay in ONE figure:
%   (1) nominal validation trajectory (n_m_val must be 1 or pick a Block)
%   (2) TRUE reachable output sets (true model + true uncertainty)
%   (3) RCSI/Gray reachable output sets
%   (4) DDRA reachable output sets (from M_AB)
%
% Usage:
%   reachset_panel_plot('.../artifacts/row_0001.mat', ...
%       'Block',1,'Dims',[1 2],'Every',1,'Order',50, ...
%       'Title','True vs RCSI/Gray vs DDRA (VAL)');

    % ---- parse name–value args ----
    p = inputParser;
    addParameter(p,'Block',1);
    addParameter(p,'Dims',[1 2]);
    addParameter(p,'Every',1);
    addParameter(p,'Order',50);
    addParameter(p,'Title','Validation overlay');
    parse(p,varargin{:});
    blk   = p.Results.Block;
    dims  = p.Results.Dims;
    Every = max(1, round(p.Results.Every));
    Kcap  = p.Results.Order;
    ttl   = p.Results.Title;

    % ---- load artifact (supports old/new schema) ----
    S = load(artfile);  % expects some of: sys_gray, sys_ddra, VAL, W_eff, M_AB, meta
    mustHave(S,'sys_ddra');   mustHave(S,'sys_gray');   mustHave(S,'VAL');

    % systems (tolerate struct {A,B,C,D,dt})
    sys_true = toLinearSysDT(S.sys_ddra);
    sys_gray = toLinearSysDT(S.sys_gray);

    % --- pull nominal y,u for chosen block (first sample if 3rd dim present)
    assert(iscell(S.VAL.y) && numel(S.VAL.y)>=blk, 'VAL.y has no block %d.', blk);
    y_nom = squeezeFirst(S.VAL.y{blk});          % [n_k x n_y]
    u_nom = fetchField(S.VAL,'u',blk,[]); 
    u_nom = squeezeFirst(u_nom);                  % [n_k x n_u] or []
    nk    = size(y_nom,1);                        % number of time steps
    if nk < 1, error('Empty validation block.'); end

    % --- R0 (fallback for old artifacts: point set at nominal x0)
    if isfield(S.VAL,'R0')
        R0 = S.VAL.R0;
    else
        assert(isfield(S.VAL,'x0') && numel(S.VAL.x0)>=blk, 'VAL.x0 missing.');
        x0_nom = squeezeFirst(S.VAL.x0{blk});
        if isrow(x0_nom), x0_nom = x0_nom.'; end
        R0 = zonotope(x0_nom);
    end

    % --- U (fallback to zero-uncertainty input set)
    if isfield(S.VAL,'U')
        Uset = S.VAL.U;
    else
        nu = size(sys_true.B,2);
        Uset = zonotope(zeros(nu,1));
    end

    % --- disturbance set (effective)
    if isfield(S,'W_eff')
        W = S.W_eff;
    elseif isfield(S,'W')
        W = S.W;
    else
        W = zonotope(0); % zero disturbance (dim inferred via augmentation)
    end

    % DDRA M_AB may be absent in old artifacts
    haveMAB = isfield(S,'M_AB');
    if haveMAB
        M_AB = S.M_AB;
    end

    % ---- build reach options ----
    options_reach = struct( ...
        'zonotopeOrder',     Kcap, ...
        'tensorOrder',       2, ...
        'errorOrder',        1, ...
        'tensorOrderOutput', 2, ...
        'verbose',           false);

    % ---- TRUE / GRAY input augmentation with W (only if dims match) ----
    U_true = maybeAugmentU(sys_true, Uset, W);
    U_gray = maybeAugmentU(sys_gray, Uset, W);

    % ---- nominal input time series (transpose to [nu x nk] and pad zeros)
    u_true_nom = nominalUFor(sys_true, u_nom);
    u_gray_nom = nominalUFor(sys_gray, u_nom);

    % ---- tFinal for CORA (discrete-time): dt*(nk-1) ----
    tFinal_true = sys_true.dt * (nk - 1);
    tFinal_gray = sys_gray.dt * (nk - 1);

    % ---- TRUE reachable OUTPUT sets via CORA ----
    params_true = struct('R0', R0, 'U', U_true, 'u', u_true_nom, 'tFinal', tFinal_true);
    R_true = reach(sys_true, params_true, options_reach);
    Y_true = R_true.timePoint.set;          % cell {k}

    % ---- GRAY reachable OUTPUT sets via CORA ----
    params_gray = struct('R0', R0, 'U', U_gray, 'u', u_gray_nom, 'tFinal', tFinal_gray);
    R_gray = reach(sys_gray, params_gray, options_reach);
    Y_gray = R_gray.timePoint.set;

    % ---- DDRA OUTPUT sets via propagation with M_AB (x-space) ----
    if haveMAB
        Xk = R0;  Xsets = cell(nk,1);
        for k = 1:nk
            Xk    = reduce(Xk,'girard',Kcap);
            Xnext = M_AB * cartProd(Xk, Uset) + W;
            Xnext = reduce(Xnext,'girard',Kcap);
            Xsets{k} = Xnext;
            Xk = Xnext;
        end
        % map to outputs y_k = C x_k (+ D u_k)
        Y_ddra = cell(nk,1);
        for k = 1:nk
            Y_ddra{k} = sys_true.C * Xsets{k};
            if hasD(sys_true)
                try
                    Y_ddra{k} = Y_ddra{k} + linearMap(sys_true.D, Uset);
                catch
                    % ignore if dims do not align
                end
            end
        end
    else
        Y_ddra = {};
    end

    % ---- Plot (single panel) ----
    figure('Color','w','Name','True vs RCSI/Gray vs DDRA (VAL)'); hold on;

    didTrue=false; didGray=false; didDDRA=false;
    for k = 1:Every:nk
        % TRUE
        if iscell(Y_true) && k <= numel(Y_true) && isSet(Y_true{k})
            label = iff(~didTrue,'True reachable','');
            plotSet(Y_true{k}, dims, [0.85 0.85 0.85], 0.6, '-', 0.5, label);
            didTrue = true;
        end
        % GRAY
        if iscell(Y_gray) && k <= numel(Y_gray) && isSet(Y_gray{k})
            label = iff(~didGray,'RCSI/Gray','');
            plotSet(Y_gray{k}, dims, [0.90 0.45 0.10], 0.25, '--', 1.0, label);
            didGray = true;
        end
        % DDRA (optional)
        if haveMAB && iscell(Y_ddra) && k <= numel(Y_ddra) && isSet(Y_ddra{k})
            label = iff(~didDDRA,'DDRA','');
            plotSet(Y_ddra{k}, dims, [0.20 0.45 0.80], 0.25, '-', 1.2, label);
            didDDRA = true;
        end
    end

    % nominal validation trajectory (first sample)
    plot(y_nom(:,dims(1)), y_nom(:,dims(2)), 'kx-', 'DisplayName','Nominal y');

    grid on; axis tight; axis equal;
    xlabel(sprintf('y_{%d}', dims(1))); ylabel(sprintf('y_{%d}', dims(2)));
    title(ttl);
    legend('Location','best');

    % =================== helpers ===================
    function mustHave(S, f)
        assert(isfield(S,f), 'Artifact missing required field "%s".', f);
    end

    function sysL = toLinearSysDT(sysAny)
        if isa(sysAny,'linearSysDT')
            sysL = sysAny; return
        end
        % tolerate plain struct {A,B,C,D,dt}
        if isstruct(sysAny) && all(isfield(sysAny,{'A','B','C','D','dt'}))
            sysL = linearSysDT(sysAny.A, sysAny.B, sysAny.C, sysAny.D, sysAny.dt);
            return
        end
        sysL = sysAny; % hope downstream works if already CORA type
    end

    function Z = fetchField(VAL, name, idx, defaultZ)
        if isfield(VAL, name) && iscell(VAL.(name)) && numel(VAL.(name))>=idx
            Z = VAL.(name){idx};
        else
            Z = defaultZ;
        end
    end

    function M = squeezeFirst(Min)
        if isempty(Min), M = Min; return, end
        if ndims(Min) == 3 && size(Min,3) > 1
            M = Min(:,:,1);
        else
            M = Min;
        end
    end

    function tf = isSet(Z)
        tf = ~isempty(Z) && (isa(Z,'contSet') || isa(Z,'zonotope') || isfield(Z,'generators'));
    end

    function out = iff(cond, a, b), if cond, out=a; else, out=b; end, end

    function tf = hasD(sysL)
        tf = (isprop(sysL,'D') || isfield(sysL,'D')) && ~isempty(sysL.D) && any(sysL.D(:)~=0);
    end

    function U_aug = maybeAugmentU(sysL, U0, Wz)
        U_aug = U0;
        % try input counts; if I/O mismatch, leave un-augmented
        try
            nu_sys = sysL.nrOfInputs;
        catch
            nu_sys = size(sysL.B,2);
        end
        nu0 = safeDim(U0);
        nw  = safeDim(Wz);
        if nu_sys == nu0 + nw && nw > 0
            U_aug = cartProd(U0, Wz);
        end
    end

    function n = safeDim(Z)
        try
            n = dim(Z);
        catch
            try, n = size(center(Z),1); catch, n = 0; end
        end
    end

    function u_nom_out = nominalUFor(sysL, u_nom_in)
        % Make params.u of size [nu x nk], padding zeros if needed
        if isempty(u_nom_in)
            nk_loc = nk;
            try, nu_sys = sysL.nrOfInputs; catch, nu_sys = size(sysL.B,2); end
            u_nom_out = zeros(nu_sys, nk_loc);
            return
        end
        u_nom_T = u_nom_in.';   % [n_u x n_k]
        try
            nu_sys = sysL.nrOfInputs;
        catch
            nu_sys = size(sysL.B,2);
        end
        diff_u = nu_sys - size(u_nom_T,1);
        if diff_u > 0
            u_nom_out = cat(1, u_nom_T, zeros(diff_u, size(u_nom_T,2)));
        elseif diff_u < 0
            u_nom_out = u_nom_T(1:nu_sys,:); % truncate if too many (defensive)
        else
            u_nom_out = u_nom_T;
        end
    end

    function plotSet(Z, d, face, alpha, ls, lw, name)
        try
            if ~isa(Z,'contSet'), Z = zonotope(Z); end
            if isempty(name)
                plot(Z, d, 'FaceColor', face, 'FaceAlpha', alpha, ...
                    'LineStyle', ls, 'LineWidth', lw, 'HandleVisibility','off');
            else
                plot(Z, d, 'FaceColor', face, 'FaceAlpha', alpha, ...
                    'LineStyle', ls, 'LineWidth', lw, 'DisplayName', name);
            end
        catch
            % ignore plotting failures (dim mismatch etc.)
        end
    end
end
