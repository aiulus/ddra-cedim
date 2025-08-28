function [ctrain_gray, cval_gray, Tvalidate_g] = gray_containment_on_VAL(configs, TS_train, TS_val, tol)
%GRAY_CONTAINMENT_ON_VAL  Containment for Gray on EXACT (x0,u) points.
% Mirrors validateReach, but:
%   - strips options.cs before reach()
%   - uses contains_interval() for membership (same as DDRA)
%
% Inputs
%   configs   : {true, gray, ...} as from gray_identify()
%   TS_train  : cell of structs with fields: initialState (nx×1), u (m×n_k×1)
%   TS_val    : same as TS_train (validation)
%   tol       : optional interval slack (default 1e-6)
%
% Outputs
%   ctrain_gray : % containment on training points (gray config only)
%   cval_gray   : % containment on validation points (gray config only)
%   Tvalidate_g : runtime (s) spent in validation loop

    if nargin < 4 || isempty(tol), tol = 1e-6; end

    % pick the *first* identified gray config (skip "true")
    if numel(configs) >= 2
        i_gray = 2;
    else
        i_gray = 1;  % fallback if only one config exists
    end
    cfg = configs{i_gray};

    % strip cs from options to match DDRA reach behavior
    if isfield(cfg, 'options') && isfield(cfg.options, 'cs')
        options_no_cs = rmfield(cfg.options, 'cs');
    else
        options_no_cs = cfg.options;
    end

    % ---- TRAIN ----
    num_in  = 0; num_all = 0;
    for m = 1:numel(TS_train)
        params_m = cfg.params;
        params_m.R0 = params_m.R0 + TS_train{m}.initialState;
        params_m.u  = TS_train{m}.u';                           % (m × n_k)
        params_m.tFinal = cfg.sys.dt * (size(params_m.u,2)-1);

        Rg = reach(cfg.sys, params_m, options_no_cs);
        Rg = Rg.timePoint.set;                                  % {1..n_k}

        hasY = isfield(TS_train{m}, 'y') && ~isempty(TS_train{m}.y);
        if ~hasY, continue; end

        Y = TS_train{m}.y;                                      % (n_k × dim_y × s)
        for k = 1:size(Y,1)
            % y_k: dim_y × s
            yk = squeeze(Y(k,:,:));
            if isvector(yk), yk = yk(:); end
            for s = 1:size(yk,2)
                if contains_interval(yk(:,s), Rg{k}, tol)
                    num_in = num_in + 1;
                end
                num_all = num_all + 1;
            end
        end
    end
    if num_all > 0
        ctrain_gray = 100 * (num_in / num_all);
    else
        ctrain_gray = NaN;  % no train measurements available
    end

    % ---- VAL ----
    t0 = tic;
    num_in  = 0; num_all = 0;
    for m = 1:numel(TS_val)
        params_m = cfg.params;
        params_m.R0 = params_m.R0 + TS_val{m}.initialState;
        params_m.u  = TS_val{m}.u';
        params_m.tFinal = cfg.sys.dt * (size(params_m.u,2)-1);

        Rg = reach(cfg.sys, params_m, options_no_cs);
        Rg = Rg.timePoint.set;

        hasY = isfield(TS_val{m}, 'y') && ~isempty(TS_val{m}.y);
        if ~hasY, continue; end

        Y = TS_val{m}.y;
        for k = 1:size(Y,1)
            yk = squeeze(Y(k,:,:));
            if isvector(yk), yk = yk(:); end
            for s = 1:size(yk,2)
                if contains_interval(yk(:,s), Rg{k}, tol)
                    num_in = num_in + 1;
                end
                num_all = num_all + 1;
            end
        end
    end
    Tvalidate_g = toc(t0);

    if num_all > 0
        cval_gray = 100 * (num_in / num_all);
    else
        cval_gray = NaN;
    end
end

% ----- shared interval-hull membership -----
function tf = contains_interval(y, Z, tol)
% Accept if y lies in the interval hull of set Z (zonotope/contSet), with slack tol.
% Works for CORA sets: if not a contSet, coerce to zonotope first.
    if nargin < 3 || isempty(tol), tol = 1e-6; end
    if ~isa(Z,'contSet'), Z = zonotope(Z); end
    Iv = interval(Z);                            % Iv.inf, Iv.sup
    y  = y(:);                                   % column
    tf = all(y >= Iv.inf - tol & y <= Iv.sup + tol);
end
