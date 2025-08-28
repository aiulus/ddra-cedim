function [ctrain_gray, cval_gray, Tvalidate_g] = gray_containment_on_VAL(configs, TS_train, TS_val, tol)
    if nargin < 4 || isempty(tol), tol = 1e-6; end

    % pick gray config
    i_gray = min(2, numel(configs));  % 2 if available, else 1
    cfg = configs{i_gray};

    % strip cs from options
    if isfield(cfg,'options') && isfield(cfg.options,'cs')
        options_no_cs = rmfield(cfg.options,'cs');
    else
        options_no_cs = cfg.options;
    end

    % normalize inputs to cell arrays of testCase
    TS_train = ensure_cell_tc(TS_train);
    TS_val   = ensure_cell_tc(TS_val);

    % ---------------- TRAIN ----------------
    num_in = 0; num_all = 0;
    for m = 1:numel(TS_train)
        tc = TS_train{m};

        % reach on this (x0,u)
        params_m = cfg.params;
        params_m.R0 = params_m.R0 + tc.initialState;
        params_m.u  = tc.u';                           % (m × n_k)
        params_m.tFinal = cfg.sys.dt * (size(params_m.u,2)-1);
        Rg = reach(cfg.sys, params_m, options_no_cs);
        Rg = Rg.timePoint.set;                         % {1..K}

        % ensure we have y; if not, simulate
        Y = ensure_y(tc, cfg.sys);
        K = min(size(Y,1), numel(Rg));                 % guard against mismatch

        for k = 1:K
            yk = squeeze(Y(k,:,:));                    % (ny × s)
            if isvector(yk), yk = yk(:); end
            for s = 1:size(yk,2)
                if contains_interval(yk(:,s), Rg{k}, tol)
                    num_in = num_in + 1;
                end
                num_all = num_all + 1;
            end
        end
    end
    ctrain_gray = pct_or_nan(num_in, num_all);

    % ---------------- VAL ----------------
    t0 = tic;
    num_in = 0; num_all = 0;
    for m = 1:numel(TS_val)
        tc = TS_val{m};

        params_m = cfg.params;
        params_m.R0 = params_m.R0 + tc.initialState;
        params_m.u  = tc.u';
        params_m.tFinal = cfg.sys.dt * (size(params_m.u,2)-1);
        Rg = reach(cfg.sys, params_m, options_no_cs);
        Rg = Rg.timePoint.set;

        Y = ensure_y(tc, cfg.sys);
        K = min(size(Y,1), numel(Rg));

        for k = 1:K
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
    cval_gray = pct_or_nan(num_in, num_all);
end

% ---- helpers ----
function TS = ensure_cell_tc(TS)
    if isempty(TS), TS = {}; end
    if ~iscell(TS), TS = {TS}; end
end

function Y = ensure_y(tc, sys)
    Y = [];
    if isfield(tc,'y') && ~isempty(tc.y)
        Y = tc.y;
    end
    if isempty(Y) || ~ismatrix(Y(:,:,1))
        % simulate outputs to fill Y; time-major (a × ny × s)
        params_sim.tFinal = sys.dt * (size(tc.u,1)-1);
        params_sim.x0 = tc.initialState;
        params_sim.u  = tc.u';
        [~,~,~,Ysim] = simulate(sys, params_sim);
        Y = reshape(Ysim, size(tc.u,1), sys.nrOfOutputs, 1);
    end
end

function p = pct_or_nan(in, all)
    if all>0, p = 100*in/all; else, p = NaN; end
end
