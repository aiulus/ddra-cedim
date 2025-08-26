%% HOW TO USE (Square / Nonlinear Sample-Size Sweep)
% What it does:
%   Uses CORA's simple NARX "Square" via custom_loadDynamics(), sweeps n_m
%   (number of unique input trajectories), and produces:
%     (i) runtime panels  and  (ii) fidelity/conservatism plots.
%
% Key knobs:
%   - sweep_grid.n_m_list (e.g., 2:1:20)
%   - Keep excitation fixed: sweep_grid.pe_list = {struct('mode','randn')}
%
% Memory-efficiency toggles:
%   cfg.lowmem.gray_check_contain = false;  % only applies to GRAY (unused here)
%   cfg.lowmem.append_csv         = true;   % stream CSV row-by-row
%   cfg.lowmem.zonotopeOrder_cap  = 50;     % reduce order in DDRA reach
%
% Outputs:
%   - CSV + plots via init_io() under experiments/results/{data,plots}/<tag>_sweeps

rng(1,'twister');

% ---------- Minimal cfg ----------
cfg = struct();
cfg.io = struct('save_tag', 'Square_sample_size_sweep');

% System & shared options (nonlinear, NARX p=1)
cfg.shared = struct();
cfg.shared.dyn    = "Square";        % <- from custom_loadDynamics()
cfg.shared.type   = "rand";          % ("standard" or "rand") preset for Square
cfg.shared.p_extr = 0.3;             % extreme-input probability during data gen

% Identification + reachability options (CORA-style)
cfg.shared.options_reach = struct( ...
    'zonotopeOrder',      100, ...
    'tensorOrder',        2, ...
    'errorOrder',         1, ...
    'tensorOrderOutput',  2, ...
    'verbose',            false);

% Conformance (RCSI) base options
cfg.shared.cs_base = struct( ...
    'robustnessMargin', 1e-9, ...
    'verbose', false, ...
    'cost', "interval", ...
    'constraints', "half");

% Data budgets (fixed base; validation ties to n_m unless overridden)
cfg.shared.n_m     = 5;    % input traj count (ident sanity only; sweep overrides)
cfg.shared.n_s     = 10;   % samples per traj
cfg.shared.n_k     = 4;    % horizon (train/ID)
cfg.shared.n_m_val = cfg.shared.n_m;  % validation: # input trajectories
cfg.shared.n_s_val = cfg.shared.n_s;
cfg.shared.n_k_val = cfg.shared.n_k;

% Black-box RCSI settings (fast-ish)
cfg.black = struct();
cfg.black.methodsBlack = ["blackCGP"];   % or ["blackGP","blackCGP"]
cfg.black.approx = struct( ...
    'p', 1, ...                               % Square has p=1
    'gp_parallel', false, ...
    'gp_pop_size', 50, ...
    'gp_num_gen', 30, ...
    'gp_func_names', {{'times','plus','square'}}, ...
    'gp_max_genes', 2, ...
    'gp_max_depth', 2, ...
    'cgp_num_gen', 5, ...
    'cgp_pop_size_base', 5, ...
    'save_res', false ...
);

% DDRA-Lipschitz noise setup
cfg.ddra = struct();
cfg.ddra.eta_w   = 1;      % number of W generators
cfg.ddra.alpha_w = 0.00;   % W scale (0 = turn off process noise)

% Shared noise policy (keep apples-to-apples if not studying noise)
cfg.shared.noise_for_black = false;   % set W=0 for RCSI if false
cfg.shared.noise_for_ddra  = true;    % DDRA uses W unless false

% Label in file names (based on first black-box method)
rcsi_lbl = rcsi_label_from_cfg(cfg);
cfg.io.save_tag = sprintf('%s_%s', cfg.io.save_tag, rcsi_lbl);

% ---------- Sweep grid ----------
sweep_grid = struct();
sweep_grid.n_m_list = 2:1:20;         % main sweep
sweep_grid.n_s_list = cfg.shared.n_s; % keep fixed
sweep_grid.n_k_list = cfg.shared.n_k; % keep fixed
sweep_grid.pe_list  = {struct('mode','randn')};  % fixed excitation

% ---------- Low-memory toggles ----------
cfg.lowmem = struct();
cfg.lowmem.gray_check_contain = true;      % no gray here, but harmless
cfg.lowmem.append_csv         = true;      % stream CSV rows
cfg.lowmem.zonotopeOrder_cap  = 50;        % DDRA reduction cap

% ---------- Run ----------
SUMMARY = run_sweeps_square_black_vs_ddraLip(cfg, sweep_grid);
SUMMARY = ensure_time_totals_square(SUMMARY);  % totals for plots (local helper)

% ---------- Plots: Runtime panels (total / learn / validation / inference) ----------
x = coerce_numeric(SUMMARY.n_m);
colors = struct('ddra',[0.23 0.49 0.77],'black',[0.20 0.55 0.30]);

f = figure('Name','Square | Runtime panels','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% 1) TOTAL
nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_total,  '-o', 'Color',colors.ddra,  'LineWidth',1.6, 'DisplayName','DDRA-Lip');
plot(x, SUMMARY.t_black_total, '-s', 'Color',colors.black, 'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl]);
xlabel('n_m (input trajectories)'); ylabel('Seconds');
title(['Total runtime vs n_m  (RCSI: ' rcsi_lbl ')']); legend('Location','best');

% 2) LEARNING
nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_build,  '-o', 'Color',colors.ddra,  'LineWidth',1.6, 'DisplayName','DDRA build');
plot(x, SUMMARY.t_black_learn, '-s', 'Color',colors.black, 'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl ' learn']);
xlabel('n_m'); ylabel('Seconds');
title('Learning / Build runtime'); legend('Location','best');

% 3) VALIDATION
nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_val, '-o', 'Color',colors.ddra,  'LineWidth',1.6, 'DisplayName','DDRA validate');
plot(x, SUMMARY.t_black_val,'-s', 'Color',colors.black, 'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl ' validate']);
xlabel('n_m'); ylabel('Seconds');
title('Validation runtime'); legend('Location','best');

% 4) INFERENCE
nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_reach_avg, '-o', 'Color',colors.ddra,  'LineWidth',1.6, 'DisplayName','DDRA reach/case (avg)');
plot(x, SUMMARY.t_black_infer,    '-s', 'Color',colors.black, 'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl ' infer']);
xlabel('n_m'); ylabel('Seconds');
title('Inference runtime'); legend('Location','best');

[plots_dir, ~] = init_io(cfg);
save_plot(f, plots_dir, ['square_runtime_panels_vs_nm_' rcsi_lbl], 'Formats', {'png','pdf'}, 'Resolution', 200);

% ---------- Plots: Fidelity & Conservatism vs n_m ----------
x_nm        = coerce_numeric(SUMMARY.n_m);
cval_ddra   = coerce_numeric(SUMMARY.cval_ddra);
cval_black  = coerce_numeric(SUMMARY.cval_black);
sizeI_ddra  = coerce_numeric(SUMMARY.sizeI_ddra);
sizeI_black = coerce_numeric(SUMMARY.sizeI_black);

f2 = figure('Name','Square | Fidelity & Conservatism','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% (1) Fidelity
nexttile; hold on; grid on;
plot(x_nm, cval_ddra,  '-o', 'Color',colors.ddra,  'LineWidth',1.6, 'DisplayName','DDRA-Lip');
plot(x_nm, cval_black, '-s', 'Color',colors.black, 'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl]);
xlabel('n_m (input trajectories)'); ylabel('Containment on validation (%)');
title(['Fidelity vs n_m  (RCSI: ' rcsi_lbl ')']); legend('Location','best');

% (2) Conservatism (proxy)
nexttile; hold on; grid on;
plot(x_nm, sizeI_ddra,  '-o', 'Color',colors.ddra,  'LineWidth',1.6, 'DisplayName','DDRA-Lip');
plot(x_nm, sizeI_black, '-s', 'Color',colors.black, 'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl]);
xlabel('n_m (input trajectories)'); ylabel('Aggregated interval size (proxy)');
title(['Conservatism vs n_m  (RCSI: ' rcsi_lbl ')']); legend('Location','best');

save_plot(f2, plots_dir, ['square_fidelity_conservatism_vs_nm_' rcsi_lbl], 'Formats', {'png','pdf'}, 'Resolution', 200);
close all force


%% ====================== Local functions ======================

function SUMMARY = run_sweeps_square_black_vs_ddraLip(cfg, grid)
% Sweep over n_m for CORA's "Square" using:
%   - Black-box RCSI (conform: GP/CGP)  and
%   - DDRA-Lipschitz (Alg. 6)
%
% Writes rows to CSV and returns SUMMARY table.

    [~, results_dir] = init_io(cfg);
    csv_path = fullfile(results_dir, 'summary.csv');

    LM = getfielddef(cfg,'lowmem',struct());
    append_csv = getfielddef(LM,'append_csv',true);

    % Sweep axes
    axes = struct();
    axes.n_m = as_list(getfielddef(grid,'n_m_list', getfielddef(cfg,'shared',struct()).n_m));
    axes.n_s = as_list(getfielddef(grid,'n_s_list', getfielddef(cfg,'shared',struct()).n_s));
    axes.n_k = as_list(getfielddef(grid,'n_k_list', getfielddef(cfg,'shared',struct()).n_k));
    axes.pe  = getfielddef(grid,'pe_list', {struct('mode','randn')});

    % Header (includes black + ddra columns)
    hdr = {'n_m','n_s','n_k', ...
           'cval_black','sizeI_black', ...
           't_black_learn','t_black_val','t_black_infer', ...
           'cval_ddra','sizeI_ddra', ...
           't_ddra_build','t_ddra_val','t_ddra_reach_avg'};

    if append_csv
        fid = fopen(csv_path,'w'); fprintf(fid,'%s\n',strjoin(hdr,',')); fclose(fid);
    else
        Ntot = numel(axes.n_m)*numel(axes.n_s)*numel(axes.n_k)*numel(axes.pe);
        cells = cell(Ntot, numel(hdr));
    end

    rowi = 0;

    % ------------------------ Main sweeps ------------------------
    for n_m = axes.n_m
      for n_s = axes.n_s
        for n_k = axes.n_k
          for ip = 1:numel(axes.pe)
            pe = axes.pe{ip};

            % -------- Build system & suites --------
            shared = getfielddef(cfg,'shared',struct());
            dyn    = getfielddef(shared,'dyn',"Square");
            type   = getfielddef(shared,'type',"standard");
            p_extr = getfielddef(shared,'p_extr',0.3);

            [sys, R0, U, ~] = custom_loadDynamics(dyn, type);  %# <- required
            options_reach = getfielddef(shared,'options_reach',struct());

            % True params and suites
            params_true = struct();
            params_true.R0     = R0;
            params_true.U      = U;
            params_true.tFinal = sys.dt * (n_k - 1);

            ts_opts = struct('p_extr', p_extr);
            params_true.testSuite       = createTestSuite(sys, params_true, n_k, n_m, n_s, ts_opts);

            n_m_train = getfielddef(cfg, {'black','train','n_m'}, 100);
            n_s_train = getfielddef(cfg, {'black','train','n_s'}, 10);
            n_k_train = getfielddef(cfg, {'black','train','n_k'}, n_k);
            params_true.testSuite_train = createTestSuite(sys, params_true, n_k_train, n_m_train, n_s_train, ts_opts);

            n_m_val = getfielddef(shared,'n_m_val', max(2,min(n_m,5)));
            n_s_val = getfielddef(shared,'n_s_val', n_s);
            n_k_val = getfielddef(shared,'n_k_val', n_k);
            params_true.testSuite_val = createTestSuite(sys, params_true, n_k_val, n_m_val, n_s_val, ts_opts);

            % ======== BLACK-BOX RCSI (conform) ========
            options = options_reach;
            options.cs = getfielddef(shared,'cs_base',struct());
            p = getfielddef(cfg, {'black','approx','p'}, getfielddef(sys,'n_p',1));
            options.approx = struct('p', p);

            apx = getfielddef(cfg,'black',struct()); apx = getfielddef(apx,'approx',struct());
            map_if(options, apx, 'gp_parallel');       map_if(options, apx, 'gp_pop_size');
            map_if(options, apx, 'gp_num_gen');        map_if(options, apx, 'gp_func_names');
            map_if(options, apx, 'gp_max_genes');      map_if(options, apx, 'gp_max_depth');
            map_if(options, apx, 'cgp_num_gen');       map_if(options, apx, 'cgp_pop_size_base');
            map_if(options, apx, 'save_res');

            % Initial guesses (centered sets)
            cR0 = center(R0);
            cU  = center(U);  cU = cU(:);
            params_id_init        = params_true;          % keep testSuite_* fields
            params_id_init.R0     = zonotope([cR0]);      % 0-radius
            params_id_init.U      = zonotope([cU, eye(numel(cU)) ones(numel(cU),1)]);
            params_id_init.n_p    = p;

            methods = getfielddef(cfg, {'black','methodsBlack'}, ["blackCGP"]);
            t_black_learn = NaN; learned_idx = 2;
            configs = cell(numel(methods)+1,1);

            % configs{1}: true
            configs{1}.sys     = sys;
            configs{1}.params  = rmfield(params_true,'testSuite');
            configs{1}.options = options_reach;
            configs{1}.name    = "true";

            for iM = 1:numel(methods)
                method = string(methods(iM));
                fprintf("Identification (BLACK %s): n_m=%d, n_s=%d, n_k=%d\n", method, n_m, n_s, n_k);
                t_id = tic;
                [configs{iM+1}.params, results] = conform(sys, params_id_init, options, method);
                if iM == 1, t_black_learn = toc(t_id); else, ~toc(t_id); end
                configs{iM+1}.sys     = results.sys;
                configs{iM+1}.options = options_reach;
                configs{iM+1}.name    = method;
            end

            % Validation (BLACK) — containment on validation suite
            t_val = tic;
            num_in = 0; num_out = 0;
            for mval = 1:length(params_true.testSuite_val)
                [~, eval_val] = validateReach(params_true.testSuite_val{mval}, configs, 1);
                num_in  = num_in  + eval_val.num_in(learned_idx);
                num_out = num_out + eval_val.num_out(learned_idx);
            end
            denom = max(1, num_in + num_out);
            cval_black   = 100 * (num_in / denom);
            t_black_val  = toc(t_val);

            % Conservatism proxy (BLACK) — size of intervals for learned model
            t_inf = tic;
            params_val = params_true;
            params_val.tFinal = sys.dt * (n_k_val - 1);
            Rlearn = reach(configs{learned_idx}.sys, params_val, options_reach);
            sizeI_black  = agg_interval_size(Rlearn);
            t_black_infer = toc(t_inf);

            % ======== DDRA-Lipschitz (Alg. 6) ========
            % Build data D from training suite
            D = build_D_from_suite(params_true.testSuite_train);
            nx = size(D.Xminus,1);

            % Disturbance W (possibly off for Black)
            if getfielddef(cfg.shared,'noise_for_ddra', true)
                eta_w  = getfielddef(cfg.ddra,'eta_w',1);
                alphaW = getfielddef(cfg.ddra,'alpha_w',0);
                W = zonotope(zeros(nx,1), alphaW*ones(nx, max(1,eta_w)));
            else
                W = zonotope(zeros(nx,1));
            end

            % Estimate Lipschitz & cover radius
            t_ddra_build_tic = tic;
            [Lvec, delta] = estimate_L_and_delta_from_data(D);
            Cnl = struct();
            Cnl.shared = struct('n_k', n_k, 'n_k_val', n_k_val);
            Cnl.lowmem = struct('zonotopeOrder_cap', getfielddef(cfg.lowmem,'zonotopeOrder_cap',100));
            Cnl.nlip   = struct('ZepsFlag', true, 'Lvec', Lvec(:), 'gamma', delta*ones(nx,1));
            t_ddra_build = toc(t_ddra_build_tic);

            % Validation with DDRA-Lip
            t_ddra_val_tic = tic;
            num_in_d = 0; num_out_d = 0; sum_sizeI_ddra = 0; t_reach_sum = 0;
            for mval = 1:length(params_true.testSuite_val)
                TS = params_true.testSuite_val{mval};
                Uk = repmat({params_true.U}, 1, n_k_val);

                t_case = tic;
                [Xsets, sizeI_case] = ddra_reach_lipschitz(params_true.R0, Uk, W, D, Cnl);
                t_reach_sum = t_reach_sum + toc(t_case);

                [nin, nout] = contain_points_in_sets(TS.y, Xsets);
                num_in_d  = num_in_d  + nin;
                num_out_d = num_out_d + nout;
                sum_sizeI_ddra = sum_sizeI_ddra + sizeI_case;
            end
            t_ddra_val       = toc(t_ddra_val_tic);
            t_ddra_reach_avg = t_reach_sum / max(1,length(params_true.testSuite_val));
            cval_ddra        = 100 * (num_in_d / max(1, num_in_d + num_out_d));
            sizeI_ddra       = sum_sizeI_ddra / max(1,length(params_true.testSuite_val));

            % -------- Row & write --------
            rowi = rowi + 1;
            row = {n_m, n_s, n_k, ...
                   cval_black, sizeI_black, ...
                   t_black_learn, t_black_val, t_black_infer, ...
                   cval_ddra, sizeI_ddra, ...
                   t_ddra_build, t_ddra_val, t_ddra_reach_avg};

            if append_csv
                fid = fopen(csv_path,'a');
                fprintf(fid, '%s\n', strjoin(cellfun(@num2str_cell, row,'uni',0), ',')); fclose(fid);
            else
                cells(rowi,:) = row;
            end

          end
        end
      end
    end

    if ~append_csv
        writecell([hdr; cells(1:rowi,:)], csv_path);
    end
    SUMMARY = readtable(csv_path);
    fprintf('Sweeps done. Rows: %d\nCSV -> %s\n', height(SUMMARY), csv_path);
end

% ---------- Totals for runtime panels (helper) ----------
function T = ensure_time_totals_square(T)
    addcol = @(tbl,name,val) iff(~ismember(name, tbl.Properties.VariableNames), ...
        @() addvars(tbl, val, 'NewVariableNames', name), ...
        @() tbl);
    iff = @(cond,a,b) (cond)*a() + (~cond)*b(); %#ok<NASGU> (trick to keep inline)

    if ~ismember('t_black_total', T.Properties.VariableNames)
        T.t_black_total = zeros(height(T),1);
    end
    if ~ismember('t_ddra_total', T.Properties.VariableNames)
        T.t_ddra_total = zeros(height(T),1);
    end
    if all(ismember({'t_black_learn','t_black_val','t_black_infer'}, T.Properties.VariableNames))
        T.t_black_total = T.t_black_learn + T.t_black_val + T.t_black_infer;
    end
    if all(ismember({'t_ddra_build','t_ddra_val','t_ddra_reach_avg'}, T.Properties.VariableNames))
        T.t_ddra_total = T.t_ddra_build + T.t_ddra_val + T.t_ddra_reach_avg;
    end
end

% =================== Shared helpers (light) ===================

function L = as_list(v), if isscalar(v), L = v; else, L = v(:)'; end, end
function v = getfielddef(S, path, def)
    if ischar(path) || isstring(path), path = {char(path)}; end
    v = S;
    for i=1:numel(path)
        if isstruct(v) && isfield(v, path{i}), v = v.(path{i}); else, v = def; return; end
    end
end
function s = num2str_cell(x)
    if ischar(x) || isstring(x), s = char(x);
    elseif islogical(x), s = char(string(x));
    else, s = num2str(x, '%.10g'); end
end
function S = agg_interval_size(R)
    S = 0;
    if ~isfield(R,'timePoint') || ~isfield(R.timePoint,'set'), return; end
    for k = 1:numel(R.timePoint.set)
        Ik = interval(R.timePoint.set{k});
        S = S + sum(abs(Ik.sup - Ik.inf), 'all');
    end
end
function x = coerce_numeric(col)
    if istable(col), col = col{:,:}; end
    if iscell(col), x = cellfun(@double, col); else, x = double(col); end
end
function map_if(options, apx, fname)
    if isfield(apx, fname)
        options.approx.(fname) = apx.(fname);
    end
end

% ------- RCSI label (1st method if any) -------
function lbl = rcsi_label_from_cfg(cfg)
    lbl = "black";
    try
        m = getfield(cfg,'black'); m = m.methodsBlack;
        if ~isempty(m), lbl = char(m(1)); end
    catch, end
end

% ===== DDRA glue (build data, estimate L/δ, containment) =====
function D = build_D_from_suite(S)
% Treat outputs y as state (Square p=1 -> x=y).
    Xm = []; Um = []; Xp = [];
    for m = 1:numel(S)
        Y = S{m}.y; U = S{m}.u;
        [~, nk, ns] = size(Y);
        for i = 1:ns
            Yi = Y(:,:,i); Ui = U(:,:,i);
            Xm = [Xm, Yi(:,1:end-1)];
            Um = [Um, Ui(:,1:end-1)];
            Xp = [Xp, Yi(:,2:end)];
        end
    end
    D = struct('Xminus', Xm, 'Uminus', Um, 'Xplus', Xp);
end

function [Lvec, delta] = estimate_L_and_delta_from_data(D)
% Per-dimension Lipschitz & covering radius from data (pragmatic).
    Z  = [D.Xminus; D.Uminus];    % z_i columns
    Xp = D.Xplus;
    n  = size(Xp,1); N = size(Z,2);
    Lvec = zeros(n,1); delta = 0;
    idx = 1:N; if N > 2000, idx = randperm(N,2000); end
    Zi  = Z(:,idx);  Xpi = Xp(:,idx);  Ni = numel(idx);
    for i = 1:Ni
        dists = vecnorm(Zi - Zi(:,i), 2, 1); dists(i) = inf;
        delta = max(delta, min(dists));
        [~,ord] = sort(dists,'descend');
        for pick = ord(1:min(10,Ni-1))
            dz = norm(Zi(:,pick) - Zi(:,i));
            if dz <= 0, continue; end
            diff = abs(Xpi(:,pick) - Xpi(:,i)) / dz;
            Lvec = max(Lvec, diff);
        end
    end
end

function [nin, nout] = contain_points_in_sets(Y, Xsets)
% Count how many observed points lie inside the DDRA set at each step.
% Y: [ny x nk x ns]; Xsets: {1..nk} with Xsets{1}=R0, compare k=2..nk
    [~, nk, ns] = size(Y);
    nin = 0; nout = 0;
    for i = 1:ns
        for k = 2:nk
            Ik = interval(Xsets{k});  % interval enclosure of X_k
            y  = Y(:,k,i);
            inside = all(y >= Ik.inf & y <= Ik.sup);
            if inside, nin = nin + 1; else, nout = nout + 1; end
        end
    end
end

function [Xsets, sizeI] = ddra_reach_lipschitz(R0, U, W, D, C)
% Alg. 6: local LS model + Z_L + optional Z_eps (from Lipschitz).
% R0: initial set (zonotope), U: single input set or cell{1..N}, W: process noise
% D: struct(Xminus,Uminus,Xplus); C.nlip: .ZepsFlag, .Lvec, .gamma
    if ~iscell(U), U = repmat({U}, 1, getfielddef(C.shared,'n_k_val', getfielddef(C.shared,'n_k',1))); end
    N  = numel(U);
    nx = size(center(R0),1);

    % LS linearization about sample means
    xS = mean(D.Xminus,2);
    uS = mean(D.Uminus,2);
    Phi = [ones(1,size(D.Xminus,2)); D.Xminus - xS; D.Uminus - uS];   % [1; x-x*; u-u*]
    Mhat = (D.Xplus) * pinv(Phi);

    % residual enclosure Z_L (component-wise box)
    res = D.Xplus - Mhat*Phi;
    rmax = max(res,[],2); rmin = min(res,[],2);
    ZL = zonotope(0.5*(rmax+rmin), diag(0.5*(rmax-rmin)));

    % Z_eps (optional): eps_i = 0.5 * L_i * gamma_i
    Zeps = zonotope(zeros(nx,1));
    if getfielddef(C.nlip,'ZepsFlag', true)
        Lvec  = getfielddef(C.nlip,'Lvec', []);
        gamma = getfielddef(C.nlip,'gamma',[]);
        if isempty(Lvec) || isempty(gamma)
            Lvec = zeros(nx,1); gamma = zeros(nx,1); % safe no-op
        end
        eps = 0.5 * abs(Lvec(:)) .* abs(gamma(:));
        Zeps = zonotope(zeros(nx,1), diag(eps));
    end

    Xsets  = cell(N+1,1); Xsets{1} = R0;
    sizeI  = 0;
    for k=1:N
        Uk = U{k};
        Xaff  = cartProd(Xsets{k} - xS, Uk - uS);
        Xnext = Mhat * aug1(Xaff) + W + ZL + Zeps;
        Xnext = reduce(Xnext, 'girard', getfielddef(C.lowmem,'zonotopeOrder_cap',100));
        Xsets{k+1} = Xnext;

        Iv = interval(Xnext); sizeI = sizeI + sum(abs(Iv.sup - Iv.inf));
    end
end

function A = aug1(S)
    c = center(S); G = generators(S);  % S in R^{n+m}
    A = zonotope([1; c], blkdiag(zeros(1,size(G,2)), G));
end
