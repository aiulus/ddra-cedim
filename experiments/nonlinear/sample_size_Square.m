%% HOW TO USE (Square / Nonlinear Sample-Size Sweep)
% What it does:
%   Uses CORA's simple NARX "Square" system and sweeps n_m (number of input
%   trajectories). Produces fidelity/conservatism plots for black-box RCSI
%   and Lipschitz-DDRA (Alg. 6), and plots runtime profiles.
%
% Outputs:
%   CSV + plots via init_io() under experiments/results/{data,plots}/<tag>_sweeps

rng(1,'twister');

% ---------- Minimal cfg ----------
cfg = struct();
cfg.io = struct('save_tag', 'Square_sample_size_sweep');

% System & shared options
cfg.shared = struct();
cfg.shared.dyn    = "Square";      % <-- NARX system (p=1)
cfg.shared.type   = "standard";    % uncertainty preset
cfg.shared.p_extr = 0.3;           % test-suite extreme input prob.

% CORA reachability options
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

% Data budgets
cfg.shared.n_m     = 10;     % identification: # input trajectories
cfg.shared.n_s     = 10;    % samples per input trajectory
cfg.shared.n_k     = 4;     % horizon (train)
cfg.shared.n_m_val = cfg.shared.n_m;   % validation: # input trajectories
cfg.shared.n_s_val = cfg.shared.n_s;
cfg.shared.n_k_val = cfg.shared.n_k;

% Black-box RCSI methods
cfg.black = struct();
cfg.black.methodsBlack = ["blackCGP"];   % or ["blackGP","blackCGP"]

% Optional: tiny/fast GPTIPS controls (read by CORA's config shim)
cfg.black.approx = struct( ...
    'p', 1, ...                               % ARX order (Square has p=1)
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

% (Optional) use "rand" to match CORA example’s randomness
rng(2,'twister');
cfg.shared.type = "rand";

% Keep explicit process noise out for side-by-side comparison
cfg.shared.noise_for_gray = false;

% --- method label for filenames
rcsi_lbl = "blackGP";
if numel(cfg.black.methodsBlack) > 0
    rcsi_lbl = char(cfg.black.methodsBlack(1));
end
cfg.io.save_tag = sprintf('%s_%s', cfg.io.save_tag, rcsi_lbl);

% ---------- Sweep grid ----------
sweep_grid = struct();
% sweep_grid.n_m_list = [2 4 8 16];   % uncomment to sweep
% defaults to cfg.shared.n_m if not set

% ---------- Low-memory toggles ----------
cfg.lowmem = struct();
cfg.lowmem.append_csv = true;               % stream CSV rows
cfg.lowmem.zonotopeOrder_cap = 100;         % used by DDRA reduction

% ---------- Run ----------
SUMMARY = run_sweeps_square_rcsi(cfg, sweep_grid);

% ---------- Plots (Fidelity + Conservatism vs n_m) ----------
x_nm        = coerce_numeric(SUMMARY.n_m);
cval_bb     = coerce_numeric(SUMMARY.cval_black);
sizeI_bb    = coerce_numeric(SUMMARY.sizeI_black);
cval_ddra   = coerce_numeric(SUMMARY.cval_ddra);
sizeI_ddra  = coerce_numeric(SUMMARY.sizeI_ddra);

colors = struct('bb',[0.2 0.55 0.3], 'ddra',[0.25 0.35 0.9]);

f = figure('Name','Square | Fidelity & Conservatism vs n_m','Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% (1) Fidelity
nexttile; hold on; grid on;
plot(x_nm, cval_bb,   '-o', 'Color',colors.bb,   'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl]);
plot(x_nm, cval_ddra, '-s', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName','DDRA-Lip');
xlabel('n_m (input trajectories)'); ylabel('Containment on validation (%)');
title(['Fidelity vs n_m  (Square)']); legend('Location','best');

% (2) Conservatism proxy
nexttile; hold on; grid on;
plot(x_nm, sizeI_bb,   '-o', 'Color',colors.bb,   'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl]);
plot(x_nm, sizeI_ddra, '-s', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName','DDRA-Lip');
xlabel('n_m (input trajectories)'); ylabel('Aggregated interval size (proxy)');
title(['Conservatism vs n_m  (Square)']); legend('Location','best');

% Save
[plots_dir, ~] = init_io(cfg);
save_plot(f, plots_dir, ['square_fid_cons_vs_nm_RCSI_vs_DDRA'], 'Formats', {'png','pdf'}, 'Resolution', 200);
close all force

% ---------- Plots (Runtime Profiles vs n_m) ----------
rt_nm            = x_nm;
t_black_learn    = coerce_numeric(SUMMARY.t_black_learn);
t_black_val      = coerce_numeric(SUMMARY.t_black_val);
t_black_infer    = coerce_numeric(SUMMARY.t_black_infer);
t_ddra_build     = coerce_numeric(SUMMARY.t_ddra_build);
t_ddra_val       = coerce_numeric(SUMMARY.t_ddra_val);
t_ddra_reach_avg = coerce_numeric(SUMMARY.t_ddra_reach_avg);

f2 = figure('Name','Square | Runtime Profiles vs n_m','Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% (1) Black-box RCSI times
nexttile; hold on; grid on;
plot(rt_nm, t_black_learn, '-o', 'LineWidth',1.6, 'Color',colors.bb, 'DisplayName','RCSI learn');
plot(rt_nm, t_black_val,   '-s', 'LineWidth',1.6, 'Color',[colors.bb*0.7 0.8], 'DisplayName','RCSI validate');
plot(rt_nm, t_black_infer, '-^', 'LineWidth',1.6, 'Color',[colors.bb*0.5 0.6], 'DisplayName','RCSI infer');
xlabel('n_m (input trajectories)'); ylabel('time (s)');
title('Runtime vs n_m (Black-box)'); legend('Location','best');

% (2) DDRA times
nexttile; hold on; grid on;
plot(rt_nm, t_ddra_build,     '-o', 'LineWidth',1.6, 'Color',colors.ddra, 'DisplayName','DDRA build');
plot(rt_nm, t_ddra_val,       '-s', 'LineWidth',1.6, 'Color',[colors.ddra*0.7 0.8], 'DisplayName','DDRA validate (total)');
plot(rt_nm, t_ddra_reach_avg, '-^', 'LineWidth',1.6, 'Color',[colors.ddra*0.5 0.6], 'DisplayName','DDRA reach/case (avg)');
xlabel('n_m (input trajectories)'); ylabel('time (s)');
title('Runtime vs n_m (DDRA)'); legend('Location','best');

save_plot(f2, plots_dir, ['square_runtime_profiles_vs_nm_RCSI_vs_DDRA'], 'Formats', {'png','pdf'}, 'Resolution', 200);
close all force


%% ====================== Local functions ======================

function SUMMARY = run_sweeps_square_rcsi(cfg, grid)
% Sweep over n_m for CORA's "Square" NARX system using black-box RCSI
% and Lipschitz-DDRA (Alg. 6). Requires ddra_reach_lipschitz() in path.

    [~, results_dir] = init_io(cfg);
    csv_path = fullfile(results_dir, 'summary.csv');

    LM = getfielddef(cfg,'lowmem',struct());
    append_csv = getfielddef(LM,'append_csv',true);

    % Axes
    axes = struct();
    axes.n_m = as_list(getfielddef(grid,'n_m_list', getfielddef(cfg,'shared',struct()).n_m));
    axes.n_s = as_list(getfielddef(grid,'n_s_list', getfielddef(cfg,'shared',struct()).n_s));
    axes.n_k = as_list(getfielddef(grid,'n_k_list', getfielddef(cfg,'shared',struct()).n_k));

    % Header (now includes DDRA columns + reach avg)
    hdr = {'n_m','n_s','n_k', ...
           'cval_black','sizeI_black', ...
           't_black_learn','t_black_val','t_black_infer', ...
           'cval_ddra','sizeI_ddra', ...
           't_ddra_build','t_ddra_val','t_ddra_reach_avg'};

    if append_csv
        fid = fopen(csv_path,'w'); fprintf(fid,'%s\n',strjoin(hdr,',')); fclose(fid);
    else
        Ntot = numel(axes.n_m)*numel(axes.n_s)*numel(axes.n_k);
        cells = cell(Ntot, numel(hdr));
    end

    rowi = 0;

    for n_m = axes.n_m
      for n_s = axes.n_s
        for n_k = axes.n_k

            % ---------- Build system ----------
            shared = getfielddef(cfg,'shared',struct());
            dyn    = getfielddef(shared,'dyn',"Square");
            type   = getfielddef(shared,'type',"standard");
            p_extr = getfielddef(shared,'p_extr',0.3);

            [sys, R0, U, ~] = custom_loadDynamics(dyn, type);
            options_reach = getfielddef(shared,'options_reach',struct());

            % Suites
            params_true = struct();
            params_true.R0     = R0;
            params_true.U      = U;
            params_true.tFinal = sys.dt * (n_k - 1);

            ts_opts = struct('p_extr', p_extr);

            % (a) Identification suite (sanity)
            params_true.testSuite = createTestSuite(sys, params_true, n_k, n_m, n_s, ts_opts);

            % (b) TRAIN suite (required)
            n_m_train = getfielddef(cfg, {'black','train','n_m'}, 100);
            n_s_train = getfielddef(cfg, {'black','train','n_s'}, 10);
            n_k_train = getfielddef(cfg, {'black','train','n_k'}, n_k);
            params_true.testSuite_train = createTestSuite(sys, params_true, n_k_train, n_m_train, n_s_train, ts_opts);

            % (c) VALIDATION suite
            n_m_val = getfielddef(shared,'n_m_val',5);
            n_s_val = getfielddef(shared,'n_s_val',n_s);
            n_k_val = getfielddef(shared,'n_k_val',n_k);
            params_true.testSuite_val = createTestSuite(sys, params_true, n_k_val, n_m_val, n_s_val, ts_opts);

            % ---------- Conformance options (Black-Box) ----------
            options = options_reach;
            options.cs = getfielddef(shared,'cs_base',struct());

            % ARX order p (Square has p=1)
            p = getfielddef(cfg, {'black','approx','p'}, 1);

            % Base approx options (CGP + GP knobs)
            options.approx = struct( ...
                'cgp_num_gen',        getfielddef(cfg, {'black','approx','cgp_num_gen'}, 5), ...
                'cgp_pop_size_base',  getfielddef(cfg, {'black','approx','cgp_pop_size_base'}, 10), ...
                'gp_parallel',        getfielddef(cfg, {'black','approx','gp_parallel'}, false), ...
                'save_res',           getfielddef(cfg, {'black','approx','save_res'}, false), ...
                'p',                  p ...
            );

            % Optional GP fine controls (read by the config_gp shim)
            apx = getfielddef(cfg,'black',struct()); apx = getfielddef(apx,'approx',struct());
            map_if(options, apx, 'gp_pop_size');      map_if(options, apx, 'gp_num_gen');
            map_if(options, apx, 'gp_max_genes');     map_if(options, apx, 'gp_max_depth');
            map_if(options, apx, 'gp_func_names');    map_if(options, apx, 'gp_fitness_goal');
            map_if(options, apx, 'gp_stall_limit');   map_if(options, apx, 'gp_tournament_size');
            map_if(options, apx, 'gp_elite_fraction');map_if(options, apx, 'gp_lexicographic');
            map_if(options, apx, 'gp_const_range');   map_if(options, apx, 'gp_seed');

            % ---------- Initial guesses (Black-Box) ----------
            cR0 = center(R0);
            cU  = center(U);  cU = cU(:);
            params_id_init        = params_true;        % keep testSuite_* fields
            params_id_init.R0     = zonotope([cR0]);    % 0-radius
            params_id_init.U      = zonotope([cU, eye(numel(cU)) ones(numel(cU),1)]);
            params_id_init.tFinal = sys.dt * (n_k_train - 1);
            params_id_init.n_p    = p;

            % ---------- Run black-box identification ----------
            methods = getfielddef(cfg, {'black','methodsBlack'}, ["blackCGP"]);
            configs = cell(numel(methods)+1,1);
            configs{1}.sys     = sys;
            configs{1}.params  = rmfield(params_true,'testSuite');
            configs{1}.options = options_reach;
            configs{1}.name    = "true";

            t_black_learn = NaN; learned_idx = 2;
            for iM = 1:numel(methods)
                method = string(methods(iM));
                fprintf("Identification with method %s  (n_m=%d, n_s=%d, n_k=%d)\n", method, n_m, n_s, n_k);
                t_id = tic;
                [configs{iM+1}.params, results] = conform(sys, params_id_init, options, method);
                if iM == 1, t_black_learn = toc(t_id); else, ~toc(t_id); end
                configs{iM+1}.sys     = results.sys;
                configs{iM+1}.options = options_reach;
                configs{iM+1}.name    = method;
            end

            % ---------- Validation (Black-Box) ----------
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

            % ---------- Conservatism proxy (Black-Box) ----------
            t_inf = tic;
            params_val = params_true;
            params_val.tFinal = sys.dt * (n_k_val - 1);
            Rlearn = reach(configs{learned_idx}.sys, params_val, options_reach);
            sizeI_black  = agg_interval_size(Rlearn);
            t_black_infer = toc(t_inf);

            % ================= DDRA-Lipschitz (Alg. 6) =================
            % Build D and estimate L, delta
            t_ddra_build_tic = tic;
            D = build_D_from_suite(params_true.testSuite_train);   % columns = samples
            nx = size(D.Xminus,1);
            [Lvec, delta] = estimate_L_and_delta_from_data(D);
            C = struct();
            C.shared = struct('n_k',      cfg.shared.n_k, ...
                              'n_k_val',  cfg.shared.n_k_val);
            C.lowmem = struct('zonotopeOrder_cap', getfielddef(cfg.lowmem,'zonotopeOrder_cap',100));
            C.nlip   = struct('ZepsFlag', true, ...
                              'Lvec', Lvec(:), ...
                              'gamma', delta * ones(nx,1));
            W = zonotope(zeros(nx,1));
            t_ddra_build = toc(t_ddra_build_tic);

            % Validation with DDRA (+ reach/case times)
            t_ddra_val_tic = tic;
            num_in_d = 0; num_out_d = 0; sum_sizeI_ddra = 0;
            t_ddra_reach_sum = 0;
            for mval = 1:length(params_true.testSuite_val)
                nsteps = cfg.shared.n_k_val;
                Uk = repmat({params_true.U}, 1, nsteps);
                R0_ddra = params_true.R0;

                t_case = tic;
                [Xsets, ~] = ddra_reach_lipschitz(R0_ddra, Uk, W, D, C);
                t_ddra_reach_sum = t_ddra_reach_sum + toc(t_case);

                % containment
                T = params_true.testSuite_val{mval};
                [nin, nout] = contain_points_in_sets(T.y, Xsets);
                num_in_d  = num_in_d  + nin;
                num_out_d = num_out_d + nout;

                % size proxy for this case
                sizeI_case = agg_interval_size(struct('timePoint',struct('set',{Xsets(2:end)})));
                sum_sizeI_ddra = sum_sizeI_ddra + sizeI_case;
            end
            t_ddra_val       = toc(t_ddra_val_tic);
            t_ddra_reach_avg = t_ddra_reach_sum / max(1,length(params_true.testSuite_val));
            denom_d   = max(1, num_in_d + num_out_d);
            cval_ddra = 100 * (num_in_d / denom_d);
            sizeI_ddra = sum_sizeI_ddra / max(1,length(params_true.testSuite_val));

            % ---------- Row & write ----------
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

    if ~append_csv
        writecell([hdr; cells(1:rowi,:)], csv_path);
    end
    SUMMARY = readtable(csv_path);
    fprintf('Sweeps done. Rows: %d\nCSV -> %s\n', height(SUMMARY), csv_path);
end

% --------- helpers ---------
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

% ===== DDRA glue (build data, estimate L/δ, containment) =====
function D = build_D_from_suite(S)
% S: cell array of test cases from createTestSuite(). Use training/val suite.
% We treat outputs y as the state (NARX p=1 -> x=y).
    Xm = []; Um = []; Xp = [];
    for m = 1:numel(S)
        Y = S{m}.y;   % [ny x nk x ns]
        U = S{m}.u;   % [nu x nk x ns]
        [~, nk, ns] = size(Y);
        for i = 1:ns
            Yi = Y(:,:,i); Ui = U(:,:,i);
            Xm = [Xm, Yi(:,1:end-1)];             % x(k)
            Um = [Um, Ui(:,1:end-1)];             % u(k)
            Xp = [Xp, Yi(:,2:end)];               % x(k+1)
        end
    end
    D = struct('Xminus', Xm, 'Uminus', Um, 'Xplus', Xp);
end

function [Lvec, delta] = estimate_L_and_delta_from_data(D)
% Per-dimension Lipschitz & covering radius from data.
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
% Y: [ny x nk x ns]; Xsets: {1..nk} with Xsets{1}=R0, so compare k=2..nk
    [~, nk, ns] = size(Y);
    nin = 0; nout = 0;
    for i = 1:ns
        for k = 2:nk
            Ik = interval(Xsets{k});  % CORA interval enclosure of X_k
            y  = Y(:,k,i);
            inside = all(y >= Ik.inf & y <= Ik.sup);
            if inside, nin = nin + 1; else, nout = nout + 1; end
        end
    end
end
