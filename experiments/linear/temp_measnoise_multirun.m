%% noise_train_meas_sweep_batch.m
% Headless multi-run driver to study *training measurement noise* levels.
% - Sweeps V_train_scale (measurement-noise generator scale) used in TRAIN only
% - Uses DDRA "meas" variant (Algorithm 4) and compares with Gray
% - Runs NREP repetitions per noise level
% - Streams per-run CSVs; aggregates to summary_agg.csv
% - Aggregate plots are try/catch guarded so CSVs always persist

clear; clc;

% ----- Controls -----
NREP        = 5;                            % repetitions per noise level
VTRAIN_LIST = [0, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2];  % TRAIN meas-noise scales
VAL_MEAS_MODE = "none";                     % 'none' (default), 'match', or numeric scale
base_tag    = 'kMSD_noise_train_meas';
do_plots    = true;

rng(1,'twister');                           % base seed; per-run uses rng(r + offset)

% ----- Base cfg (headless per-run) -----
cfg = struct();
cfg.io = struct('save_tag', base_tag, ...
                'plot_mode', "offline", ...
                'make_reach_plot', false, ...
                'save_artifacts', false);

% Shared
cfg.shared = struct();
cfg.shared.dyn   = "k-Mass-SD";
cfg.shared.type  = "standard";
cfg.shared.p_extr = 0.3;
cfg.shared.options_reach = struct('zonotopeOrder',100,'verbose',false);
cfg.shared.cs_base = struct('robustnessMargin',1e-9,'verbose',false, ...
                            'cost',"interval",'constraints',"half");

% Data budgets (modest)
cfg.shared.n_m     = 30;   cfg.shared.n_s     = 5;  cfg.shared.n_k     = 6;
cfg.shared.n_m_val = 2;    cfg.shared.n_s_val = cfg.shared.n_s;
cfg.shared.n_k_val = cfg.shared.n_k;

% Unified process-noise policy (state disturbance W)
cfg.shared.noise_for_gray = true;
cfg.shared.noise_for_ddra = true;
cfg.shared.use_noise      = true;

% DDRA settings
cfg.ddra = struct('eta_w',1,'alpha_w',0.01, ...
                  'allow_ridge',true,'lambda',1e-8, ...
                  'ridge_gamma',1.0,'ridge_policy',"MAB", ...
                  'variant',"meas");          % <- measurement-noise aware variant

% Gray
cfg.gray = struct('methodsGray', ["graySeq"]);

% Efficiency toggles (CSV streaming; no set storage)
cfg.lowmem = struct('gray_check_contain', true, ...
                    'store_ddra_sets',   false, ...
                    'append_csv',        true, ...
                    'zonotopeOrder_cap', 50);

% Metrics
cfg.metrics = struct();
cfg.metrics.enhanced   = true;
cfg.metrics.directional = struct('enable', true, 'Nd', 64, 'seed', 12345);
cfg.metrics.safety = struct('enable', false, 'H', [], 'h', [], ...
                            'k_set', [], 'mode', 'per_block', ...
                            'tau_sweep', linspace(0,0.2,21), ...
                            'reach_reduce', 80);

cfg.io.base_dir    = fileparts(fileparts(mfilename('fullpath')));
cfg.allow_parallel = false;

% Sweep grid (keep alpha_W etc. fixed here; we sweep V_train externally)
sweep_grid = struct();
sweep_grid.D_list        = 2;
sweep_grid.alpha_w_list  = cfg.ddra.alpha_w;
sweep_grid.n_m_list      = cfg.shared.n_m;
sweep_grid.n_s_list      = cfg.shared.n_s;
sweep_grid.n_k_list      = cfg.shared.n_k;
sweep_grid.pe_list       = {struct('mode','randn','deterministic',true)};

% Label used in folder names
rcsi_lbl = rcsi_label_from_cfg(cfg);

% ----- Preflight: if safety enabled, auto-fit H/h to ny -----
try
  if isfield(cfg,'metrics') && isfield(cfg.metrics,'safety') && cfg.metrics.safety.enable
      [~, baseC0] = init_sweep_axes(cfg, sweep_grid);
      [sys_cora, ~, ~, ~] = build_true_system(baseC0);
      ny = size(sys_cora.C,1);
      H = getfield(cfg.metrics.safety,'H',[]);
      h = getfield(cfg.metrics.safety,'h',[]);
      if isempty(H) || size(H,2)~=ny || isempty(h) || size(h,1)~=size(H,1)
          cfg.metrics.safety.H = [ eye(ny); -eye(ny) ];
          cfg.metrics.safety.h = 10*ones(2*ny,1);
      end
  end
catch ME
  fprintf('Safety preflight skipped: %s\n', ME.message);
end

% ----- Run over VTRAIN_LIST × NREP -----
run_tables = cell(numel(VTRAIN_LIST), NREP);
res_dirs   = strings(numel(VTRAIN_LIST), NREP);

for iv = 1:numel(VTRAIN_LIST)
    vscale = VTRAIN_LIST(iv);

    for r = 1:NREP
        cfg_r = cfg;
        cfg_r.io.save_tag = sprintf('%s_%s_vtrain%.4g_rep%02d', base_tag, rcsi_lbl, vscale, r);

        % TRAIN measurement noise: enable + scaled eye(nx) at data gen time
        % (The run_sweeps bridge allows either cfg.ddra.noise.meas.* or cfg.data.train/val.meas.*)
        cfg_r.ddra.noise = struct('meas', struct('enable', true, 'gen_scale', vscale));
        % Also mirror under cfg.data for immediate clarity (bridge will no-op if already set)
        cfg_r.data.train = struct('meas', struct('enable', true, 'gen_scale', vscale));

        % VAL measurement noise policy
        switch string(VAL_MEAS_MODE)
            case "none"
                cfg_r.data.val = struct('meas', struct('enable', false));
            case "match"
                cfg_r.data.val = struct('meas', struct('enable', true, 'gen_scale', vscale));
            otherwise
                % numeric scale
                if isnumeric(VAL_MEAS_MODE)
                    cfg_r.data.val = struct('meas', struct('enable', true, 'gen_scale', VAL_MEAS_MODE));
                else
                    cfg_r.data.val = struct('meas', struct('enable', false));
                end
        end

        % Enforce headless per-run
        cfg_r.io.plot_mode       = "offline";
        cfg_r.io.make_reach_plot = false;
        cfg_r.io.save_artifacts  = false;

        % Different RNG per (iv,r) without destroying determinism
        rng(uint32( 1 + 7919*iv + 104729*r ), 'twister');

        fprintf('\n=== TRAIN meas-noise sweep — v_train=%.4g — rep %d/%d — tag: %s ===\n', ...
                vscale, r, NREP, cfg_r.io.save_tag);

        S = run_sweeps(cfg_r, sweep_grid);
        S = ensure_time_totals(S);

        % Tag the current training noise level into the table (so aggregation has the key)
        S.v_train = vscale*ones(height(S),1);

        run_tables{iv, r} = S;

        % Copy per-run CSV for traceability; also write a version including v_train
        try
            [~, rd] = init_io(cfg_r);
            res_dirs(iv, r) = string(rd);
            src = fullfile(rd,'summary.csv');
            if isfile(src)
                copyfile(src, fullfile(rd, sprintf('summary_%s.csv', cfg_r.io.save_tag)));
                % Write an augmented CSV with v_train column for convenience
                writetable(S, fullfile(rd, sprintf('summary_%s_with_v.csv', cfg_r.io.save_tag)));
            end
        catch ME
            warning('Per-run CSV copy/write failed (v=%.4g, rep %d): %s', vscale, r, ME.message);
        end
    end
end

% ----- Aggregate across runs (mean & std by v_train, alpha_w, D, n_m, n_s, n_k) -----
Tall = table();
for iv = 1:numel(VTRAIN_LIST)
  for r = 1:NREP
    T = run_tables{iv,r};
    if ~isempty(T)
        T.rep = r*ones(height(T),1);
        Tall  = [Tall; T];          
    end
  end
end
if isempty(Tall), error('All repetitions failed—nothing to aggregate.'); end

% Grouping keys (present columns only)
keys = {};
for key = ["v_train","alpha_w","D","n_m","n_s","n_k"]
    if ismember(key, string(Tall.Properties.VariableNames))
        keys{end+1} = char(key); 
    end
end

% Numeric variables to aggregate (exclude keys + rep)
isNum   = varfun(@isnumeric, Tall, 'OutputFormat','uniform');
numVars = string(Tall.Properties.VariableNames(isNum));
numVars = setdiff(numVars, [string(keys), "rep"]);

Gm  = groupsummary(Tall, keys, "mean", cellstr(numVars));
Gs  = groupsummary(Tall, keys, "std",  cellstr(numVars));
AGG = outerjoin(Gm, Gs, 'Keys', keys, 'MergeKeys', true, 'Type','left');

% Flatten mean_/std_ prefixes into *_mean / *_std
for v = numVars
    mn = "mean_" + v;  sd = "std_" + v;
    if ismember(mn, string(AGG.Properties.VariableNames))
        AGG.(v + "_mean") = AGG.(mn);  AGG.(mn) = [];
    end
    if ismember(sd, string(AGG.Properties.VariableNames))
        AGG.(v + "_std")  = AGG.(sd);  AGG.(sd) = [];
    end
end

agg_dir = fullfile(cfg.io.base_dir, 'experiments','results','data', ...
                   sprintf('%s_%s_agg', base_tag, rcsi_lbl));
if ~exist(agg_dir,'dir'), mkdir(agg_dir); end

agg_csv = fullfile(agg_dir, 'summary_agg.csv');
writetable(AGG, agg_csv);
fprintf('\nWrote aggregated CSV: %s\n', agg_csv);

% ----- Optional aggregate plots (guarded; CSV already saved) -----
if do_plots
try
    % fix x ordering
    AGG = sortrows(AGG, "v_train");
    x   = AGG.v_train;

    M = @(name) AGG.(name + "_mean");
    % S = @(name) AGG.(name + "_std");  % kept for future error-bar variants

    % Fidelity / Conservatism vs v_train
    figure('Color','w','Name','Fidelity & Conservatism vs training meas-noise'); 
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile; hold on; grid on; title('Fidelity vs v_{train}');
    plot(x, M("cval_ddra"), '-o','LineWidth',1.6, 'DisplayName','DDRA (meas)');
    plot(x, M("cval_gray"), '-s','LineWidth',1.6, 'DisplayName','GraySeq');
    xlabel('v_{train} scale'); ylabel('Containment on validation (%)'); legend('Location','best');

    nexttile; hold on; grid on; title('Conservatism proxy vs v_{train}');
    plot(x, M("sizeI_ddra"), '-o','LineWidth',1.6, 'DisplayName','DDRA (meas)');
    plot(x, M("sizeI_gray"), '-s','LineWidth',1.6, 'DisplayName','GraySeq');
    xlabel('v_{train} scale'); ylabel('\Sigma interval widths'); legend('Location','best');

    try
        saveas(gcf, fullfile(agg_dir, 'fidcons_vs_vtrain_mean.png'));
        exportgraphics(gcf, fullfile(agg_dir, 'fidcons_vs_vtrain_mean.pdf'),'ContentType','vector');
    catch, end

    % Shape-aware (if present)
    hasDir = ismember("dir_eps_med_mean", string(AGG.Properties.VariableNames)) || ...
             ismember("hout_p90_mean",    string(AGG.Properties.VariableNames)) || ...
             ismember("haus_sym_med_mean",string(AGG.Properties.VariableNames));
    if hasDir
        figure('Color','w','Name','Shape-aware vs training noise');
        tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

        if ismember("dir_eps_med_mean", string(AGG.Properties.VariableNames))
            nexttile; hold on; grid on; title('dir\_eps median');
            plot(x, M("dir_eps_med"), '-o','LineWidth',1.6); yline(1,'--');
            xlabel('v_{train} scale'); ylabel('ratio');
        end
        if ismember("hout_p90_mean", string(AGG.Properties.VariableNames))
            nexttile; hold on; grid on; title('H_{out} p90');
            plot(x, M("hout_p90"), '-o','LineWidth',1.6);
            xlabel('v_{train} scale'); ylabel('gap');
        end
        if ismember("haus_sym_med_mean", string(AGG.Properties.VariableNames))
            nexttile; hold on; grid on; title('Hausdorff (sym) median');
            plot(x, M("haus_sym_med"), '-o','LineWidth',1.6);
            xlabel('v_{train} scale'); ylabel('distance');
        end

        try
            saveas(gcf, fullfile(agg_dir, 'shapeaware_vs_vtrain_mean.png'));
            exportgraphics(gcf, fullfile(agg_dir, 'shapeaware_vs_vtrain_mean.pdf'),'ContentType','vector');
        catch, end
    end

    fprintf('Aggregate plots saved under: %s\n', agg_dir);
catch ME
    fprintf('Aggregate plotting failed: %s\n(Your CSVs are already saved at %s)\n', ME.message, agg_dir);
end
end
