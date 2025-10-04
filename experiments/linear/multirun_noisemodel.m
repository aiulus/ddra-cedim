%% noise_add_sweep_batch.m
% Headless multi-run driver for additive *measurement* noise sweep (alpha_W).
% - Single config per run; alpha_W is the TRAIN noise radius when meas.enable=true
% - Inference noise has its own toggle/magnitude
% - Gray's extra disturbance is zeroed automatically when meas.enable=true
% - Aggregates CSV like other runners

clear; clc;

%% ---------- Per-run knobs (set these, then run) ----------
NREP = 5;                                % repetitions (keep 1 for speed)
base_tag = 'singlerun_kMSD_noise_add';   % save_tag prefix
do_plots = true;

% Measurement-noise controls (per *run*)
MEAS_ENABLE      = true;     % ON: use measurement-noise experiment mode
MEAS_TRAIN_ON    = true;     % ON: inject bounded noise into regressors during TRAIN
MEAS_INFER_ON    = false;    % ON: add bounded regressor noise during INFERENCE
MEAS_ALPHA_INFER = 0.00;     % inference noise radius (state units)

% DDRA variant: "meas" (measurement-noise-aware) or "std"
DDRA_VARIANT = "meas";

% Sweep list (alpha_W = TRAIN noise radius when MEAS_ENABLE=true)
alpha_list = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.75, 1];

rng(1,'twister');

%% ---------- Base cfg ----------
cfg = struct();
cfg.io = struct('save_tag', base_tag, ...
                'plot_mode', "offline", ...
                'make_reach_plot', false, ...
                'save_artifacts', false);

cfg.shared = struct();
cfg.shared.dyn   = "k-Mass-SD";
cfg.shared.type  = "standard";
cfg.shared.p_extr = 0.3;
cfg.shared.options_reach = struct('zonotopeOrder',100,'verbose',false);
cfg.shared.cs_base = struct('robustnessMargin',1e-9,'verbose',false, ...
                            'cost',"interval",'constraints',"half");

% Budgets
cfg.shared.n_m = 10;
cfg.shared.n_s = 2;  
cfg.shared.n_k = 10;
cfg.shared.n_m_val = 5;
cfg.shared.n_s_val = 2; 
cfg.shared.n_k_val = 10;

% Unified noise policy (unchanged)
cfg.shared.noise_for_gray = true;
cfg.shared.noise_for_ddra = true;
cfg.shared.use_noise      = true;

% DDRA + Gray defaults
cfg.ddra = struct('eta_w',1,'alpha_w',0.01, ...
                  'allow_ridge',true,'lambda',1e-8, ...
                  'ridge_gamma',1.0,'ridge_policy',"MAB", ...
                  'variant', string(DDRA_VARIANT));
cfg.gray = struct('methodsGray', ["graySeq"], 'alpha_w', cfg.ddra.alpha_w, 'eta_w', cfg.ddra.eta_w);

cfg.lowmem = struct('gray_check_contain', true, 'store_ddra_sets', false, ...
                    'append_csv', true, 'zonotopeOrder_cap', 25);

cfg.metrics = struct();
cfg.metrics.enhanced = true;
cfg.metrics.directional = struct('enable', true, 'Nd', 64, 'seed', 12345);
cfg.metrics.safety = struct('enable', false, 'H', [], 'h', [], ...
                            'k_set', [], 'mode', 'per_block', ...
                            'tau_sweep', linspace(0,0.2,21), ...
                            'reach_reduce', 80);

cfg.io.base_dir    = fileparts(fileparts(mfilename('fullpath')));
cfg.allow_parallel = false;

% Per-run measurement-noise settings (consumed by run_sweeps)
cfg.meas = struct('enable',          MEAS_ENABLE, ...
                  'inject_on_train', MEAS_TRAIN_ON, ...
                  'infer_on',        MEAS_INFER_ON, ...
                  'alpha_infer',     MEAS_ALPHA_INFER);

% Sweep grid: alpha_W is interpreted in run_sweeps as TRAIN noise radius when meas.enable=true
sweep_grid = struct();
sweep_grid.D_list        = 2;
sweep_grid.alpha_w_list  = alpha_list;
sweep_grid.n_m_list      = cfg.shared.n_m;
sweep_grid.n_s_list      = cfg.shared.n_s;
sweep_grid.n_k_list      = cfg.shared.n_k;
sweep_grid.pe_list       = {struct('mode','randn','deterministic',true)};

rcsi_lbl = rcsi_label_from_cfg(cfg);

%% ---------- Run (like other batch scripts) ----------
run_tables = cell(NREP,1); res_dirs = strings(NREP,1);

for r = 1:NREP
    cfg_r = cfg;
    cfg_r.io.save_tag = sprintf('%s_%s_rep%02d', base_tag, rcsi_lbl, r);
    rng(r,'twister');

    cfg_r.io.plot_mode       = "offline";
    cfg_r.io.make_reach_plot = false;
    cfg_r.io.save_artifacts  = false;

    fprintf('\n=== Measurement-noise sweep — repetition %d/%d — tag: %s ===\n', r, NREP, cfg_r.io.save_tag);
    S = run_sweeps(cfg_r, sweep_grid);
    S = ensure_time_totals(S);

    run_tables{r} = S;

    [~, rd] = init_io(cfg_r);
    res_dirs(r) = string(rd);
    try
        src = fullfile(rd,'summary.csv');
        if isfile(src)
            copyfile(src, fullfile(rd, sprintf('summary_%s.csv', cfg_r.io.save_tag)));
        end
    catch ME
        warning('Could not copy per-run CSV (rep %d): %s', r, ME.message);
    end
end

%% ---------- Aggregate (unchanged pattern) ----------
Tall = table();
for r = 1:NREP
    if ~isempty(run_tables{r})
        T = run_tables{r}; T.rep = r*ones(height(T),1);
        Tall = [Tall; T]; %#ok<AGROW>
    end
end
if isempty(Tall), error('All repetitions failed—nothing to aggregate.'); end

% If you later run mixed variants in one go, add 'ddra_variant' (and even meas flags)
keys = {};
for key = ["alpha_w","D","n_m","n_s","n_k"]  % add "ddra_variant" here if mixing
    if ismember(key, string(Tall.Properties.VariableNames))
        keys{end+1} = char(key);
    end
end

isNum   = varfun(@isnumeric, Tall, 'OutputFormat','uniform');
numVars = string(Tall.Properties.VariableNames(isNum));
numVars = setdiff(numVars, [string(keys), "rep"]);

Gm = groupsummary(Tall, keys, "mean", cellstr(numVars));
Gs = groupsummary(Tall, keys, "std",  cellstr(numVars));
AGG = outerjoin(Gm, Gs, 'Keys', keys, 'MergeKeys', true, 'Type','left');

for v = numVars
    mn = "mean_" + v;  sd = "std_" + v;
    if ismember(mn, string(AGG.Properties.VariableNames)), AGG.(v + "_mean") = AGG.(mn);  AGG.(mn) = []; end
    if ismember(sd, string(AGG.Properties.VariableNames)), AGG.(v + "_std")  = AGG.(sd);  AGG.(sd) = []; end
end

[ddra_lbl, gray_lbl, hdr] = legend_bits_for_batch(cfg, sweep_grid, 'alpha_w', rcsi_lbl);
colors = struct('ddra',[0.23 0.49 0.77],'gray',[0.85 0.33 0.10]);

agg_dir = fullfile(cfg.io.base_dir,'experiments','results','data', ...
                   sprintf('%s_%s_agg', base_tag, rcsi_lbl));
if ~exist(agg_dir,'dir'), mkdir(agg_dir); end
agg_csv = fullfile(agg_dir,'summary_agg.csv');
writetable(AGG, agg_csv);
fprintf('\nWrote aggregated CSV: %s\n', agg_csv);

if do_plots
  try
    x = AGG.alpha_w;  M = @(name) AGG.(name + "_mean");
    figure('Color','w','Name','Fidelity & Conservatism (mean)'); 
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    nexttile; hold on; grid on;
    title(sprintf('%s — Fidelity vs \\alpha_W', hdr));
    plot(x, M("cval_ddra"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("cval_gray"), '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('\alpha_W'); ylabel('Containment on validation (%)'); legend('Location','best');

    nexttile; hold on; grid on;
    title(sprintf('%s — Conservatism vs \\alpha_W', hdr));
    plot(x, M("sizeI_ddra"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("sizeI_gray"), '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('\alpha_W'); ylabel('\Sigma interval widths'); legend('Location','best');
  catch ME
    fprintf('Aggregate plotting failed: %s\n(Your CSVs are already saved at %s)\n', ME.message, agg_dir);
  end
end
