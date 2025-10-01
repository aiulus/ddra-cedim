%% noise_add_sweep_batch.m
% Headless multi-run driver for the additive process-noise sweep (alpha_W).
% - Runs NREP repetitions with unique save_tag and RNG
% - Streams per-run CSVs (no per-run plotting/artifacts)
% - Aggregates to summary_agg.csv
% - Optional aggregate plots are wrapped in try/catch so CSVs persist

clear; clc;

NREP = 5;                      % how many repetitions
base_tag = 'kMSD_noise_add';    % base save_tag prefix
do_plots = true;                % aggregate plots at the end (safe-guarded)

rng(1,'twister');               % base seed; per-run we set rng(r)

%% ----- Build base cfg + sweep grid (no plotting during each run) -----
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
cfg.shared.n_m = 10;  
cfg.shared.n_s = 2;  
cfg.shared.n_k = 10;
cfg.shared.n_m_val = 10;
cfg.shared.n_s_val = 2;
cfg.shared.n_k_val = cfg.shared.n_k;

% Unified noise policy
cfg.shared.noise_for_gray = true;
cfg.shared.noise_for_ddra = true;
cfg.shared.use_noise      = true;

% DDRA noise template (eta fixed; alpha_w swept)
cfg.ddra = struct('eta_w',1,'alpha_w',0.01, ...
                  'allow_ridge',true,'lambda',1e-8, ...
                  'ridge_gamma',1.0,'ridge_policy',"MAB");

% Gray method(s)
cfg.gray = struct('methodsGray', ["graySeq"]);

% Efficiency toggles (stream CSV; no set storage)
cfg.lowmem = struct('gray_check_contain', true, ...
                    'store_ddra_sets',   false, ...
                    'append_csv',        true, ...
                    'zonotopeOrder_cap', 25);

% Metrics
cfg.metrics = struct();
cfg.metrics.enhanced = true;                 % richer, cheap metrics
cfg.metrics.directional = struct('enable', true, 'Nd', 64, 'seed', 12345);

% Optional safety metric block (disabled by default for speed)
% To enable ROC-style conservatism eval, set enable=true; we auto-fill H/h dims.
cfg.metrics.safety = struct('enable', false, 'H', [], 'h', [], ...
                            'k_set', [], 'mode', 'per_block', ...
                            'tau_sweep', linspace(0,0.2,21), ...
                            'reach_reduce', 80);

cfg.io.base_dir    = fileparts(fileparts(mfilename('fullpath')));
cfg.allow_parallel = false;

% Sweep grid: vary additive process noise scale α_W
sweep_grid = struct();
sweep_grid.D_list        = 2;
sweep_grid.alpha_w_list  = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.75, 1];
sweep_grid.n_m_list      = cfg.shared.n_m;
sweep_grid.n_s_list      = cfg.shared.n_s;
sweep_grid.n_k_list      = cfg.shared.n_k;
sweep_grid.pe_list       = {struct('mode','randn','deterministic',true)};

% Label used in folder names
rcsi_lbl = rcsi_label_from_cfg(cfg);

%% ----- Preflight: if safety enabled, auto-fit H/h to ny -----
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
  fpintf('Safety preflight skipped: %s', ME.message);
end

%% ----- Run NREP times (unique save_tag + RNG), keep per-run CSVs -----
run_tables = cell(NREP,1);
res_dirs   = strings(NREP,1);

for r = 1:NREP
    cfg_r = cfg;
    cfg_r.io.save_tag = sprintf('%s_%s_rep%02d', base_tag, rcsi_lbl, r);
    rng(r,'twister');

    % enforce headless per-run
    cfg_r.io.plot_mode       = "offline";
    cfg_r.io.make_reach_plot = false;
    cfg_r.io.save_artifacts  = false;

    fprintf('\n=== Additive noise sweep — repetition %d/%d — tag: %s ===\n', r, NREP, cfg_r.io.save_tag);
    S = run_sweeps(cfg_r, sweep_grid);
    S = ensure_time_totals(S);

    run_tables{r} = S;

    % Remember results dir; also copy summary.csv to a unique name
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

%% ----- Aggregate across runs (mean & std by alpha_w, D, n_m, n_s, n_k) -----
Tall = table();
for r = 1:NREP
    if ~isempty(run_tables{r})
        T = run_tables{r}; T.rep = r*ones(height(T),1);
        Tall = [Tall; T]; %#ok<AGROW>
    end
end
if isempty(Tall), error('All repetitions failed—nothing to aggregate.'); end

% Grouping keys (any sweep identifiers present)
keys = {};
for key = ["alpha_w","D","n_m","n_s","n_k"]
    if ismember(key, string(Tall.Properties.VariableNames))
        keys{end+1} = char(key); %#ok<SAGROW>
    end
end

% Numeric variables to aggregate (exclude keys + rep)
isNum = varfun(@isnumeric, Tall, 'OutputFormat','uniform');
numVars = string(Tall.Properties.VariableNames(isNum));
numVars = setdiff(numVars, [string(keys), "rep"]);

Gm = groupsummary(Tall, keys, "mean", cellstr(numVars));
Gs = groupsummary(Tall, keys, "std",  cellstr(numVars));

AGG = outerjoin(Gm, Gs, 'Keys', keys, 'MergeKeys', true, 'Type','left');
for v = numVars
    mn = "mean_" + v;  sd = "std_" + v;
    if ismember(mn, string(AGG.Properties.VariableNames))
        AGG.(v + "_mean") = AGG.(mn);  AGG.(mn) = [];
    end
    if ismember(sd, string(AGG.Properties.VariableNames))
        AGG.(v + "_std")  = AGG.(sd);  AGG.(sd) = [];
    end
end

[ddra_lbl, gray_lbl, hdr] = legend_bits_for_batch(cfg, sweep_grid, 'alpha_w', rcsi_lbl);
colors = struct('ddra',[0.23 0.49 0.77],'gray',[0.85 0.33 0.10]);

agg_dir = fullfile(cfg.io.base_dir, 'experiments','results','data', ...
                   sprintf('%s_%s_agg', base_tag, rcsi_lbl));
if ~exist(agg_dir,'dir'), mkdir(agg_dir); end
agg_csv = fullfile(agg_dir, 'summary_agg.csv');
writetable(AGG, agg_csv);
fprintf('\nWrote aggregated CSV: %s\n', agg_csv);

%% ----- Optional aggregate plots (safe-guarded; data already saved) -----
if do_plots
try
    x = AGG.alpha_w;
    M = @(name) AGG.(name + "_mean");
    S = @(name) AGG.(name + "_std"); % (std unused in lines; kept for extensions)

    % Fidelity / Conservatism vs alpha_W
    figure('Color','w','Name','Fidelity & Conservatism (mean \pm std)'); 
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
    
    try
        saveas(gcf, fullfile(agg_dir, 'fidcons_vs_alphaW_mean.png'));
        exportgraphics(gcf, fullfile(agg_dir, 'fidcons_vs_alphaW_mean.pdf'),'ContentType','vector');
    catch, end


    % Shape-aware (if present)
    hasDir = ismember("dir_eps_med_mean", string(AGG.Properties.VariableNames)) || ...
             ismember("hout_p90_mean",    string(AGG.Properties.VariableNames)) || ...
             ismember("haus_sym_med_mean",string(AGG.Properties.VariableNames));
    if hasDir
        figure('Color','w','Name','Shape-aware (mean \pm std)');
        tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

        if ismember("dir_eps_med_mean", string(AGG.Properties.VariableNames))
            nexttile; hold on; grid on; title('dir\_eps median');
            plot(x, M("dir_eps_med"), '-o','LineWidth',1.6); yline(1,'--');
            xlabel('\alpha_W'); ylabel('ratio');
        end
        if ismember("hout_p90_mean", string(AGG.Properties.VariableNames))
            nexttile; hold on; grid on; title('H_{out} p90');
            plot(x, M("hout_p90"), '-o','LineWidth',1.6);
            xlabel('\alpha_W'); ylabel('gap');
        end
        if ismember("haus_sym_med_mean", string(AGG.Properties.VariableNames))
            nexttile; hold on; grid on; title('Hausdorff (sym) median');
            plot(x, M("haus_sym_med"), '-o','LineWidth',1.6);
            xlabel('\alpha_W'); ylabel('distance');
        end

        try
            saveas(gcf, fullfile(agg_dir, 'shapeaware_vs_alphaW_mean.png'));
            exportgraphics(gcf, fullfile(agg_dir, 'shapeaware_vs_alphaW_mean.pdf'),'ContentType','vector');
        catch, end
    end

    fprintf('Aggregate plots saved under: %s\n', agg_dir);
catch ME
    fprintf('Aggregate plotting failed: %s\n(Your CSVs are already saved at %s)', ME.message, agg_dir);
end
end
