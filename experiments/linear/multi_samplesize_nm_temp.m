% Headless multi-run driver: run trajectory-length sweep NREP times,
% aggregate CSVs, then (optionally) plot mean±std. Plotting is inside try/catch.

clear; clc;

NREP = 5;                           % number of repetitions
base_tag = 'kMSD_sample_size_nominalcount_sweep';
do_plots = true;                     % aggregate plots at the end (safe-guarded)

rng(1,'twister');                    % base seed; we still vary it per run

%% ----- Build a base cfg + sweep grid (no plotting inside each run) -----
cfg = struct();
cfg.io = struct('save_tag', base_tag, 'plot_mode', "offline", ...
                'make_reach_plot', false, 'save_artifacts', false);

% System & shared options
cfg.shared = struct();
cfg.shared.dyn  = "k-Mass-SD";
cfg.shared.type = "standard";
cfg.shared.p_extr = 0.3;

% Reach options
cfg.shared.options_reach = struct('zonotopeOrder',100,'verbose',false);

% Conformance (Gray)
cfg.shared.cs_base = struct('robustnessMargin',1e-9,'verbose',false, ...
                            'cost',"interval",'constraints',"half");

% Data budgets
cfg.shared.n_m = 10;   % train blocks
cfg.shared.n_s = 5;
cfg.shared.n_k = 10;

cfg.shared.n_m_val = 5;  % val blocks
cfg.shared.n_s_val = 5;
cfg.shared.n_k_val = 5;

% DDRA
cfg.ddra = struct('eta_w',1,'alpha_w',1.50, ...
                  'allow_ridge',true,'lambda',1e-8,'ridge_gamma',1.0,'ridge_policy',"MAB");

% Gray
cfg.gray = struct('methodsGray',["graySeq"]);

% Noise policy + metrics
cfg.shared.noise_for_gray = true;
cfg.shared.noise_for_ddra = true;
cfg.shared.use_noise      = true;
cfg.metrics.enhanced      = true;     % keep the richer metrics
cfg.metrics.directional   = struct('enable', true, 'Nd', 64, 'seed', 12345);

% Safety metric (we will auto-fix H/h dims below if needed)
cfg.metrics.safety = struct( ...
  'enable', false, ...
  'H', [], ...                
  'h', [], ...
  'k_set', [], ...
  'mode', 'per_block', ...
  'tau_sweep', linspace(0,0.2,21), ...
  'reach_reduce', 80);

% Memory/IO toggles (stream CSVs; don't keep sets)
cfg.lowmem = struct('gray_check_contain', true, ...
                    'store_ddra_sets',   false, ...
                    'append_csv',        true, ...
                    'zonotopeOrder_cap', 25);

cfg.allow_parallel = false;
cfg.io.base_dir = fileparts(fileparts(mfilename('fullpath'))); 

% Label used in folder names
rcsi_lbl = rcsi_label_from_cfg(cfg);

% Sweep: vary trajectory length n_k (headless)
sweep_grid = struct();
sweep_grid.D_list        = 2;
sweep_grid.alpha_w_list  = cfg.ddra.alpha_w;
sweep_grid.n_m_list      = 4:2:20;          
sweep_grid.n_s_list      = 2;           
sweep_grid.n_k_list      = 10;    
sweep_grid.pe_list       = { struct('mode','randn','order',4,'strength',1,'deterministic',true) };

%% ----- Preflight: fill safety H/h to match ny if enabled -----
try
    if isfield(cfg,'metrics') && isfield(cfg.metrics,'safety') && cfg.metrics.safety.enable
        % Build one system to learn ny
        [axes0, baseC0] = init_sweep_axes(cfg, sweep_grid);
        [sys_cora, ~, ~, ~] = build_true_system(baseC0);
        ny = size(sys_cora.C,1);
        H = getfield(cfg.metrics.safety,'H',[]);
        h = getfield(cfg.metrics.safety,'h',[]);
        if isempty(H) || size(H,2)~=ny || isempty(h) || size(h,1)~=size(H,1)
            % Axis-aligned box: |y_i| <= 10 for all outputs
            cfg.metrics.safety.H = [ eye(ny); -eye(ny) ];
            cfg.metrics.safety.h = [ 10*ones(ny,1); 10*ones(ny,1) ];
        end
    end
catch ME
    fprintf('Preflight safety H/h setup skipped: %s', ME.message);
end

%% ----- Run NREP times (unique save_tag & RNG each), collect results -----
run_tables = cell(NREP,1);
res_dirs   = strings(NREP,1);

for r = 1:NREP
    cfg_r = cfg;
    cfg_r.io.save_tag = sprintf('%s_%s_rep%02d', base_tag, rcsi_lbl, r);
    rng(r,'twister');                               % different seed each rep

    % IMPORTANT: ensure no plotting during runs
    cfg_r.io.plot_mode      = "offline";
    cfg_r.io.make_reach_plot= false;
    cfg_r.io.save_artifacts = false;

    fprintf('\n=== Running repetition %d/%d — tag: %s ===\n', r, NREP, cfg_r.io.save_tag);
    S = run_sweeps(cfg_r, sweep_grid);
    S = ensure_time_totals(S);

    run_tables{r} = S;

    % Remember where summary.csv was written (same logic as run_sweeps)
    [~, rd] = init_io(cfg_r);
    res_dirs(r) = string(rd);

    % Convenience copy with explicit name (optional)
    try
        src = fullfile(rd,'summary.csv');
        if isfile(src)
            copyfile(src, fullfile(rd, sprintf('summary_%s.csv', cfg_r.io.save_tag)));
        end
    catch ME
        warning('Could not copy per-run CSV for rep %d: %s', r, ME.message);
    end
    % caps = write_captions_from_cfg(cfg);
end

%% ----- Aggregate across runs (mean & std by n_k, n_m, n_s, D, alpha_w) -----
% Stack and add a 'rep' column
Tall = table();
for r = 1:NREP
    if ~isempty(run_tables{r})
        T = run_tables{r};
        T.rep = r*ones(height(T),1);
        Tall = [Tall; T]; 
    end
end
if isempty(Tall)
    error('All repetitions failed — nothing to aggregate.');
end

% Grouping keys (keep all sweep identifiers that might vary)
keys = {};
for key = ["n_k","n_m","n_s","D","alpha_w"]
    if ismember(key, string(Tall.Properties.VariableNames))
        keys{end+1} = char(key); 
    end
end

% Numeric variables to aggregate (exclude the keys and 'rep')
isNum = varfun(@isnumeric, Tall, 'OutputFormat','uniform');
numVars = string(Tall.Properties.VariableNames(isNum));
numVars = setdiff(numVars, [string(keys), "rep"]);

Gm = groupsummary(Tall, keys, "mean", cellstr(numVars));
Gs = groupsummary(Tall, keys, "std",  cellstr(numVars));

% Merge mean & std; simplify names (mean_x -> x_mean, std_x -> x_std)
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

[ddra_lbl, gray_lbl, hdr] = legend_bits_for_batch(cfg, sweep_grid, 'n_m', rcsi_lbl);
colors = struct('ddra',[0.23 0.49 0.77],'gray',[0.85 0.33 0.10]);

% Write aggregated CSV before any plotting 
agg_dir = fullfile(cfg.io.base_dir, 'experiments','results','data', ...
                   sprintf('%s_%s_agg', base_tag, rcsi_lbl));
if ~exist(agg_dir,'dir'), mkdir(agg_dir); end
agg_csv = fullfile(agg_dir, 'summary_agg.csv');
writetable(AGG, agg_csv);
fprintf('\nWrote aggregated CSV: %s\n', agg_csv);

%% ----- Optional aggregate plots (safe-guarded) -----
if do_plots
try
    % Expect n_m to vary
    x = AGG.n_m;
    % Helpers to pull mean/std safely
    M = @(name) AGG.(name + "_mean");
    S = @(name) AGG.(name + "_std");

    % a) Runtime panels
    figure('Color','w','Name','Runtime (mean \pm std)');
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
    
    nexttile; hold on; grid on;
    plot(x, M("t_ddra_total"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("t_gray_total"), '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('n_m'); ylabel('seconds'); title(sprintf('%s — Total runtime vs n_m', hdr)); legend('Location','best');
    
    nexttile; hold on; grid on;
    plot(x, M("t_ddra_learn"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("t_gray_learn"), '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('n_m'); ylabel('seconds'); title(sprintf('%s — Learning runtime', hdr)); legend('Location','best');
    
    nexttile; hold on; grid on;
    plot(x, M("t_ddra_check"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("t_gray_val"),   '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('n_m'); ylabel('seconds'); title(sprintf('%s — Validation/Check', hdr)); legend('Location','best');
    
    nexttile; hold on; grid on;
    plot(x, M("t_ddra_infer"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("t_gray_infer"), '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('n_m'); ylabel('seconds'); title(sprintf('%s — Inference', hdr)); legend('Location','best');
    
    try
        saveas(gcf, fullfile(agg_dir,'runtime_panels_mean.png'));
        exportgraphics(gcf, fullfile(agg_dir,'runtime_panels_mean.pdf'),'ContentType','vector');
    catch, end
    
    % b) Fidelity / Conservatism
    figure('Color','w','Name','Fidelity & Conservatism (mean \pm std)');
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    
    nexttile; hold on; grid on;
    plot(x, M("cval_ddra"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("cval_gray"), '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('n_m'); ylabel('containment (%)'); title(sprintf('%s — Fidelity vs n_m', hdr)); legend('Location','best');
    
    nexttile; hold on; grid on;
    plot(x, M("sizeI_ddra"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("sizeI_gray"), '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('n_m'); ylabel('\Sigma interval widths'); title(sprintf('%s — Conservatism vs n_m', hdr)); legend('Location','best');
    
    try
        saveas(gcf, fullfile(agg_dir,'fidelity_conservatism_mean.png'));
        exportgraphics(gcf, fullfile(agg_dir,'fidelity_conservatism_mean.pdf'),'ContentType','vector');
    catch, end


    fprintf('Aggregate plots saved under: %s\n', agg_dir);
catch ME
    fprintf('Aggregate plotting failed: %s\n(Your CSVs are already saved at %s)', ME.message, agg_dir);
end
end
