%% State-Dimension Scalability — Multi-run batch (DDRA vs RCSI/Gray)
% - Runs NREP repetitions with unique save_tag & RNG
% - Streams per-run CSVs (no per-run plotting/artifacts)
% - Aggregates to summary_agg.csv
% - Optional aggregate plots are wrapped in try/catch so CSVs persist

clear; clc;

NREP     = 1;                       % number of repetitions
base_tag = 'scalability_kMSD_dim_sweep';        % base save_tag
do_plots = true;                    % aggregate plots at the end (guarded)

rng(1,'twister');                   % base seed; per-run we set rng(r)

%% ----- Build base cfg + sweep grid (no plotting inside each run) -----
cfg = struct();
cfg.io = struct('save_tag', base_tag, ...
                'plot_mode', "offline", ...
                'make_reach_plot', false, ...
                'save_artifacts', false);

% Shared system options (entry-point defaults)
cfg.shared = struct();
cfg.shared.dyn   = "k-Mass-SD";
cfg.shared.type  = "standard";
cfg.shared.p_extr = 0.3;
cfg.shared.options_reach = struct('zonotopeOrder',100,'verbose',false);
cfg.shared.cs_base = struct('robustnessMargin',1e-9,'verbose',false, ...
                            'cost',"interval",'constraints',"half");

% Data budgets (fixed)
cfg.shared.n_m     = 10;
cfg.shared.n_s     = 5;
cfg.shared.n_k     = 15;
cfg.shared.n_m_val = 10;
cfg.shared.n_s_val = cfg.shared.n_s;
cfg.shared.n_k_val = cfg.shared.n_k;

% DDRA noise setup (fixed)
cfg.ddra = struct('eta_w',1, 'alpha_w',0.01, ...
                  'allow_ridge',false,'lambda',1e-8, ...
                  'ridge_gamma',1.0,'ridge_policy',"MAB");

% Gray methods (fast)
cfg.gray = struct('methodsGray', ["graySeq"]);

% Unified noise policy for this study (as in entry-point: overall off)
cfg.shared.noise_for_gray = false;   % Gray W = 0
cfg.shared.noise_for_ddra = true;    % DDRA would use W unless disabled by overall flag
cfg.shared.use_noise      = false;   % override -> both effectively W = 0

% Memory/IO
cfg.lowmem = struct('gray_check_contain', false, ...
                    'store_ddra_sets',   false, ...
                    'append_csv',        true, ...
                    'zonotopeOrder_cap', 25);

% Metrics (keep lean for scalability study)
cfg.metrics = struct('enhanced', false);

cfg.allow_parallel = false;
cfg.io.base_dir    = fileparts(fileparts(mfilename('fullpath'))); 

% Label used in folder names
rcsi_lbl = rcsi_label_from_cfg(cfg);

% Sweep grid: vary only D (others fixed)
sweep_grid = struct();
sweep_grid.D_list       = 2:2:10;     % <-- scalability axis
sweep_grid.alpha_w_list = cfg.ddra.alpha_w;
sweep_grid.n_m_list     = cfg.shared.n_m;
sweep_grid.n_s_list     = cfg.shared.n_s;
sweep_grid.n_k_list     = cfg.shared.n_k;
Lcap = cfg.shared.n_k - 2;    % with n_k=10 -> 8 (safe margin for the checker)
sweep_grid.pe_list = { struct('mode','randn','strength',1,'deterministic',true, ...
                              'order_gate', Lcap, 'force_order', true) };


%% ----- Preflight: if safety enabled, auto-fit H/h to ny (keeps dims consistent) -----
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

%% ----- Run NREP times (unique save_tag & RNG each), collect results -----
run_tables = cell(NREP,1);
res_dirs   = strings(NREP,1);

for r = 1:NREP
    cfg_r = cfg;
    cfg_r.io.save_tag = sprintf('%s_%s_rep%02d', base_tag, rcsi_lbl, r);
    rng(r,'twister');                               % different seed each rep

    % enforce headless per-run
    cfg_r.io.plot_mode      = "offline";
    cfg_r.io.make_reach_plot= false;
    cfg_r.io.save_artifacts = false;

    fprintf('\n=== Dim sweep — repetition %d/%d — tag: %s ===\n', r, NREP, cfg_r.io.save_tag);
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
end

%% ----- Aggregate across runs (mean & std by D, alpha_w, n_m, n_s, n_k) -----
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
for key = ["D","alpha_w","n_m","n_s","n_k"]
    if ismember(key, string(Tall.Properties.VariableNames))
        keys{end+1} = char(key); 
    end
end

% Numeric variables to aggregate (exclude keys + 'rep')
isNum   = varfun(@isnumeric, Tall, 'OutputFormat','uniform');
numVars = string(Tall.Properties.VariableNames(isNum));
numVars = setdiff(numVars, [string(keys), "rep"]);

Gm = groupsummary(Tall, keys, "mean", cellstr(numVars));
Gs = groupsummary(Tall, keys, "std",  cellstr(numVars));

% Merge mean & std; rename mean_x -> x_mean, std_x -> x_std
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

% Legend strings + header
[ddra_lbl, gray_lbl, hdr] = legend_bits_for_batch(cfg, sweep_grid, 'D', rcsi_lbl);
colors = struct('ddra',[0.23 0.49 0.77],'gray',[0.85 0.33 0.10]);

% Write aggregated CSV
agg_dir = fullfile(cfg.io.base_dir, 'experiments','results','data', ...
                   sprintf('%s_%s_agg', base_tag, rcsi_lbl));
if ~exist(agg_dir,'dir'), mkdir(agg_dir); end
agg_csv = fullfile(agg_dir, 'summary_agg.csv');
writetable(AGG, agg_csv);
fprintf('\nWrote aggregated CSV: %s\n', agg_csv);

%% ----- Optional aggregate plots (safe-guarded; data already saved) -----
if do_plots
try
    x = AGG.D;
    M = @(name) AGG.(name + "_mean");
    S = @(name) AGG.(name + "_std"); 

    % Runtime panels
    figure('Color','w','Name','Runtime (mean \pm std)');
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
    nexttile; hold on; grid on;
    plot(x, M("t_ddra_total"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("t_gray_total"), '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('D'); ylabel('seconds'); title(sprintf('%s — Total vs D', hdr)); legend('Location','best');

    nexttile; hold on; grid on;
    plot(x, M("t_ddra_learn"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("t_gray_learn"), '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('D'); ylabel('seconds'); title(sprintf('%s — Learning', hdr)); legend('Location','best');

    nexttile; hold on; grid on;
    plot(x, M("t_ddra_check"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("t_gray_val"),   '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('D'); ylabel('seconds'); title(sprintf('%s — Validation/Check', hdr)); legend('Location','best');

    nexttile; hold on; grid on;
    plot(x, M("t_ddra_infer"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("t_gray_infer"), '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('D'); ylabel('seconds'); title(sprintf('%s — Inference', hdr)); legend('Location','best');

    try
        saveas(gcf, fullfile(agg_dir,'runtime_panels_vs_D_mean.png'));
        exportgraphics(gcf, fullfile(agg_dir,'runtime_panels_vs_D_mean.pdf'),'ContentType','vector');
    catch, end

    % Fidelity / Conservatism
    figure('Color','w','Name','Fidelity & Conservatism (mean \pm std)');
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile; hold on; grid on;
    plot(x, M("cval_ddra"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("cval_gray"), '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('D'); ylabel('containment (%)'); title(sprintf('%s — Fidelity vs D', hdr)); legend('Location','best');

    nexttile; hold on; grid on;
    plot(x, M("sizeI_ddra"), '-o','LineWidth',1.6,'Color',colors.ddra,'DisplayName',ddra_lbl);
    plot(x, M("sizeI_gray"), '-s','LineWidth',1.6,'Color',colors.gray,'DisplayName',gray_lbl);
    xlabel('D'); ylabel('\Sigma interval widths'); title(sprintf('%s — Conservatism vs D', hdr)); legend('Location','best');

    try
        saveas(gcf, fullfile(agg_dir,'fidcons_vs_D_mean.png'));
        exportgraphics(gcf, fullfile(agg_dir,'fidcons_vs_D_mean.pdf'),'ContentType','vector');
    catch, end

    % Shape-aware (if present)
    hasDir = ismember("dir_eps_med_mean", string(AGG.Properties.VariableNames)) || ...
             ismember("hout_p90_mean",    string(AGG.Properties.VariableNames)) || ...
             ismember("haus_sym_med_mean",string(AGG.Properties.VariableNames));
    if hasDir
        figure('Color','w','Name','Shape-aware (mean \pm std)');
        tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

        if ismember("dir_eps_med_mean", string(AGG.Properties.VariableNames))
            nexttile; hold on; grid on; title('dir\_eps median');
            plot(x, M("dir_eps_med"), '-o','LineWidth',1.6); yline(1,'--');
            xlabel('D'); ylabel('ratio');
        end
        if ismember("hout_p90_mean", string(AGG.Properties.VariableNames))
            nexttile; hold on; grid on; title('H_{out} p90');
            plot(x, M("hout_p90"), '-o','LineWidth',1.6);
            xlabel('D'); ylabel('gap');
        end
        if ismember("haus_sym_med_mean", string(AGG.Properties.VariableNames))
            nexttile; hold on; grid on; title('Hausdorff (sym) median');
            plot(x, M("haus_sym_med"), '-o','LineWidth',1.6);
            xlabel('D'); ylabel('distance');
        end
        if ismember("mw_gray_mean_mean", string(AGG.Properties.VariableNames)) || ...
           ismember("mw_ddra_mean_mean", string(AGG.Properties.VariableNames))
            nexttile; hold on; grid on; title('Mean width (PRE-output)');
            if ismember("mw_gray_mean_mean", string(AGG.Properties.VariableNames))
                plot(x, M("mw_gray_mean"), '-s','LineWidth',1.6,'DisplayName','Gray');
            end
            if ismember("mw_ddra_mean_mean", string(AGG.Properties.VariableNames))
                plot(x, M("mw_ddra_mean"), '-o','LineWidth',1.6,'DisplayName','DDRA');
            end
            xlabel('D'); ylabel('mean width'); legend('Location','best');
        end

        try
            saveas(gcf, fullfile(agg_dir,'shapeaware_vs_D_mean.png'));
            exportgraphics(gcf, fullfile(agg_dir,'shapeaware_vs_D_mean.pdf'),'ContentType','vector');
        catch, end
    end

    fprintf('Aggregate plots saved under: %s\n', agg_dir);
catch ME
    fprintf('Aggregate plotting failed: %s\n(Your CSVs are already saved at %s)\n', ME.message, agg_dir);
end
end
