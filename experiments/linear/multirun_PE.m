%% Persistency-of-Excitation Sweep — Multi-run batch (DDRA vs RCSI/Gray)
% - Runs NREP repetitions with unique save_tag and RNG
% - Shared datasets per run; unified noise policy
% - Aggregates across reps (mean/std) with grouping by pe_mode & pe_order
% - Optional plots: fidelity & conservatism vs PE order (per PE mode)
% - Optional TikZ export for LaTeX

clear; clc;

NREP     = 5;                              % how many repetitions
base_tag = 'kMSD_pe_sweep';                % base save_tag prefix
do_plots = true;                           % aggregate plots at the end

rng(1,'twister');                          % base seed; per-run we set rng(r)

%% ----- Build base cfg + sweep grid (no plotting during each run) -----
cfg = struct();
cfg.io = struct('save_tag', base_tag, ...
                'plot_mode', "offline", ...
                'make_reach_plot', false, ...
                'save_artifacts', false);

% --- LaTeX/TikZ export toggles (optional) ---
cfg.io.export_tikz        = true;                    % set false to disable
cfg.io.latex_interpreter  = true;                    % LaTeX tick/legend/text in MATLAB
cfg.io.tikz               = struct('width','\figW','height','\figH');  % PGFPlots size macros

% Shared system/options
cfg.shared = struct();
cfg.shared.dyn   = "k-Mass-SD";
cfg.shared.type  = "standard";
cfg.shared.p_extr = 0.3;
cfg.shared.options_reach = struct('zonotopeOrder',100,'verbose',false);
cfg.shared.cs_base = struct('robustnessMargin',1e-9,'verbose',false, ...
                            'cost',"interval",'constraints',"half");
cfg.shared.cs_base.derivRecomputation = false;  % linear systems: skip re-derives

% Data budgets (train/val)
cfg.shared.n_m = 10; 
cfg.shared.n_s = 5;
cfg.shared.n_k = 20;
cfg.shared.n_m_val = 5;
cfg.shared.n_s_val = 5;  
cfg.shared.n_k_val = 5;

% Unified noise policy
cfg.shared.noise_for_gray = true;
cfg.shared.noise_for_ddra = true;
cfg.shared.use_noise      = true;

% DDRA (incl. ridge guard & variant tag for legends)
cfg.ddra = struct('eta_w',1,'alpha_w',0.01, ...
                  'allow_ridge',false,'lambda',1e-8, ...
                  'ridge_gamma',1.0,'ridge_policy',"MAB", ...
                  'variant',"std");        % "std" or "meas"

% Gray methods
cfg.gray = struct('methodsGray', ["graySeq"]);

% Metrics (lean; speed)
cfg.metrics = struct();
cfg.metrics.enhanced = true;
cfg.metrics.directional = struct('enable', false, 'Nd', 64, 'seed', 12345);
cfg.metrics.safety = struct('enable', false, 'H', [], 'h', [], ...
                            'k_set', [], 'mode', 'per_block', ...
                            'tau_sweep', linspace(0,0.2,21), ...
                            'reach_reduce', 80);

% Efficiency toggles
cfg.lowmem = struct('gray_check_contain', true, ...
                    'store_ddra_sets',   false, ...
                    'append_csv',        true, ...
                    'zonotopeOrder_cap', 25);

cfg.allow_parallel = false;
cfg.io.base_dir    = fileparts(fileparts(mfilename('fullpath'))); 

% Label used in folder names
rcsi_lbl = rcsi_label_from_cfg(cfg);

% --- Sweep grid: choose your PE orders & shapes here ---
PE_orders = 1:8;                                  
% Two example shapes: 'randn' and 'sinwavecs' 
pe_randn = arrayfun(@(L) struct('mode','randn',    'order',L,'strength',1,  'deterministic',true), PE_orders, 'uni',0);
pe_sin   = arrayfun(@(L) struct('mode','sinwavecs','order',L,'strength',100,'deterministic',true), PE_orders, 'uni',0);

sweep_grid = struct();
sweep_grid.D_list        = 2;
sweep_grid.alpha_w_list  = cfg.ddra.alpha_w;
sweep_grid.n_m_list      = cfg.shared.n_m;
sweep_grid.n_s_list      = cfg.shared.n_s;
sweep_grid.n_k_list      = cfg.shared.n_k;
sweep_grid.pe_list = arrayfun(@(L) struct('mode','sinwave', ...
    'order',L, 'strength',100, 'deterministic',true), PE_orders, 'uni',0);

% Header/legend color palette
colors = struct('ddra',[0.23 0.49 0.77], 'gray',[0.85 0.33 0.10]);

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
  fprintf('Safety preflight skipped: %s\n', ME.message);
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

    fprintf('\n=== PE sweep — repetition %d/%d — tag: %s ===\n', r, NREP, cfg_r.io.save_tag);
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

%% ----- Aggregate across runs (mean & std by pe_mode, pe_order, D, n_m, n_s, n_k, alpha_w) -----
Tall = table();
for r = 1:NREP
    if ~isempty(run_tables{r})
        T = run_tables{r}; T.rep = r*ones(height(T),1);
        Tall = [Tall; T]; 
    end
end
if isempty(Tall), error('All repetitions failed—nothing to aggregate.'); end

% Grouping keys (whatever is present)
keys = {};
for key = ["pe_mode","pe_order","D","n_m","n_s","n_k","alpha_w"]
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

% Flatten mean/std names: mean_x -> x_mean, std_x -> x_std
for v = numVars
    mn = "mean_" + v;  sd = "std_" + v;
    if ismember(mn, string(AGG.Properties.VariableNames))
        AGG.(v + "_mean") = AGG.(mn);  AGG.(mn) = [];
    end
    if ismember(sd, string(AGG.Properties.VariableNames))
        AGG.(v + "_std")  = AGG.(sd);  AGG.(sd) = [];
    end
end

% Legend strings + header (exclude the swept axis 'pe')
try
    [ddra_lbl, gray_lbl, hdr] = legend_bits_for_batch(cfg, sweep_grid, 'pe', rcsi_lbl);
catch
    % Fallback header if the helper isn't available
    sys_name = char(cfg.shared.dyn);
    D = getfield(sweep_grid,'D_list',2); if numel(D)>1, D=D(1); end
    ddra_lbl = sprintf('DDRA (variant=%s)', string(getfield(cfg.ddra,'variant',"std")));
    gray_lbl = sprintf('RCSI-%s', rcsi_lbl);
    hdr = sprintf('%s (D=%g, n_m=%g, n_s=%g, n_k=%g)', sys_name, D, ...
        sweep_grid.n_m_list(1), sweep_grid.n_s_list(1), sweep_grid.n_k_list(1));
end

% Write aggregated CSV before plotting
agg_dir = fullfile(cfg.io.base_dir, 'experiments','results','data', ...
                   sprintf('%s_%s_agg', base_tag, rcsi_lbl));
if ~exist(agg_dir,'dir'), mkdir(agg_dir); end
agg_csv = fullfile(agg_dir, 'summary_agg.csv');
writetable(AGG, agg_csv);
fprintf('\nWrote aggregated CSV: %s\n', agg_csv);

%% ----- Optional aggregate plots (safe-guarded; data already saved) -----
if do_plots
try
    % LaTeX interpreters (optional)
    if isfield(cfg,'io') && getfield(cfg.io,'latex_interpreter', false)
        set(groot,'defaultTextInterpreter','latex');
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
    end

    M = @(name) AGG.(name + "_mean");

    % Prepare modes present
    modes = unique(string(AGG.pe_mode));
    if isempty(modes), error('AGG has no pe_mode column to plot.'); end

    % Marker set per mode (cycled)
    mk = {'o','s','^','d','v','>','<','p','h','x','+'};

    % ----- Fidelity & Conservatism vs PE order (per mode) -----
    figure('Color','w','Name','Fidelity & Conservatism (mean \pm std)');
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % (1) Fidelity
    nexttile; hold on; grid on;
    title(sprintf('%s — Fidelity vs $L$', hdr),'Interpreter','latex');
    xlabel('$L$','Interpreter','latex'); ylabel('Containment on validation (%)');
    leg1 = {};
    for i = 1:numel(modes)
        m = modes(i);
        mask = string(AGG.pe_mode)==m;
        [ord,idx] = sort(double(AGG.pe_order(mask)));
        xL  = ord;
        yDd = M("cval_ddra"); yDd = yDd(mask); yDd = yDd(idx);
        yGr = M("cval_gray"); yGr = yGr(mask); yGr = yGr(idx);

        plot(xL, yDd, ['-' mk{mod(i-1,numel(mk))+1}], 'LineWidth',1.6, ...
            'Color',colors.ddra, 'DisplayName', sprintf('%s (%s)', ddra_lbl, m));
        plot(xL, yGr, ['-' mk{mod(i-1,numel(mk))+1}], 'LineWidth',1.6, ...
            'Color',colors.gray, 'DisplayName', sprintf('%s (%s)', gray_lbl, m));
    end
    legend('Location','best');

    % (2) Conservatism
    nexttile; hold on; grid on;
    title(sprintf('%s — Conservatism vs $L$', hdr),'Interpreter','latex');
    xlabel('$L$','Interpreter','latex'); ylabel('\Sigma interval widths');
    for i = 1:numel(modes)
        m = modes(i);
        mask = string(AGG.pe_mode)==m;
        [ord,idx] = sort(double(AGG.pe_order(mask)));
        xL  = ord;
        yDd = M("sizeI_ddra"); yDd = yDd(mask); yDd = yDd(idx);
        yGr = M("sizeI_gray"); yGr = yGr(mask); yGr = yGr(idx);

        plot(xL, yDd, ['-' mk{mod(i-1,numel(mk))+1}], 'LineWidth',1.6, ...
            'Color',colors.ddra, 'DisplayName', sprintf('%s (%s)', ddra_lbl, m));
        plot(xL, yGr, ['-' mk{mod(i-1,numel(mk))+1}], 'LineWidth',1.6, ...
            'Color',colors.gray, 'DisplayName', sprintf('%s (%s)', gray_lbl, m));
    end
    legend('Location','best');

    % Save figures + optional TikZ
    try
        saveas(gcf, fullfile(agg_dir,'fidcons_vs_L_mean.png'));
        exportgraphics(gcf, fullfile(agg_dir,'fidcons_vs_L_mean.pdf'),'ContentType','vector');
    catch, end
    maybe_export_tikz(gcf, fullfile(agg_dir,'fidcons_vs_L_mean.tex'), cfg);

    fprintf('Aggregate plots saved under: %s\n', agg_dir);
catch ME
    fprintf('Aggregate plotting failed: %s\n(Your CSVs are already saved at %s)\n', ME.message, agg_dir);
end
end

%% ----- Helpers: TikZ export (optional) -----
function maybe_export_tikz(fig, outfile, cfg)
% Export current figure to PGFPlots/TikZ if matlab2tikz is available.
    try
        if ~isfield(cfg,'io') || ~isfield(cfg.io,'export_tikz') || ~cfg.io.export_tikz
            return;
        end
        hasM2T = exist('matlab2tikz','file')==2;
        if ~hasM2T, return; end

        try, cleanfigure; catch, end  % optional, if available

        width  = get_tikz_dim(cfg,'width','\textwidth');
        height = get_tikz_dim(cfg,'height','0.6\textwidth');

        matlab2tikz('filename', outfile, ...
            'showInfo', false, ...
            'standalone', false, ...
            'width',  width, ...
            'height', height, ...
            'interpretTickLabelsAsTex', true);

        fprintf('TikZ exported: %s\n', outfile);
    catch ME
        fprintf('TikZ export failed: %s', ME.message);
    end
end

function v = get_tikz_dim(cfg, key, defaultVal)
    v = defaultVal;
    if isfield(cfg,'io') && isfield(cfg.io,'tikz') && isfield(cfg.io.tikz, key)
        val = cfg.io.tikz.(key);
        if ~isempty(val), v = val; end
    end
end
