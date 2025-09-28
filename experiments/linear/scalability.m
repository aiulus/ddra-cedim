%% HOW TO USE — State-Dimension Scalability (DDRA vs RCSI/Gray)
rng(1,'twister');

% ---------- Minimal cfg ----------
cfg = struct();
cfg.io = struct('save_tag', 'kMSD_dim_sweep');

% System & shared options
cfg.shared = struct();
cfg.shared.dyn   = "k-Mass-SD";
cfg.shared.type  = "standard";
cfg.shared.p_extr = 0.3;

% Identification + reachability options 
cfg.shared.options_reach = struct('zonotopeOrder',100,'verbose',false);

% Conformance (Gray) base options
cfg.shared.cs_base = struct( ...
    'robustnessMargin', 1e-9, ...
    'verbose', false, ...
    'cost', "interval", ...
    'constraints', "half");

% Data budgets (fixed)
cfg.shared.n_m     = 3;
cfg.shared.n_s     = 2;
cfg.shared.n_k     = 20;
cfg.shared.n_m_val = 1;
cfg.shared.n_s_val = cfg.shared.n_s;
cfg.shared.n_k_val = cfg.shared.n_k;

% DDRA noise setup (fixed)
cfg.ddra = struct();
cfg.ddra.eta_w   = 1;
cfg.ddra.alpha_w = 0.01;
cfg.ddra.allow_ridge   = false;
cfg.ddra.lambda        = 1e-8;
cfg.ddra.ridge_gamma   = 1.0;
cfg.ddra.ridge_policy  = "MAB";

% Gray methods (keep simple/fast)
cfg.gray = struct();
cfg.gray.methodsGray = ["grayLS"];

rcsi_lbl = rcsi_label_from_cfg(cfg);                
cfg.io.save_tag = sprintf('%s_%s', cfg.io.save_tag, rcsi_lbl);  

% ---------- Sweep grid (vary only D) ----------
sweep_grid = struct();
sweep_grid.D_list       = 2:2:10;     
sweep_grid.alpha_w_list = cfg.ddra.alpha_w;
sweep_grid.n_m_list     = cfg.shared.n_m;
sweep_grid.n_s_list     = cfg.shared.n_s;
sweep_grid.n_k_list     = cfg.shared.n_k;
sweep_grid.pe_list      = { struct('mode','randn','strength',1,'deterministic',true) };

% Memory efficiency toggles
cfg.lowmem = struct();
cfg.lowmem.gray_check_contain = false;
cfg.lowmem.store_ddra_sets    = false;
cfg.lowmem.append_csv         = true;
cfg.lowmem.zonotopeOrder_cap  = 50;

% Unified noise policy for this study (no noise by default)
cfg.shared.noise_for_gray = false;   % Gray/RCSI W=0
cfg.shared.noise_for_ddra = true;    % DDRA W≠0 unless you set false
cfg.shared.use_noise      = false;   % overall

% plotting / IO
cfg.io.plot_mode        = "offline";
cfg.io.make_reach_plot  = true;
cfg.io.plot_rows        = [1];
cfg.io.plot_dims        = [1 2];
cfg.io.plot_every_k     = 1;
cfg.io.save_artifacts   = false;     % set true if you want artifact-derived plots
cfg.metrics             = struct('enhanced', false);
cfg.io.base_dir         = fileparts(fileparts(mfilename('fullpath')));
cfg.allow_parallel      = false;
cfg.shared.pe_min_policy = 'none';
cfg.shared.pe_verbose    = true;

%% ---------- Run ----------
SUMMARY = run_sweeps(cfg, sweep_grid);
SUMMARY = ensure_time_totals(SUMMARY);

%% ---------- Axis selection (here: D varies) ----------
var_axis = 'D';
x        = coerce_numeric(SUMMARY.D);
xlab     = 'Dimension D';
hdr      = sprintf('%s — Dimensional scalability', char(cfg.shared.dyn));

% Colors & labels
colors = struct('ddra',[0.23 0.49 0.77],'gray',[0.85 0.33 0.10]);
ddra_lbl = 'DDRA';
rcsi_lbl_full = 'RCSI-Gray';

%% ---------- Runtime panels ----------
[plots_dir, results_dir] = init_io(cfg);

f = figure('Name','Runtime panels','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_total, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName',ddra_lbl);
plot(x, SUMMARY.t_gray_total, '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName',rcsi_lbl_full);
xlabel(xlab); ylabel('Seconds'); title([hdr ' — Total']); legend('Location','best');

nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_learn, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName',ddra_lbl);
plot(x, SUMMARY.t_gray_learn, '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName',rcsi_lbl_full);
xlabel(xlab); ylabel('Seconds'); title('Learning'); legend('Location','best');

nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_check, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName',ddra_lbl);
plot(x, SUMMARY.t_gray_val,   '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName',rcsi_lbl_full);
xlabel(xlab); ylabel('Seconds'); title('Validation/Check'); legend('Location','best');

nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_infer, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName',ddra_lbl);
plot(x, SUMMARY.t_gray_infer, '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName',rcsi_lbl_full);
xlabel(xlab); ylabel('Seconds'); title('Inference'); legend('Location','best');

save_plot(f, plots_dir, sprintf('runtime_panels_vs_%s_%s', var_axis, rcsi_lbl), ...
    'Formats', {'png','pdf'}, 'Resolution', 200);

%% ---------- Fidelity / Conservatism ----------
cval_ddra  = coerce_numeric(SUMMARY.cval_ddra);
cval_gray  = coerce_numeric(SUMMARY.cval_gray);
sizeI_ddra = coerce_numeric(SUMMARY.sizeI_ddra);
sizeI_gray = coerce_numeric(SUMMARY.sizeI_gray);

f2 = figure('Name','Fidelity & Conservatism','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
plot(x, cval_ddra, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName',ddra_lbl);
plot(x, cval_gray, '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName',rcsi_lbl_full);
xlabel(xlab); ylabel('Containment on validation (%)');
title([hdr ' — Fidelity']); legend('Location','best');

nexttile; hold on; grid on;
plot(x, sizeI_ddra, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName',ddra_lbl);
plot(x, sizeI_gray, '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName',rcsi_lbl_full);
xlabel(xlab); ylabel('Aggregated interval size'); 
title([hdr ' — Conservatism']); legend('Location','best');

save_plot(f2, plots_dir, sprintf('fidelity_conservatism_vs_%s_%s', var_axis, rcsi_lbl), ...
    'Formats', {'png','pdf'}, 'Resolution', 200);

%% ---------- Shape-aware metrics across sweep (if present) ----------
dir_med  = coerce_numeric(getcol(SUMMARY,'dir_eps_med'));
dir_p90  = coerce_numeric(getcol(SUMMARY,'dir_eps_p90'));
hout_p90 = coerce_numeric(getcol(SUMMARY,'hout_p90'));
haus_med = coerce_numeric(getcol(SUMMARY,'haus_sym_med'));
mw_g     = coerce_numeric(getcol(SUMMARY,'mw_gray_mean'));
mw_d     = coerce_numeric(getcol(SUMMARY,'mw_ddra_mean'));

f3 = figure('Name','Shape-aware metrics','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
plot(x, dir_med, '-o', 'LineWidth',1.6, 'DisplayName','dir\_eps median');
plot(x, dir_p90, '-s', 'LineWidth',1.6, 'DisplayName','dir\_eps p90');
yline(1,'--'); xlabel(xlab); ylabel('ratio'); 
title([hdr ' — Directional ratio (Gray/True)']); legend('Location','best');

nexttile; hold on; grid on;
plot(x, hout_p90, '-o', 'LineWidth',1.6, 'DisplayName','H_{out} p90');
xlabel(xlab); ylabel('support gap'); 
title([hdr ' — H_{out} (p90)']); legend('Location','best');

nexttile; hold on; grid on;
plot(x, haus_med, '-o', 'LineWidth',1.6, 'DisplayName','Hausdorff (sym) median');
xlabel(xlab); ylabel('distance'); 
title([hdr ' — Symmetric Hausdorff']); legend('Location','best');

nexttile; hold on; grid on;
plot(x, mw_g, '-s', 'LineWidth',1.6, 'DisplayName','Gray mean width');
plot(x, mw_d, '-o', 'LineWidth',1.6, 'DisplayName','DDRA mean width');
xlabel(xlab); ylabel('mean width'); 
title([hdr ' — Mean width (PRE-output)']); legend('Location','best');

save_plot(f3, plots_dir, sprintf('shape_aware_vs_%s_%s', var_axis, rcsi_lbl), ...
    'Formats', {'png','pdf'}, 'Resolution', 200);

%% ---------- Per-step panels from CSV (if available) ----------
perstep_path = fullfile(results_dir, 'summary_perstep.csv');
if exist(perstep_path,'file')
    P = readtable(perstep_path);
    SUMMARY.row = (1:height(SUMMARY))';
    [~,ix] = max(x); 
    row_pick = SUMMARY.row(ix);
    K = P(P.row==row_pick, :);

    f4 = figure('Name','Per-step width, coverage, and ratio','Color','w');
    tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

    nexttile; hold on; grid on;
    plot(K.k, K.wid_ddra, '-o', 'LineWidth',1.6, 'DisplayName','DDRA');
    plot(K.k, K.wid_gray, '-s', 'LineWidth',1.6, 'DisplayName','Gray');
    xlabel('step k'); ylabel('interval width (sum)'); 
    title([hdr sprintf(' — Per-step width @ %s=%g', var_axis, x(ix))]); legend('Location','best');

    nexttile; hold on; grid on;
    plot(K.k, K.cov_ddra, '-o', 'LineWidth',1.6, 'DisplayName','DDRA');
    plot(K.k, K.cov_gray, '-s', 'LineWidth',1.6, 'DisplayName','Gray');
    xlabel('step k'); ylabel('containment (%)'); 
    title([hdr sprintf(' — Per-step coverage @ %s=%g', var_axis, x(ix))]); legend('Location','best');

    nexttile; hold on; grid on;
    if any(strcmpi(K.Properties.VariableNames,'ratio_gray_true'))
        plot(K.k, K.ratio_gray_true, '-d', 'LineWidth',1.6, 'DisplayName','Gray/True');
        xlabel('step k'); ylabel('ratio'); 
        title([hdr sprintf(' — Per-step size ratio @ %s=%g', var_axis, x(ix))]);
        yline(1,'--'); legend('Location','best');
    else
        text(0.5,0.5,'ratio\_gray\_true not available','HorizontalAlignment','center'); axis off
    end

    save_plot(f4, plots_dir, sprintf('perstep_row%03d_%s_%s', row_pick, var_axis, rcsi_lbl), ...
        'Formats', {'png','pdf'}, 'Resolution', 200);
end

%% ---------- Artifact-derived helper plots (if artifacts exist) ----------
artifact_dir = fullfile(results_dir, 'artifacts');
if exist(artifact_dir,'dir')
    SUMMARY.row = (1:height(SUMMARY))';
    [~,ix] = max(x);
    row_pick = SUMMARY.row(ix);
    art_path = find_artifact_by_row(artifact_dir, row_pick);
    if ~isempty(art_path)
        savebase = fullfile(plots_dir, sprintf('row%04d_artifacts_%s_%s', row_pick, var_axis, rcsi_lbl));
        plot_perstep_fidelity_sizes(art_path, ...
            'Dims', cfg.io.plot_dims, ...
            'Reduce', getfielddef(cfg.lowmem,'zonotopeOrder_cap', 60), ...
            'SaveBase', savebase, 'Show', true);
        nkv = unique(coerce_numeric(SUMMARY.n_k)); 
        if numel(nkv)==1, nkv = nkv(1); else, nkv = cfg.shared.n_k_val; end
        Klist = unique([1, round(nkv/2), nkv]);
        overlay_reach_sets(art_path, ...
            'Dims', cfg.io.plot_dims, 'Klist', Klist, ...
            'Reduce', getfielddef(cfg.lowmem,'zonotopeOrder_cap', 60), ...
            'SaveBase', savebase, 'Show', true);
    end
end

close all force

% -------- helpers ----------
function v = getcol(T, name)
    vars = string(T.Properties.VariableNames);
    idx  = find(strcmpi(vars, string(name)), 1);
    if ~isempty(idx), v = T.(vars(idx));
    else, v = nan(height(T),1);
    end
end

function p = find_artifact_by_row(artifact_dir, rownum)
    p = "";
    pat = sprintf('row_%04d', rownum);
    L = dir(fullfile(artifact_dir, '*.mat'));
    for k = 1:numel(L)
        if contains(L(k).name, pat), p = string(fullfile(L(k).folder, L(k).name)); break; end
    end
    if p=="" && ~isempty(L), p = string(fullfile(L(1).folder, L(1).name)); end
end
