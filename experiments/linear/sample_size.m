%% HOW TO USE (Sample-Size Sweep)
% What it does:
%   Varies the number of unique input trajectories n_m (and can vary n_s, n_k)
%   at fixed dimension, producing (i) runtime panels and (ii) fidelity/conservatism plots.
%
% Key knobs:
%   - sweep_grid.n_m_list (e.g., [2 4 8 16 32])
%   - Fix dimension via sweep_grid.D_list = <D>;
%   - Keep excitation fixed: sweep_grid.pe_list = {struct('mode','randn')};
%
% Memory-efficiency toggles:
%   cfg.lowmem.gray_check_contain = false;
%   cfg.lowmem.store_ddra_sets    = false;
%   cfg.lowmem.append_csv         = true;
%   cfg.lowmem.zonotopeOrder_cap  = 50;
%
% Outputs:
%   - CSV + plots stored via init_io() under experiments/results/{data,plots}/<tag>_sweeps

rng(1,'twister');

% ---------- Minimal cfg ----------
cfg = struct();
cfg.io = struct('save_tag', 'kMSD_sample_size_sweep');

% System & shared options
cfg.shared = struct();
cfg.shared.dyn = "k-Mass-SD";
cfg.shared.type = "standard";      % same uncertainty preset across D
cfg.shared.p_extr = 0.3;

% Identification + reachability options (CORA-style)
cfg.shared.options_reach = struct( ...
    'zonotopeOrder',      100, ...
    'tensorOrder',        2, ...
    'errorOrder',         1, ...
    'tensorOrderOutput',  2, ...
    'verbose',            false);

% Conformance (Gray) base options
cfg.shared.cs_base = struct( ...
    'robustnessMargin', 1e-9, ...
    'verbose', false, ...
    'cost', "interval", ...
    'constraints', "half");

% Data budgets (fixed)
cfg.shared.n_m = 10;    % input traj count
cfg.shared.n_s = 20;   % samples per traj
cfg.shared.n_k = 4;    % horizon (train)
cfg.shared.n_m_val = 2;    % val traj count

cfg.shared.n_s_val = cfg.shared.n_s;
cfg.shared.n_k_val = cfg.shared.n_k;

% DDRA noise setup (fixed)
cfg.ddra = struct();
cfg.ddra.eta_w   = 1;      % number of W generators
cfg.ddra.alpha_w = 0.00;   % W scale

% Gray methods (keep simple/fast)
cfg.gray = struct();
cfg.gray.methodsGray = ["grayLS"];

rcsi_lbl = rcsi_label_from_cfg(cfg);                
cfg.io.save_tag = sprintf('%s_%s', cfg.io.save_tag, rcsi_lbl);  

% --- shared noise policy to keep DDRA/Gray comparable when not studying noise
cfg.shared.noise_for_gray = false;   % if false => Gray/RCSI runs with W = 0
cfg.shared.noise_for_ddra = true;    % DDRA uses W unless this is set to false


% ---------- Sweep grid ----------
sweep_grid = struct();
sweep_grid.D_list       = [2];
sweep_grid.alpha_w_list = cfg.ddra.alpha_w;  % keep W fixed
sweep_grid.n_m_list = [2 4 8 16 32 64 128];
sweep_grid.n_m_list = [2 16 32 64];
sweep_grid.n_s_list = 5;
sweep_grid.n_k_list = 10;
sweep_grid.pe_list = {struct('mode','randn','order',2,'deterministic',true,'strength',1)}; % keep excitation mode fixed

% New: Memory efficiency toggles
cfg.lowmem = struct();
cfg.lowmem.gray_check_contain = true;   % don’t do expensive Gray containment
cfg.lowmem.store_ddra_sets    = true;   % don’t keep DDRA sets; compute metrics on the fly
cfg.lowmem.append_csv         = true;    % stream CSV row-by-row; don’t keep a giant table
cfg.lowmem.zonotopeOrder_cap  = 50;      % optional: lower order to shrink sets in memory

% ---------- Run ----------
SUMMARY = run_sweeps(cfg, sweep_grid);
SUMMARY = ensure_time_totals(SUMMARY);

% ---------- Quick plots: 4 panels (total / learn / validation / inference) ----------
x = coerce_numeric(SUMMARY.n_m);
colors = struct('ddra',[0.23 0.49 0.77],'gray',[0.85 0.33 0.10]);

f = figure('Name','Runtime panels','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% 1) TOTAL
nexttile; hold on;
plot(x, SUMMARY.t_ddra_total, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName','DDRA');
plot(x, SUMMARY.t_gray_total, '-s', 'Color',colors.gray, 'LineWidth',1.6, ...
     'DisplayName', ['RCSI-' rcsi_lbl]);
xlabel('n_m (number of input trajectories)'); ylabel('Seconds');
title(['Total runtime vs n_m  (RCSI: ' rcsi_lbl ')']); grid on; legend('Location','best');

% 2) LEARNING
nexttile; hold on;
plot(x, SUMMARY.t_ddra_learn, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName','DDRA learn');
plot(x, SUMMARY.t_gray_learn, '-s', 'Color',colors.gray, 'LineWidth',1.6, ...
     'DisplayName', ['RCSI-' rcsi_lbl ' learn']);
xlabel('n_m'); ylabel('Seconds');
title(['Learning runtime  (RCSI: ' rcsi_lbl ')']); grid on; legend('Location','best');

% 3) VALIDATION
nexttile; hold on;
plot(x, SUMMARY.t_ddra_check, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName','DDRA check');
plot(x, SUMMARY.t_gray_val,   '-s', 'Color',colors.gray, 'LineWidth',1.6, ...
     'DisplayName', ['RCSI-' rcsi_lbl ' validate']);
xlabel('n_m'); ylabel('Seconds');
title(['Validation / Check runtime  (RCSI: ' rcsi_lbl ')']); grid on; legend('Location','best');

% 4) INFERENCE
nexttile; hold on;
plot(x, SUMMARY.t_ddra_infer, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName','DDRA infer');
plot(x, SUMMARY.t_gray_infer, '-s', 'Color',colors.gray, 'LineWidth',1.6, ...
     'DisplayName', ['RCSI-' rcsi_lbl ' infer']);
xlabel('n_m'); ylabel('Seconds');
title(['Inference runtime  (RCSI: ' rcsi_lbl ')']); grid on; legend('Location','best');

% Save with method in filename
[plots_dir, ~] = init_io(cfg);
save_plot(f, plots_dir, ['runtime_panels_vs_nm_' rcsi_lbl], 'Formats', {'png','pdf'}, 'Resolution', 200);

% -------- Fidelity / Conservatism panels (vs n_m) --------
x_nm        = coerce_numeric(SUMMARY.n_m);
cval_ddra   = coerce_numeric(SUMMARY.cval_ddra);
cval_gray   = coerce_numeric(SUMMARY.cval_gray);
sizeI_ddra  = coerce_numeric(SUMMARY.sizeI_ddra);
sizeI_gray  = coerce_numeric(SUMMARY.sizeI_gray);

f2 = figure('Name','Fidelity & Conservatism','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% (1) Fidelity
nexttile; hold on;
plot(x_nm, cval_ddra, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName','DDRA');
plot(x_nm, cval_gray, '-s', 'Color',colors.gray, 'LineWidth',1.6, ...
     'DisplayName', ['RCSI-' rcsi_lbl]);
xlabel('n_m (number of input trajectories)'); ylabel('Containment on validation (%)');
title(['Fidelity vs n_m  (RCSI: ' rcsi_lbl ')']); grid on; legend('Location','best');

% (2) Conservatism
nexttile; hold on;
plot(x_nm, sizeI_ddra, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName','DDRA');
plot(x_nm, sizeI_gray, '-s', 'Color',colors.gray, 'LineWidth',1.6, ...
     'DisplayName', ['RCSI-' rcsi_lbl]);
xlabel('n_m (number of input trajectories)'); ylabel('Aggregated interval size (proxy)');
title(['Conservatism vs n_m  (RCSI: ' rcsi_lbl ')']); grid on; legend('Location','best');

% Save with method in filename
[plots_dir, ~] = init_io(cfg);
save_plot(f2, plots_dir, ['fidelity_conservatism_vs_nm_' rcsi_lbl], 'Formats', {'png','pdf'}, 'Resolution', 200);

close all force