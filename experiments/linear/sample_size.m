%% HOW TO USE — Sample-Size Sweep (DDRA vs RCSI/Gray)
% What this does
%   Sweeps the number of distinct nominal input trajectories n_m (optionally n_s, n_k)
%   at fixed dimension, producing:
%     (i) runtime panels (total / learn / validation / inference)
%     (ii) fidelity (containment %) and conservatism (interval-size proxy) plots.
%
% Equal-setting evaluation protocol 
%   • Shared datasets: We generate TRAIN/VAL datasets once via ddra_generate_data and
%     pass the exact same (x0,u,y) sequences to both DDRA and Gray.
%   • Unified noise policy: A single switch controls whether both methods use W≠0 or W=0.
%       - Effective rule: use_noise = shared.noise_for_gray && shared.noise_for_ddra
%       - If false → Gray runs with W=0 and DDRA uses W=0 (hard zero zonotope).
%   • Ridge guard (DDRA): 
%       cfg.ddra.allow_ridge=false (default) → rank-deficient Z → the sweep point is skipped
%       (row marked as skipped in CSV). If true → ridge is used and uncertainty is widened.
%   • Metrics are identical across methods on the SAME points:
%       Containment: point-in-interval-hull of the OUTPUT reachable set (contains_interval).
%       Size proxy: aggregated output-interval width across VAL (not state size).
%   • Set reduction policy: consistent Girard reduction with cap
%       (cfg.lowmem.zonotopeOrder_cap applies in streaming; options_reach.zonotopeOrder otherwise).
%
%% Key knobs
%   • Dimension:  sweep_grid.D_list = [D]; e.g., 2 or 5
%   • Sample size axes: sweep_grid.n_m_list (e.g., [2 4 8 16 32 64])
%                       optionally n_s_list, n_k_list
%   • Excitation (fixed recommended): sweep_grid.pe_list = {struct('mode','randn',...)}
%
% Memory / IO toggles (safe defaults shown)
%   cfg.lowmem.gray_check_contain = true;    % light, interval-based containment on VAL
%   cfg.lowmem.store_ddra_sets    = true;    % store sets; set false to use streaming path
%   cfg.lowmem.append_csv         = true;    % stream rows to CSV on the go
%   cfg.lowmem.zonotopeOrder_cap  = 50;      % cap order for memory friendliness
%
%% Outputs
%   CSV:   experiments/results/data/<save_tag>_sweeps/summary.csv
%   Plots: experiments/results/plots/<save_tag>_sweeps/*.png|pdf


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
cfg.shared.n_s = 2;   % samples per traj
cfg.shared.n_k = 4;    % horizon (train)
cfg.shared.n_m_val = 2;    % val traj count

cfg.shared.n_s_val = cfg.shared.n_s;
cfg.shared.n_k_val = cfg.shared.n_k;

% DDRA noise setup (fixed)
cfg.ddra = struct();
cfg.ddra.eta_w   = 1;      % number of W generators
cfg.ddra.alpha_w = 0.00;   % W scale
% --- DDRA ridge guard (defaults) ---
cfg.ddra.allow_ridge   = false;   % if false and rank-deficient -> skip point
cfg.ddra.lambda        = 1e-8;    % ridge lambda when allowed
cfg.ddra.ridge_gamma   = 1.0;     % scale for added uncertainty
cfg.ddra.ridge_policy  = "MAB";   % "MAB" (add generator to M_AB) or "W" (inflate W)

% Gray methods (keep simple/fast)
cfg.gray = struct();
cfg.gray.methodsGray = ["grayLS"];

rcsi_lbl = rcsi_label_from_cfg(cfg);                
cfg.io.save_tag = sprintf('%s_%s', cfg.io.save_tag, rcsi_lbl);  

% --- shared noise policy to keep DDRA/Gray comparable when not studying noise
cfg.shared.noise_for_gray = false;   % if false => Gray/RCSI runs with W = 0
cfg.shared.noise_for_ddra = true;    % DDRA uses W unless this is set to false
cfg.shared.use_noise = false;   % or true


% ---------- Sweep grid ----------
sweep_grid = struct();
sweep_grid.D_list       = [5];
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

%% ---------- Run ----------
SUMMARY = run_sweeps(cfg, sweep_grid);
SUMMARY = ensure_time_totals(SUMMARY);

disp([SUMMARY.cval_gray SUMMARY.cval_ddra])
disp([SUMMARY.sizeI_gray SUMMARY.sizeI_ddra])

%% ---------- Visualization: 4 panels (total / learn / validation / inference) ----------
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