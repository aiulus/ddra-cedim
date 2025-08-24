rng(1,'twister');

% ---------- Minimal cfg ----------
cfg = struct();
cfg.io = struct('save_tag', 'kMSD_dim_sweep');

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
cfg.shared.n_m = 3;    % input traj count
cfg.shared.n_s = 20;   % samples per traj
cfg.shared.n_k = 4;    % horizon (train)
cfg.shared.n_m_val = 2;    % val traj count

cfg.shared.n_s_val = cfg.shared.n_s;
cfg.shared.n_k_val = cfg.shared.n_k;

% DDRA noise setup (fixed)
cfg.ddra = struct();
cfg.ddra.eta_w   = 1;      % number of W generators
cfg.ddra.alpha_w = 0.01;   % W scale

% Gray methods (keep simple/fast)
cfg.gray = struct();
cfg.gray.methodsGray = ["graySeq"];

% ---------- Sweep grid ----------
sweep_grid = struct();
sweep_grid.D_list       = [2];
sweep_grid.alpha_w_list = cfg.ddra.alpha_w;  % keep W fixed
sweep_grid.n_m_list = [2 4 8 16 32];
sweep_grid.n_m_list = [2 4];
sweep_grid.n_s_list = 1;
sweep_grid.n_k_list = 10;
sweep_grid.pe_list = {struct('mode','randn')};  % keep excitation mode fixed

% ---------- Run ----------
SUMMARY = run_sweeps(cfg, sweep_grid);

% ---------- Quick plots: 4 panels (total / learn / validation / inference) ----------
% ensure totals exist
if ~ismember('t_ddra_total', SUMMARY.Properties.VariableNames)
    SUMMARY.t_ddra_total = SUMMARY.t_ddra_learn + SUMMARY.t_ddra_check + SUMMARY.t_ddra_infer;
end
if ~ismember('t_gray_total', SUMMARY.Properties.VariableNames)
    SUMMARY.t_gray_total = SUMMARY.t_gray_learn + SUMMARY.t_gray_val   + SUMMARY.t_gray_infer;
end

x = SUMMARY.n_m;                       
colors = struct('ddra',[0.23 0.49 0.77],'gray',[0.85 0.33 0.10]);

f = figure('Name','Runtime panels','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% 1) TOTAL
nexttile; hold on;
plot(x, SUMMARY.t_ddra_total, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName','DDRA');
plot(x, SUMMARY.t_gray_total, '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName','Gray');
xlabel('n_m (number of input trajectories)'); ylabel('Seconds');
title('Total runtime vs n_m'); grid on; legend('Location','best');

% 2) LEARNING
nexttile; hold on;
plot(x, SUMMARY.t_ddra_learn, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName','DDRA learn');
plot(x, SUMMARY.t_gray_learn, '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName','Gray learn');
xlabel('n_m'); ylabel('Seconds');
title('Learning runtime'); grid on; legend('Location','best');

% 3) VALIDATION  (DDRA “check” vs Gray validation)
nexttile; hold on;
plot(x, SUMMARY.t_ddra_check, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName','DDRA check');
plot(x, SUMMARY.t_gray_val,   '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName','Gray validate');
xlabel('n_m'); ylabel('Seconds');
title('Validation / Check runtime'); grid on; legend('Location','best');

% 4) INFERENCE
nexttile; hold on;
plot(x, SUMMARY.t_ddra_infer, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName','DDRA infer');
plot(x, SUMMARY.t_gray_infer, '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName','Gray infer');
xlabel('n_m'); ylabel('Seconds');
title('Inference runtime'); grid on; legend('Location','best');

% save into experiments/results/plots/
[plots_dir, ~] = init_io(cfg);  
out_png = fullfile(plots_dir, 'runtime_panels_vs_nm.png');
exportgraphics(f, out_png, 'Resolution', 200);

disp(['Saved runtime panels -> ' out_png]);
