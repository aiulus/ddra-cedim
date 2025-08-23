rng(1,'twister');

% ---------- Minimal cfg ----------
cfg = struct();
cfg.io = struct('save_tag', 'kMSD_dim_sweep');

% System & shared options
cfg.shared = struct();
cfg.shared.dyn   = "k-Mass-SD";
cfg.shared.type  = "standard";      % same uncertainty preset across D
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
cfg.shared.n_m     = 3;    % input traj count
cfg.shared.n_s     = 20;   % samples per traj
cfg.shared.n_k     = 4;    % horizon (train)
cfg.shared.n_m_val = 2;    % val traj count

%cfg.shared.n_m     = 2;    % input traj count
%cfg.shared.n_s     = 2;   % samples per traj
%cfg.shared.n_k     = 2;    % horizon (train)
%cfg.shared.n_m_val = 2;    % val traj count

cfg.shared.n_s_val = cfg.shared.n_s;
cfg.shared.n_k_val = cfg.shared.n_k;

% DDRA noise setup (fixed)
cfg.ddra = struct();
cfg.ddra.eta_w   = 1;      % number of W generators
cfg.ddra.alpha_w = 0.01;   % W scale

% Gray methods (keep simple/fast)
cfg.gray = struct();
cfg.gray.methodsGray = ["graySeq"];

% ---------- Sweep grid (vary only D) ----------
sweep_grid = struct();
%sweep_grid.D_list       = [2 3 4 5 6];       % dimensional scalability
%sweep_grid.D_list       = [2 3 4 5 6 7 8 9 10];
sweep_grid.D_list       = [2 3 4 5];
sweep_grid.alpha_w_list = cfg.ddra.alpha_w;  % keep W fixed
sweep_grid.n_m_list     = cfg.shared.n_m;
sweep_grid.n_s_list     = cfg.shared.n_s;
sweep_grid.n_k_list     = cfg.shared.n_k;
sweep_grid.pe_list      = {struct('mode','randn')};  % keep excitation mode fixed

% ---------- Run ----------
SUMMARY = run_sweeps(cfg, sweep_grid);

% ---------- Quick plots ----------
% Note: results saved to results/<save_tag>_sweeps/summary.csv
figure; 
tiledlayout(1,2);

% Containment vs D
nexttile;
plot(SUMMARY.D, SUMMARY.cval_ddra, '-o', 'DisplayName', 'DDRA'); hold on;
plot(SUMMARY.D, SUMMARY.cval_gray, '-s', 'DisplayName', 'GraySeq');
xlabel('Dimension D'); ylabel('Containment on validation (%)');
title('Fidelity vs Dimension'); grid on; legend('Location','best');

% Size proxy vs D
nexttile;
plot(SUMMARY.D, SUMMARY.sizeI_ddra, '-o', 'DisplayName', 'DDRA'); hold on;
plot(SUMMARY.D, SUMMARY.sizeI_gray, '-s', 'DisplayName', 'GraySeq');
xlabel('Dimension D'); ylabel('Aggregated interval size (proxy)');
title('Conservatism proxy vs Dimension'); grid on; legend('Location','best');

disp('Done. CSV is in results/kMSD_dim_sweep_sweeps/summary.csv');
