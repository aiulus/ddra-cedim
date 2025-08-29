%% HOW TO USE — State-Dimension Scalability (DDRA vs RCSI/Gray)
% What this does
%   Sweeps the state dimension D on the k-Mass-Spring-Damper chain and compares:
%     (i) fidelity (containment %) (ii) conservatism (interval-size proxy) (iii) runtime
%
% Equal-setting evaluation protocol 
%   • Shared datasets: DDRA and Gray use the same (x0,u,y) sequences for TRAIN/VAL.
%   • Unified noise policy: use_noise = shared.noise_for_gray && shared.noise_for_ddra.
%   • Ridge guard: if cfg.ddra.allow_ridge=false and Z is rank-deficient → skip & mark row.
%   • Metrics and reduction: identical containment/size metrics and common Girard cap.
%
%% Key knobs
%   • Dimension sweep: sweep_grid.D_list = [2 3 4 5 ...]
%   • Fixed budgets: cfg.shared.n_m, n_s, n_k (+ *_val)
%   • Excitation: sweep_grid.pe_list = {struct('mode','randn')} (simple default)
%
% Memory / IO toggles
%   cfg.lowmem.gray_check_contain = false;   % lightweight interval-based checks
%   cfg.lowmem.store_ddra_sets    = false;   % streaming DDRA metrics
%   cfg.lowmem.append_csv         = true;    % stream rows to disk
%   cfg.lowmem.zonotopeOrder_cap  = 50;      % reduce memory footprint
%
%% Outputs
%   CSV:   experiments/results/data/<save_tag>_sweeps/summary.csv
%   Plots: experiments/results/plots/<save_tag>_sweeps/*.png|pdf


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

% ---------- Sweep grid (vary only D) ----------
sweep_grid = struct();
%sweep_grid.D_list = [2 3 4 5 6];       % dimensional scalability
%sweep_grid.D_list = [2 3 4 5 6 7 8 9 10];
sweep_grid.D_list = 2:2:10;
sweep_grid.alpha_w_list = cfg.ddra.alpha_w;  % keep W fixed
sweep_grid.n_m_list = cfg.shared.n_m;
sweep_grid.n_s_list = cfg.shared.n_s;
sweep_grid.n_k_list = cfg.shared.n_k;
PE_orders = [4];
sweep_grid.pe_list = [ ...
    arrayfun(@(L) struct('mode','sinWave','order',L,'strength',1,'deterministic',true), PE_orders, 'uni',0) ...
];


% New: Memory efficiency toggles
cfg.lowmem = struct();
cfg.lowmem.gray_check_contain = false;   % don’t do expensive Gray containment
cfg.lowmem.store_ddra_sets    = false;   % don’t keep DDRA sets; compute metrics on the fly
cfg.lowmem.append_csv         = true;    % stream CSV row-by-row; don’t keep a giant table
cfg.lowmem.zonotopeOrder_cap  = 50;      % optional: lower order to shrink sets in memory

% --- shared noise policy to keep DDRA/Gray comparable when not studying noise
cfg.shared.noise_for_gray = false;   % if false => Gray/RCSI runs with W = 0
cfg.shared.noise_for_ddra = true;    % DDRA uses W unless this is set to false
cfg.shared.use_noise = false;   % or true

% --- plotting mode + knobs
cfg.io.plot_mode     = "offline";   % "online" | "offline" | "both"
cfg.io.make_reach_plot = true;      % online: produce figures during run
cfg.io.plot_rows       = [1];       % which sweep rows to plot 
cfg.io.plot_dims       = [1 2];     % output dims
cfg.io.plot_every_k    = 1;         % plot every k-th step (declutter)
cfg.io.save_artifacts  = true;      % offline: keep .mat files for peeking/plotting

cfg.io.base_dir = fileparts(fileparts(mfilename('fullpath'))); % or hard-code
cfg.allow_parallel = false;  % keep serial


%% ---------- Run ----------
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

% ---------- Quick plots ----------
% Note: results saved to experiments/results/data/<tag>_sweeps/summary.csv
[fplots_dir, ~] = init_io(cfg);   % where we'll save the figures

% --- Figure 1: Fidelity / Conservatism vs D ---
f_fid = figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Containment vs D (Fidelity)
nexttile;
plot(SUMMARY.D, SUMMARY.cval_ddra, '-o', 'DisplayName', 'DDRA'); hold on;
plot(SUMMARY.D, SUMMARY.cval_gray, '-s', 'DisplayName', 'GraySeq');
xlabel('Dimension D'); ylabel('Containment on validation (%)');
title('Fidelity vs Dimension'); grid on; legend('Location','best');

% Size proxy vs D (Conservatism)
nexttile;
plot(SUMMARY.D, SUMMARY.sizeI_ddra, '-o', 'DisplayName', 'DDRA'); hold on;
plot(SUMMARY.D, SUMMARY.sizeI_gray, '-s', 'DisplayName', 'GraySeq');
xlabel('Dimension D'); ylabel('Aggregated interval size (proxy)');
title('Conservatism proxy vs Dimension'); grid on; legend('Location','best');

% Save Figure 1
fname1 = sprintf('%s_fidcons_vs_D', cfg.io.save_tag);
if exist('save_plot','file')
    save_plot(f_fid, fplots_dir, fname1, 'Formats', {'png','pdf'}, 'Resolution', 200);
else
    exportgraphics(f_fid, fullfile(fplots_dir, [fname1 '.png']), 'Resolution', 200);
    exportgraphics(f_fid, fullfile(fplots_dir, [fname1 '.pdf']), 'ContentType','vector');
end

disp(['Done. CSV is in ', fullfile(fplots_dir, '..', 'data', [cfg.io.save_tag '_sweeps']), filesep, 'summary.csv']);

%% -------- Runtime plots (vs D) --------
% Compute totals if not already present
if ~ismember('t_ddra_total', SUMMARY.Properties.VariableNames)
    SUMMARY.t_ddra_total = SUMMARY.t_ddra_learn + SUMMARY.t_ddra_check + SUMMARY.t_ddra_infer;
end
if ~ismember('t_gray_total', SUMMARY.Properties.VariableNames)
    SUMMARY.t_gray_total = SUMMARY.t_gray_learn + SUMMARY.t_gray_val   + SUMMARY.t_gray_infer;
end

% One consolidated figure: total + 3 phase panels
f_run = figure('Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% (1) Total runtime
nexttile;
plot(SUMMARY.D, SUMMARY.t_ddra_total, '-o', 'LineWidth',1.5, 'DisplayName','DDRA total'); hold on;
plot(SUMMARY.D, SUMMARY.t_gray_total, '-s', 'LineWidth',1.5, 'DisplayName','Gray total');
xlabel('Dimension D'); ylabel('Seconds'); title('Total runtime'); grid on; legend('Location','best');

% (2) Learning phase
nexttile;
plot(SUMMARY.D, SUMMARY.t_ddra_learn, '-o', 'LineWidth',1.5, 'DisplayName','DDRA learn'); hold on;
plot(SUMMARY.D, SUMMARY.t_gray_learn, '-s', 'LineWidth',1.5, 'DisplayName','Gray learn');
xlabel('Dimension D'); ylabel('Seconds'); title('Learning phase'); grid on; legend('Location','best');

% (3) Validation phase
nexttile;
plot(SUMMARY.D, SUMMARY.t_ddra_check, '-o', 'LineWidth',1.5, 'DisplayName','DDRA check'); hold on;
plot(SUMMARY.D, SUMMARY.t_gray_val,  '-s', 'LineWidth',1.5, 'DisplayName','Gray validate');
xlabel('Dimension D'); ylabel('Seconds'); title('Validation phase'); grid on; legend('Location','best');

% (4) Inference phase
nexttile;
plot(SUMMARY.D, SUMMARY.t_ddra_infer, '-o', 'LineWidth',1.5, 'DisplayName','DDRA infer'); hold on;
plot(SUMMARY.D, SUMMARY.t_gray_infer, '-s', 'LineWidth',1.5, 'DisplayName','Gray infer');
xlabel('Dimension D'); ylabel('Seconds'); title('Inference phase'); grid on; legend('Location','best');

% Save Figure 2
fname2 = sprintf('%s_runtime_vs_D', cfg.io.save_tag);
if exist('save_plot','file')
    save_plot(f_run, fplots_dir, fname2, 'Formats', {'png','pdf'}, 'Resolution', 200);
else
    exportgraphics(f_run, fullfile(fplots_dir, [fname2 '.png']), 'Resolution', 200);
    exportgraphics(f_run, fullfile(fplots_dir, [fname2 '.pdf']), 'ContentType','vector');
end

close all force
