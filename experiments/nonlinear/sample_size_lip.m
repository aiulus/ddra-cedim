%% HOW TO USE (Lipschitz / Nonlinear Sample-Size Sweep)
% What it does:
%   Uses the scalable Lipschitz MSD system ("kLipMSD") and sweeps n_m
%   (number of distinct input trajectories). Produces fidelity/conservatism
%   plots for RCSI black-box identification (GP/CGP). Runtime is recorded too.
%
% Key knobs:
%   - sweep_grid.n_m_list (e.g., [2 4 8 16 32])
%   - Fix dimension via sweep_grid.D_list = <D>;
%   - Keep excitation fixed: sweep_grid.pe_list = {struct('mode','randn')};
%
% Memory-efficiency toggles:
%   cfg.lowmem.append_csv = true;       % stream rows to CSV
%
% Outputs:
%   - CSV + plots via init_io() under experiments/results/{data,plots}/<tag>_sweeps

rng(1,'twister');

% ---------- Minimal cfg ----------
cfg = struct();
cfg.io = struct('save_tag', 'kLipMSD_sample_size_sweep');

% System & shared options
cfg.shared = struct();
cfg.shared.dyn   = "kLipMSD";      % <-- nonlinear Lipschitz system
cfg.shared.type  = "standard";     % uncertainty preset 
cfg.shared.p_extr = 0.3;           % test-suite extreme input prob.

% CORA reachability options (used by validateReach / reach)
cfg.shared.options_reach = struct( ...
    'zonotopeOrder',      100, ...
    'tensorOrder',        2, ...
    'errorOrder',         1, ...
    'tensorOrderOutput',  2, ...
    'verbose',            false);

% Conformance (RCSI) base options
cfg.shared.cs_base = struct( ...
    'robustnessMargin', 1e-9, ...
    'verbose', false, ...
    'cost', "interval", ...
    'constraints', "half");

% Data budgets
cfg.shared.n_m     = 3;   % identification: # input trajectories
cfg.shared.n_s     = 10;  % samples per input trajectory
cfg.shared.n_k     = 6;   % horizon (train)
cfg.shared.n_m_val = cfg.shared.n_m;   % validation: # input trajectories
cfg.shared.n_s_val = cfg.shared.n_s;
cfg.shared.n_k_val = cfg.shared.n_k;

% Black-box RCSI methods 
cfg.black = struct();
cfg.black.methodsBlack = ["blackCGP"];   % or ["blackGP","blackCGP"]
%cfg.black.approx = tiny_preset(); % Comment out if full GPTIPS run wanted

cfg.black.approx = struct( ... % Also Comment out if full GPTIPS run wanted
    'p', 1, ...                       
    'gp_parallel', false, ...
    'gp_pop_size', 40, ...
    'gp_num_gen', 6, ...
    'gp_func_names', {{'times','plus','tanh'}}, ... % or {'times','plus','square'}
    'gp_max_genes', 2, ...
    'gp_max_depth', 3, ...
    'cgp_num_gen', 2, ...
    'cgp_pop_size_base', 20, ...
    'save_res', false, ...
    'verbose', true);

% (Optional) lighter CGP for speed
%cfg.black.approx = struct('cgp_num_gen', 5, 'cgp_pop_size_base', 10);

% Keep explicit process noise out for side-by-side comparison
cfg.shared.noise_for_gray = false;

% --- method label for filenames
rcsi_lbl = "blackGP";                 % if multiple, plots use the first
if numel(cfg.black.methodsBlack) > 0
    rcsi_lbl = char(cfg.black.methodsBlack(1));
end
cfg.io.save_tag = sprintf('%s_%s', cfg.io.save_tag, rcsi_lbl);

% ---------- Sweep grid ----------
sweep_grid = struct();
%sweep_grid.D_list       = 4;           % dimension of the MSD chain (q,v -> 2D states)
%sweep_grid.n_m_list     = [2 4 8 16];  % <-- sweep here

sweep_grid.D_list       = 2;           % dimension of the MSD chain (q,v -> 2D states)
%sweep_grid.n_m_list     = [2 4];  % <-- sweep here

%sweep_grid.n_s_list     = cfg.shared.n_s;
%sweep_grid.n_k_list     = cfg.shared.n_k;
%sweep_grid.pe_list      = {struct('mode','randn','order',2,'strength',1,'deterministic',true)};

% ---------- Low-memory toggles ----------
cfg.lowmem = struct();
cfg.lowmem.append_csv = true;          % stream CSV rows

% ---------- Run ----------
SUMMARY = run_sweeps_lip_rcsi(cfg, sweep_grid);

% ---------- Plots (Fidelity + Conservatism vs n_m) ----------
x_nm        = coerce_numeric(SUMMARY.n_m);
cval_bb     = coerce_numeric(SUMMARY.cval_black);
sizeI_bb    = coerce_numeric(SUMMARY.sizeI_black);

colors = struct('bb',[0.2 0.55 0.3]);

f = figure('Name','Lipschitz | Fidelity & Conservatism vs n_m','Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% (1) Fidelity (containment on validation)
nexttile; hold on; grid on;
plot(x_nm, cval_bb, '-o', 'Color',colors.bb, 'LineWidth',1.6, ...
     'DisplayName', ['RCSI-' rcsi_lbl]);
xlabel('n_m (input trajectories)'); ylabel('Containment on validation (%)');
title(['Fidelity vs n_m  (RCSI: ' rcsi_lbl ', kLipMSD)']); legend('Location','best');

% (2) Conservatism proxy (aggregated interval size)
nexttile; hold on; grid on;
plot(x_nm, sizeI_bb, '-s', 'Color',colors.bb, 'LineWidth',1.6, ...
     'DisplayName', ['RCSI-' rcsi_lbl]);
xlabel('n_m (input trajectories)'); ylabel('Aggregated interval size (proxy)');
title(['Conservatism vs n_m  (RCSI: ' rcsi_lbl ', kLipMSD)']); legend('Location','best');

% Save
[plots_dir, ~] = init_io(cfg);
save_plot(f, plots_dir, ['lip_fid_cons_vs_nm_' rcsi_lbl], 'Formats', {'png','pdf'}, 'Resolution', 200);

close all force

%% >>> Helpers that define smaller GPTIPS runs for faster debugging
function approx = tiny_preset()
    approx = struct( ...
      'p', 1, 'gp_num_gen', 5, 'cgp_num_gen', 2, 'cgp_pop_size_base', 20, ...
      'save_res', false, 'verbose', true, ...
      'gp', struct('pop_size', 40,'num_gen', 6,'num_runs', 1, ...
                   'tournament_size', 5,'elite_fraction', 0.1, ...
                   'max_depth', 3,'max_nodes', 60,'max_genes', 2, ...
                   'const_range',[-1 1],'functions',{{'PLUS','MINUS','TIMES','TANH'}}, ...
                   'use_rdivide', false, 'complexity','expressional', ...
                   'lexicographic', true,'parallel', false, ...
                   'fitness_goal', 5e-3,'stall_limit', 3,'seed', 1));
end

function approx = fast_preset()
    approx = struct( ...
      'p', 1, 'gp_num_gen', 12, 'cgp_num_gen', 4, 'cgp_pop_size_base', 40, ...
      'save_res', false, 'verbose', true, ...
      'gp', struct('pop_size', 80,'num_gen', 12,'num_runs', 1, ...
                   'tournament_size', 7,'elite_fraction', 0.15, ...
                   'max_depth', 4,'max_nodes', 120,'max_genes', 3, ...
                   'const_range',[-2 2],'functions',{{'PLUS','MINUS','TIMES','TANH','SQUARE'}}, ...
                   'use_rdivide', false, 'complexity','expressional', ...
                   'lexicographic', true,'parallel', false, ...
                   'fitness_goal', 1e-3,'stall_limit', 5,'seed', 1));
end
