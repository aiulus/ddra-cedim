%% HOW TO USE (Square / Nonlinear Sample-Size Sweep)
% What it does:
%   Uses CORA's simple NARX "Square" via custom_loadDynamics(), sweeps n_m
%   (number of unique input trajectories), and produces:
%     (i) runtime panels  and  (ii) fidelity/conservatism plots.
%
% Key knobs:
%   - sweep_grid.n_m_list (e.g., 2:1:20)
%   - Keep excitation fixed: sweep_grid.pe_list = {struct('mode','randn')}
%
% Memory-efficiency toggles:
%   cfg.lowmem.gray_check_contain = false;  % only applies to GRAY (unused here)
%   cfg.lowmem.append_csv         = true;   % stream CSV row-by-row
%   cfg.lowmem.zonotopeOrder_cap  = 50;     % reduce order in DDRA reach
%
% Outputs:
%   - CSV + plots via init_io() under experiments/results/{data,plots}/<tag>_sweeps

rng(1,'twister');

% ---------- Minimal cfg ----------
cfg = struct();
cfg.io = struct('save_tag', 'Square_sample_size_sweep');

% System & shared options (nonlinear, NARX p=1)
cfg.shared = struct();
cfg.shared.dyn    = "Square";        % <- from custom_loadDynamics()
cfg.shared.type   = "rand";          % ("standard" or "rand") preset for Square
cfg.shared.p_extr = 0.3;             % extreme-input probability during data gen

% Identification + reachability options (CORA-style)
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

% Data budgets (fixed base; validation ties to n_m unless overridden)
cfg.shared.n_m     = 5;    % input traj count (ident sanity only; sweep overrides)
cfg.shared.n_s     = 10;   % samples per traj
cfg.shared.n_k     = 4;    % horizon (train/ID)
cfg.shared.n_m_val = cfg.shared.n_m;  % validation: # input trajectories
cfg.shared.n_s_val = cfg.shared.n_s;
cfg.shared.n_k_val = cfg.shared.n_k;

% Black-box RCSI settings (fast-ish)
cfg.black = struct();
cfg.black.methodsBlack = ["blackCGP"];   % or ["blackGP","blackCGP"]
cfg.black.approx = struct( ...
    'p', 1, ...                               % Square has p=1
    'gp_parallel', false, ...
    'gp_pop_size', 5, ...
    'gp_num_gen', 3, ...
    'gp_func_names', {{'times','plus','square'}}, ...
    'gp_max_genes', 2, ...
    'gp_max_depth', 2, ...
    'cgp_num_gen', 5, ...
    'cgp_pop_size_base', 5, ...
    'save_res', false ...
);

% DDRA-Lipschitz noise setup
cfg.ddra = struct();
cfg.ddra.eta_w   = 1;      % number of W generators
cfg.ddra.alpha_w = 0.00;   % W scale (0 = turn off process noise)

% Shared noise policy (keep apples-to-apples if not studying noise)
cfg.shared.noise_for_black = false;   % set W=0 for RCSI if false
cfg.shared.noise_for_ddra  = true;    % DDRA uses W unless false

% Label in file names (based on first black-box method)
rcsi_lbl = rcsi_label_from_cfg(cfg);
cfg.io.save_tag = sprintf('%s_%s', cfg.io.save_tag, rcsi_lbl);

% ---------- Sweep grid ----------
sweep_grid = struct();
sweep_grid.n_m_list = 2:1:20;         % main sweep
sweep_grid.n_m_list = 2;    
sweep_grid.n_s_list = cfg.shared.n_s; % keep fixed
sweep_grid.n_k_list = cfg.shared.n_k; % keep fixed
sweep_grid.pe_list  = {struct('mode','randn')};  % fixed excitation

% ---------- Low-memory toggles ----------
cfg.lowmem = struct();
cfg.lowmem.gray_check_contain = true;      % no gray here, but harmless
cfg.lowmem.append_csv         = true;      % stream CSV rows
cfg.lowmem.zonotopeOrder_cap  = 50;        % DDRA reduction cap

% ---------- Run ----------
SUMMARY = run_sweeps_square_black_vs_ddraLip(cfg, sweep_grid);
SUMMARY = ensure_time_totals_square(SUMMARY);  % totals for plots (local helper)

% ---------- Plots: Runtime panels (total / learn / validation / inference) ----------
x = coerce_numeric(SUMMARY.n_m);
colors = struct('ddra',[0.23 0.49 0.77],'black',[0.20 0.55 0.30]);

f = figure('Name','Square | Runtime panels','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% 1) TOTAL
nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_total,  '-o', 'Color',colors.ddra,  'LineWidth',1.6, 'DisplayName','DDRA-Lip');
plot(x, SUMMARY.t_black_total, '-s', 'Color',colors.black, 'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl]);
xlabel('n_m (input trajectories)'); ylabel('Seconds');
title(['Total runtime vs n_m  (RCSI: ' rcsi_lbl ')']); legend('Location','best');

% 2) LEARNING
nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_build,  '-o', 'Color',colors.ddra,  'LineWidth',1.6, 'DisplayName','DDRA build');
plot(x, SUMMARY.t_black_learn, '-s', 'Color',colors.black, 'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl ' learn']);
xlabel('n_m'); ylabel('Seconds');
title('Learning / Build runtime'); legend('Location','best');

% 3) VALIDATION
nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_val, '-o', 'Color',colors.ddra,  'LineWidth',1.6, 'DisplayName','DDRA validate');
plot(x, SUMMARY.t_black_val,'-s', 'Color',colors.black, 'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl ' validate']);
xlabel('n_m'); ylabel('Seconds');
title('Validation runtime'); legend('Location','best');

% 4) INFERENCE
nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_reach_avg, '-o', 'Color',colors.ddra,  'LineWidth',1.6, 'DisplayName','DDRA reach/case (avg)');
plot(x, SUMMARY.t_black_infer,    '-s', 'Color',colors.black, 'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl ' infer']);
xlabel('n_m'); ylabel('Seconds');
title('Inference runtime'); legend('Location','best');

[plots_dir, ~] = init_io(cfg);
save_plot(f, plots_dir, ['square_runtime_panels_vs_nm_' rcsi_lbl], 'Formats', {'png','pdf'}, 'Resolution', 200);

% ---------- Plots: Fidelity & Conservatism vs n_m ----------
x_nm        = coerce_numeric(SUMMARY.n_m);
cval_ddra   = coerce_numeric(SUMMARY.cval_ddra);
cval_black  = coerce_numeric(SUMMARY.cval_black);
sizeI_ddra  = coerce_numeric(SUMMARY.sizeI_ddra);
sizeI_black = coerce_numeric(SUMMARY.sizeI_black);

f2 = figure('Name','Square | Fidelity & Conservatism','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% (1) Fidelity
nexttile; hold on; grid on;
plot(x_nm, cval_ddra,  '-o', 'Color',colors.ddra,  'LineWidth',1.6, 'DisplayName','DDRA-Lip');
plot(x_nm, cval_black, '-s', 'Color',colors.black, 'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl]);
xlabel('n_m (input trajectories)'); ylabel('Containment on validation (%)');
title(['Fidelity vs n_m  (RCSI: ' rcsi_lbl ')']); legend('Location','best');

% (2) Conservatism (proxy)
nexttile; hold on; grid on;
plot(x_nm, sizeI_ddra,  '-o', 'Color',colors.ddra,  'LineWidth',1.6, 'DisplayName','DDRA-Lip');
plot(x_nm, sizeI_black, '-s', 'Color',colors.black, 'LineWidth',1.6, 'DisplayName',['RCSI-' rcsi_lbl]);
xlabel('n_m (input trajectories)'); ylabel('Aggregated interval size (proxy)');
title(['Conservatism vs n_m  (RCSI: ' rcsi_lbl ')']); legend('Location','best');

save_plot(f2, plots_dir, ['square_fidelity_conservatism_vs_nm_' rcsi_lbl], 'Formats', {'png','pdf'}, 'Resolution', 200);
close all force



% =================== Helpers ===================

function [Xsets, sizeI] = ddra_reach_lipschitz(R0, U, W, D, C)
% Alg. 6: local LS model + Z_L + optional Z_eps (from Lipschitz).
% R0: initial set (zonotope), U: single input set or cell{1..N}, W: process noise
% D: struct(Xminus,Uminus,Xplus); C.nlip: .ZepsFlag, .Lvec, .gamma
    if ~iscell(U), U = repmat({U}, 1, getfielddef(C.shared,'n_k_val', getfielddef(C.shared,'n_k',1))); end
    N  = numel(U);
    nx = size(center(R0),1);

    % LS linearization about sample means
    xS = mean(D.Xminus,2);
    uS = mean(D.Uminus,2);
    Phi = [ones(1,size(D.Xminus,2)); D.Xminus - xS; D.Uminus - uS];   % [1; x-x*; u-u*]
    Mhat = (D.Xplus) * pinv(Phi);

    % residual enclosure Z_L (component-wise box)
    res = D.Xplus - Mhat*Phi;
    rmax = max(res,[],2); rmin = min(res,[],2);
    ZL = zonotope(0.5*(rmax+rmin), diag(0.5*(rmax-rmin)));

    % Z_eps (optional): eps_i = 0.5 * L_i * gamma_i
    Zeps = zonotope(zeros(nx,1));
    if getfielddef(C.nlip,'ZepsFlag', true)
        Lvec  = getfielddef(C.nlip,'Lvec', []);
        gamma = getfielddef(C.nlip,'gamma',[]);
        if isempty(Lvec) || isempty(gamma)
            Lvec = zeros(nx,1); gamma = zeros(nx,1); % safe no-op
        end
        eps = 0.5 * abs(Lvec(:)) .* abs(gamma(:));
        Zeps = zonotope(zeros(nx,1), diag(eps));
    end

    Xsets  = cell(N+1,1); Xsets{1} = R0;
    sizeI  = 0;
    for k=1:N
        Uk = U{k};
        Xaff  = cartProd(Xsets{k} - xS, Uk - uS);
        Xnext = Mhat * aug1(Xaff) + W + ZL + Zeps;
        Xnext = reduce(Xnext, 'girard', getfielddef(C.lowmem,'zonotopeOrder_cap',100));
        Xsets{k+1} = Xnext;

        Iv = interval(Xnext); sizeI = sizeI + sum(abs(Iv.sup - Iv.inf));
    end
end


