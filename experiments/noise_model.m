%% Additive noise profile sweep: Fidelity & Conservatism
rng(1,'twister');

cfg = struct();
cfg.io = struct('save_tag','kMSD_noise_add');

% Shared
cfg.shared = struct();
cfg.shared.dyn   = "k-Mass-SD";
cfg.shared.type  = "standard";
cfg.shared.p_extr = 0.3;
cfg.shared.options_reach = struct('zonotopeOrder',100,'tensorOrder',2,'errorOrder',1,'tensorOrderOutput',2,'verbose',false);
cfg.shared.cs_base = struct('robustnessMargin',1e-9,'verbose',false,'cost',"interval",'constraints',"half");

% Data budgets (modest to keep it light; adjust as needed)
cfg.shared.n_m = 3;
cfg.shared.n_s = 5;
cfg.shared.n_k = 6;
cfg.shared.n_m_val = 2;
cfg.shared.n_s_val = cfg.shared.n_s;
cfg.shared.n_k_val = cfg.shared.n_k;

% DDRA noise (eta fixed; we will sweep alpha_w)
cfg.ddra = struct('eta_w',1,'alpha_w',0.01);

% Gray
cfg.gray = struct('methodsGray', ["graySeq"]);

% Memory-efficiency toggles (optional)
cfg.lowmem = struct();
cfg.lowmem.gray_check_contain = true;   % turn off expensive Gray containment check speedup
cfg.lowmem.store_ddra_sets    = false;   % stream DDRA inference size & fidelity
cfg.lowmem.append_csv         = true;    % stream rows to CSV
cfg.lowmem.zonotopeOrder_cap  = 50;      % optionally cap zonotope order during streaming

% Sweep grid: vary additive noise scale
sweep_grid = struct();
sweep_grid.D_list        = 2;            % fix state dimension for this study
sweep_grid.alpha_w_list  = [0, 0.005, 0.01, 0.02, 0.05];
sweep_grid.alpha_w_list  = [0, 0.005];
sweep_grid.n_m_list      = cfg.shared.n_m;
sweep_grid.n_s_list      = cfg.shared.n_s;
sweep_grid.n_k_list      = cfg.shared.n_k;
sweep_grid.pe_list       = {struct('mode','randn')}; % keep excitation fixed

% Run
SUMMARY = run_sweeps(cfg, sweep_grid);

% ---- Plots: Fidelity / Conservatism vs noise scale ----
[plots_dir, ~] = init_io(cfg);

f = figure('Color','w'); tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% (1) Fidelity vs alpha_w
nexttile; hold on; grid on; title('Fidelity vs noise scale \alpha_W');
plot(SUMMARY.alpha_w, SUMMARY.cval_ddra, '-o', 'DisplayName','DDRA');
plot(SUMMARY.alpha_w, SUMMARY.cval_gray, '-s', 'DisplayName','GraySeq');
xlabel('\alpha_W'); ylabel('Containment on validation (%)'); legend('Location','best');

% (2) Conservatism proxy vs alpha_w
nexttile; hold on; grid on; title('Conservatism proxy vs \alpha_W');
plot(SUMMARY.alpha_w, SUMMARY.sizeI_ddra, '-o', 'DisplayName','DDRA');
plot(SUMMARY.alpha_w, SUMMARY.sizeI_gray, '-s', 'DisplayName','GraySeq');
xlabel('\alpha_W'); ylabel('Aggregated interval size');

% Save
tag = cfg.io.save_tag;
fname = sprintf('%s_fidcons_vs_alphaW', tag);
if exist('save_plot','file')
    save_plot(f, plots_dir, fname, 'Formats', {'png','pdf'}, 'Resolution', 200);
else
    exportgraphics(f, fullfile(plots_dir, [fname '.png']), 'Resolution', 200);
    exportgraphics(f, fullfile(plots_dir, [fname '.pdf']), 'ContentType','vector');
end

disp('Additive noise sweep done.');
close all force
