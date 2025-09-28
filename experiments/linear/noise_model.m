%% Additive process-noise sweep: Fidelity & Conservatism (DDRA vs RCSI/Gray)
rng(1,'twister');

cfg = struct();
cfg.io = struct('save_tag','kMSD_noise_add');

% Shared
cfg.shared = struct();
cfg.shared.dyn   = "k-Mass-SD";
cfg.shared.type  = "standard";
cfg.shared.p_extr = 0.3;
cfg.shared.options_reach = struct('zonotopeOrder',100,'verbose',false);
cfg.shared.cs_base = struct('robustnessMargin',1e-9,'verbose',false,'cost',"interval",'constraints',"half");

% Data budgets (modest to keep light)
cfg.shared.n_m = 30;  cfg.shared.n_s = 5;  cfg.shared.n_k = 6;
cfg.shared.n_m_val = 2;
cfg.shared.n_s_val = cfg.shared.n_s; cfg.shared.n_k_val = cfg.shared.n_k;

% Noise policy — EXPLICIT for this experiment
% (We want both methods to use identical process noise W where applicable)
cfg.shared.noise_for_gray = true;
cfg.shared.noise_for_ddra = true;
cfg.shared.use_noise      = true;   % overall

% DDRA noise (eta fixed; sweep alpha_w)
cfg.ddra = struct('eta_w',1,'alpha_w',0.01);
cfg.ddra.allow_ridge  = true;     % ok to allow ridge here
cfg.ddra.lambda       = 1e-8;
cfg.ddra.ridge_gamma  = 1.0;
cfg.ddra.ridge_policy = "MAB";

% Gray
cfg.gray = struct('methodsGray', ["graySeq"]);

% Efficiency toggles
cfg.lowmem = struct();
cfg.lowmem.gray_check_contain = true;
cfg.lowmem.store_ddra_sets    = false;
cfg.lowmem.append_csv         = true;
cfg.lowmem.zonotopeOrder_cap  = 50;

% Sweep grid: vary additive process noise scale α_W
sweep_grid = struct();
sweep_grid.D_list        = 2;   % fix D for this study
sweep_grid.alpha_w_list  = [0, 0.005, 0.01, 0.02, 0.05];
sweep_grid.n_m_list      = cfg.shared.n_m;
sweep_grid.n_s_list      = cfg.shared.n_s;
sweep_grid.n_k_list      = cfg.shared.n_k;
sweep_grid.pe_list       = {struct('mode','randn','deterministic',true)};

% Optional plots config
cfg.io.plot_mode       = "offline";
cfg.io.plot_dims       = [1 2];
cfg.io.save_artifacts  = false;   % set true to enable artifact-derived plots
cfg.io.base_dir        = fileparts(fileparts(mfilename('fullpath')));
cfg.allow_parallel     = false;
cfg.metrics.enhanced   = true;    % collect richer metrics if cheap enough

% Run
SUMMARY = run_sweeps(cfg, sweep_grid);
SUMMARY = ensure_time_totals(SUMMARY);

% Plots: Fidelity / Conservatism vs α_W
[plots_dir, results_dir] = init_io(cfg);

f = figure('Color','w'); tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on; grid on; title('Fidelity vs noise scale \alpha_W');
plot(SUMMARY.alpha_w, SUMMARY.cval_ddra, '-o', 'DisplayName','DDRA');
plot(SUMMARY.alpha_w, SUMMARY.cval_gray, '-s', 'DisplayName','GraySeq');
xlabel('\alpha_W'); ylabel('Containment on validation (%)'); legend('Location','best');

nexttile; hold on; grid on; title('Conservatism proxy vs \alpha_W');
plot(SUMMARY.alpha_w, SUMMARY.sizeI_ddra, '-o', 'DisplayName','DDRA');
plot(SUMMARY.alpha_w, SUMMARY.sizeI_gray, '-s', 'DisplayName','GraySeq');
xlabel('\alpha_W'); ylabel('Aggregated interval size');

save_plot(f, plots_dir, sprintf('%s_fidcons_vs_alphaW', cfg.io.save_tag), ...
    'Formats', {'png','pdf'}, 'Resolution', 200);

% If present, also show shape-aware summaries vs α_W
dir_med  = coerce_numeric(getcol(SUMMARY,'dir_eps_med'));
dir_p90  = coerce_numeric(getcol(SUMMARY,'dir_eps_p90'));
hout_p90 = coerce_numeric(getcol(SUMMARY,'hout_p90'));
haus_med = coerce_numeric(getcol(SUMMARY,'haus_sym_med'));
if any(~isnan(dir_med) | ~isnan(hout_p90) | ~isnan(haus_med))
    f2 = figure('Color','w'); tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    nexttile; hold on; grid on;
    plot(SUMMARY.alpha_w, dir_med, '-o'); yline(1,'--');
    xlabel('\alpha_W'); ylabel('ratio'); title('dir\_eps median');
    nexttile; hold on; grid on;
    plot(SUMMARY.alpha_w, hout_p90, '-o');
    xlabel('\alpha_W'); ylabel('gap'); title('H_{out} p90');
    nexttile; hold on; grid on;
    plot(SUMMARY.alpha_w, haus_med, '-o');
    xlabel('\alpha_W'); ylabel('distance'); title('Hausdorff (sym) median');
    save_plot(f2, plots_dir, sprintf('%s_shapeaware_vs_alphaW', cfg.io.save_tag), ...
        'Formats', {'png','pdf'}, 'Resolution', 200);
end

close all force

% --- helper ---
function v = getcol(T, name)
    vars = string(T.Properties.VariableNames);
    idx  = find(strcmpi(vars, string(name)), 1);
    if ~isempty(idx), v = T.(vars(idx));
    else, v = nan(height(T),1);
    end
end
