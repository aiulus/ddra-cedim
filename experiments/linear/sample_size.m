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
cfg.io = struct('save_tag', 'kMSD_sample_size_trajlength_sweep');

% System & shared options
cfg.shared = struct();
cfg.shared.dyn = "k-Mass-SD";
cfg.shared.type = "standard";      % same uncertainty preset across D
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
cfg.shared.n_m = 10;    % input traj count
cfg.shared.n_s = 5;   % samples per traj
cfg.shared.n_k = 10;    % horizon (train)

cfg.shared.n_m_val = 5;    % val traj count
cfg.shared.n_s_val = 5;
cfg.shared.n_k_val = 5;

% DDRA noise setup (fixed)
cfg.ddra = struct();
cfg.ddra.eta_w   = 1;      % number of W generators
cfg.ddra.alpha_w = 1.50;   % W scale
% --- DDRA ridge guard (defaults) ---
cfg.ddra.allow_ridge   = true;   % if false and rank-deficient -> skip point
cfg.ddra.lambda        = 1e-8;    % ridge lambda when allowed
cfg.ddra.ridge_gamma   = 1.0;     % scale for added uncertainty
cfg.ddra.ridge_policy  = "MAB";   % "MAB" (add generator to M_AB) or "W" (inflate W)
% --- choose which DDRA variant to run ---
cfg.ddra.variant = "std";      % "std" (standard) | "meas" (measurement-noise-aware)

% (optional) inject state measurement noise during TRAIN/VAL generation
cfg.data.train.meas = struct('enable', true, 'gen_scale', 0.01);  % or provide V directly
cfg.data.val.meas   = struct('enable', true, 'gen_scale', 0.01);  % symmetric VAL


% Gray methods (keep simple/fast)
cfg.gray = struct();
cfg.gray.methodsGray = ["graySeq"];

rcsi_lbl = rcsi_label_from_cfg(cfg);                
cfg.io.save_tag = sprintf('%s_%s', cfg.io.save_tag, rcsi_lbl);  

% --- shared noise policy to keep DDRA/Gray comparable when not studying noise
cfg.shared.noise_for_gray = true;   % false <=> W = 0
cfg.shared.noise_for_ddra = true;   % false <=> W = 0
cfg.shared.use_noise = true;   % or true
cfg.metrics.enhanced = true; % computes a wider range of evaluation metrics
%% TODO - wire as a bool for the extended plots below
cfg.metrics.directional = struct('enable', true, 'Nd', 64, 'seed', 12345);

% Optional safety verification task to measure conservatism
cfg.metrics.safety = struct( ...
  'enable', false, ...
  'H', ones(2,1), ...          % (#ineq x ny)
  'h', 10*ones(2,1), ...          % (#ineq x 1)
  'k_set', [], ...     % default = 1:n_k_val
  'mode', 'per_block', ... % or 'per_step'
  'tau_sweep', linspace(0, 0.2, 21), ... % spec tightening in output units
  'reach_reduce', 80 ... % optional zonotope order for true/gray reach during safety eval
);

% ---------- Sweep grid ----------
sweep_grid = struct();
sweep_grid.D_list       = 2;
sweep_grid.alpha_w_list = cfg.ddra.alpha_w;  % keep W fixed
sweep_grid.n_m_list = 10;
sweep_grid.n_s_list = 4;
sweep_grid.n_k_list = 4:4:100;
sweep_grid.pe_list = { struct('mode','randn','order', 4, 'strength',1,'deterministic',true) };

% New: Memory efficiency toggles
cfg.lowmem = struct();
cfg.lowmem.gray_check_contain = true;   % don’t do expensive Gray containment
cfg.lowmem.store_ddra_sets    = false;   % don’t keep DDRA sets; compute metrics on the fly
cfg.lowmem.append_csv         = true;    % stream CSV row-by-row; don’t keep a giant table
cfg.lowmem.zonotopeOrder_cap  = 25;      % optional: lower order to shrink sets in memory

% --- plotting mode + knobs
cfg.io.plot_mode     = "offline";   % "online" | "offline" | "both"
cfg.io.make_reach_plot = true;      % online: produce figures during run
cfg.io.plot_rows       = [1];       % which sweep rows to plot 
cfg.io.plot_dims       = [1 2];     % output dims
cfg.io.plot_every_k    = 1;         % plot every k-th step (declutter)
cfg.io.save_artifacts  = false;      % offline: keep .mat files for peeking/plotting
cfg.io.plot_mode = "offline";              % keep headless unless you opt in

cfg.io.base_dir = fileparts(fileparts(mfilename('fullpath'))); % or hard-code
cfg.allow_parallel = false;  % keep serial

cfg.shared.pe_min_policy = 'none';   % honor the L's in PE_orders exactly
cfg.shared.pe_verbose    = true;     % (optional) show requested->effective L

cfg.data = struct();
cfg.data.inject_proc_train = true;      % add process noise during TRAIN sims
cfg.data.inject_proc_val   = true;      % add process noise during VAL sims
cfg.data.proc_W_scale      = 1.0;       % scale *relative to W_used* (or give explicit set)
cfg.data.inject_meas_noise = false;     % optional: add measurement noise V
cfg.data.V_meas            = [];        % if [], no meas noise; else zonotope(cV, GV)

%% ---------- Run ----------
SUMMARY = run_sweeps(cfg, sweep_grid);
SUMMARY = ensure_time_totals(SUMMARY);

disp([SUMMARY.cval_gray SUMMARY.cval_ddra])
disp([SUMMARY.sizeI_gray SUMMARY.sizeI_ddra])

%% ---------- Visualization: auto-detect axis (n_m / n_k / n_s) ----------
vary_nm = numel(unique(sweep_grid.n_m_list)) > 1;
vary_nk = numel(unique(sweep_grid.n_k_list)) > 1;
vary_ns = numel(unique(sweep_grid.n_s_list)) > 1;
assert(sum([vary_nm,vary_nk,vary_ns])==1, ...
    'Expected exactly one of n_m, n_k, n_s to vary.');

if vary_nm
    var_axis = 'n_m';
    x        = coerce_numeric(SUMMARY.n_m);
    xlab     = 'n_m (distinct input trajectories)';
    lg_info  = sprintf('(n_k=%g, n_s=%g)', sweep_grid.n_k_list(1), sweep_grid.n_s_list(1));
elseif vary_nk
    var_axis = 'n_k';
    x        = coerce_numeric(SUMMARY.n_k);
    xlab     = 'n_k (trajectory length)';
    lg_info  = sprintf('(n_m=%g, n_s=%g)', sweep_grid.n_m_list(1), sweep_grid.n_s_list(1));
else
    var_axis = 'n_s';
    x        = coerce_numeric(SUMMARY.n_s);
    xlab     = 'n_s (samples per trajectory)';
    lg_info  = sprintf('(n_m=%g, n_k=%g)', sweep_grid.n_m_list(1), sweep_grid.n_k_list(1));
end

% Header bits: system, dimension, PE-order
sys_name = char(cfg.shared.dyn);
D        = sweep_grid.D_list(1);
L = NaN;
try
    if iscell(sweep_grid.pe_list),  L = sweep_grid.pe_list{1}.order; end
    if isstruct(sweep_grid.pe_list) && isfield(sweep_grid.pe_list,'order')
        L = sweep_grid.pe_list(1).order;
    end
catch, end
hdr = sprintf('%s (D=%d, PE-order=%g)', sys_name, D, L);

% Colors & label tag
colors = struct('ddra',[0.23 0.49 0.77],'gray',[0.85 0.33 0.10]);
ddra_lbl = ['DDRA ' lg_info];
rcsi_lbl_full = ['RCSI-' rcsi_lbl ' ' lg_info];


%% ---------- Runtime panels (total / learn / validation / inference) ----------
f = figure('Name','Runtime panels','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% 1) TOTAL
nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_total, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName',ddra_lbl);
plot(x, SUMMARY.t_gray_total, '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName',rcsi_lbl_full);
xlabel(xlab); ylabel('Seconds');
title(sprintf('%s — Total runtime vs %s', hdr, var_axis)); legend('Location','best');

% 2) LEARNING
nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_learn, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName',ddra_lbl);
plot(x, SUMMARY.t_gray_learn, '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName',rcsi_lbl_full);
xlabel(xlab); ylabel('Seconds');
title(sprintf('%s — Learning runtime', hdr)); legend('Location','best');

% 3) VALIDATION
nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_check, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName',ddra_lbl);
plot(x, SUMMARY.t_gray_val,   '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName',rcsi_lbl_full);
xlabel(xlab); ylabel('Seconds');
title(sprintf('%s — Validation/Check runtime', hdr)); legend('Location','best');

% 4) INFERENCE
nexttile; hold on; grid on;
plot(x, SUMMARY.t_ddra_infer, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName',ddra_lbl);
plot(x, SUMMARY.t_gray_infer, '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName',rcsi_lbl_full);
xlabel(xlab); ylabel('Seconds');
title(sprintf('%s — Inference runtime', hdr)); legend('Location','best');

[plots_dir, ~] = init_io(cfg);
save_plot(f, plots_dir, sprintf('runtime_panels_vs_%s_%s', var_axis, rcsi_lbl), ...
    'Formats', {'png','pdf'}, 'Resolution', 200);

%% ---------- Fidelity / Conservatism panels ----------
cval_ddra  = coerce_numeric(SUMMARY.cval_ddra);
cval_gray  = coerce_numeric(SUMMARY.cval_gray);
sizeI_ddra = coerce_numeric(SUMMARY.sizeI_ddra);
sizeI_gray = coerce_numeric(SUMMARY.sizeI_gray);

f2 = figure('Name','Fidelity & Conservatism','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% (1) Fidelity
nexttile; hold on; grid on;
plot(x, cval_ddra, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName',ddra_lbl);
plot(x, cval_gray, '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName',rcsi_lbl_full);
xlabel(xlab); ylabel('Containment on validation (%)');
title(sprintf('%s — Fidelity vs %s', hdr, var_axis)); legend('Location','best');

% (2) Conservatism
nexttile; hold on; grid on;
plot(x, sizeI_ddra, '-o', 'Color',colors.ddra, 'LineWidth',1.6, 'DisplayName',ddra_lbl);
plot(x, sizeI_gray, '-s', 'Color',colors.gray, 'LineWidth',1.6, 'DisplayName',rcsi_lbl_full);
xlabel(xlab); ylabel('Aggregated interval size (proxy)');
title(sprintf('%s — Conservatism vs %s', hdr, var_axis)); legend('Location','best');

save_plot(f2, plots_dir, sprintf('fidelity_conservatism_vs_%s_%s', var_axis, rcsi_lbl), ...
    'Formats', {'png','pdf'}, 'Resolution', 200);

close all force

% Reload paths for reading per-step logs
[plots_dir, results_dir] = init_io(cfg);

%% ---------- Shape-aware metrics across sweep (NEW) ----------
dir_med = coerce_numeric(getcol(SUMMARY,'dir_eps_med'));
dir_p90 = coerce_numeric(getcol(SUMMARY,'dir_eps_p90'));
hout_p90= coerce_numeric(getcol(SUMMARY,'hout_p90'));
haus_med= coerce_numeric(getcol(SUMMARY,'haus_sym_med'));
mw_g    = coerce_numeric(getcol(SUMMARY,'mw_gray_mean'));
mw_d    = coerce_numeric(getcol(SUMMARY,'mw_ddra_mean'));

f3 = figure('Name','Shape-aware metrics','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% (1) Directional ratio (Gray vs True)
nexttile; hold on; grid on;
plot(x, dir_med, '-o', 'LineWidth',1.6, 'DisplayName','dir\_eps median');
plot(x, dir_p90, '-s', 'LineWidth',1.6, 'DisplayName','dir\_eps p90');
yline(1,'--'); xlabel(xlab); ylabel('ratio'); 
title(sprintf('%s — Directional ratio (Gray/True)', hdr)); legend('Location','best');

% (2) Hausdorff_outer p90 (one-sided gap)
nexttile; hold on; grid on;
plot(x, hout_p90, '-o', 'LineWidth',1.6, 'DisplayName','H_{out} p90');
xlabel(xlab); ylabel('support gap'); 
title(sprintf('%s — H_{out} (p90)', hdr)); legend('Location','best');

% (3) Symmetric Hausdorff (median)
nexttile; hold on; grid on;
plot(x, haus_med, '-o', 'LineWidth',1.6, 'DisplayName','Hausdorff (sym) median');
xlabel(xlab); ylabel('distance'); 
title(sprintf('%s — Symmetric Hausdorff', hdr)); legend('Location','best');

% (4) Mean width (Gray vs DDRA)
nexttile; hold on; grid on;
plot(x, mw_g, '-s', 'LineWidth',1.6, 'DisplayName','Gray mean width');
plot(x, mw_d, '-o', 'LineWidth',1.6, 'DisplayName','DDRA mean width');
xlabel(xlab); ylabel('mean width'); 
title(sprintf('%s — Mean width (PRE-output)', hdr)); legend('Location','best');

save_plot(f3, plots_dir, sprintf('shape_aware_vs_%s_%s', var_axis, rcsi_lbl), ...
    'Formats', {'png','pdf'}, 'Resolution', 200);


%% ---------- Per-step panels for a selected sweep row (NEW) ----------
% Read per-step CSV if it exists
perstep_path = fullfile(results_dir, 'summary_perstep.csv');
if exist(perstep_path,'file')
    P = readtable(perstep_path);

    % Map SUMMARY rows -> per-step 'row' indices
    SUMMARY.row = (1:height(SUMMARY))';
    % Choose a representative row to visualize (last x by default)
    [~,ix] = max(x); 
    row_pick = SUMMARY.row(ix);
    K = P(P.row==row_pick, :);

    f4 = figure('Name','Per-step width, coverage, and ratio','Color','w');
    tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

    % (1) Width vs step
    nexttile; hold on; grid on;
    plot(K.k, K.wid_ddra, '-o', 'LineWidth',1.6, 'DisplayName','DDRA');
    plot(K.k, K.wid_gray, '-s', 'LineWidth',1.6, 'DisplayName','Gray');
    xlabel('step k'); ylabel('interval width (sum)'); 
    title(sprintf('%s — Per-step width @ %s=%g', hdr, var_axis, x(ix))); legend('Location','best');

    % (2) Coverage vs step
    nexttile; hold on; grid on;
    plot(K.k, K.cov_ddra, '-o', 'LineWidth',1.6, 'DisplayName','DDRA');
    plot(K.k, K.cov_gray, '-s', 'LineWidth',1.6, 'DisplayName','Gray');
    xlabel('step k'); ylabel('containment (%)'); 
    title(sprintf('%s — Per-step coverage @ %s=%g', hdr, var_axis, x(ix))); legend('Location','best');

    % (3) Gray/True size ratio vs step
    nexttile; hold on; grid on;
    if any(strcmpi(K.Properties.VariableNames,'ratio_gray_true'))
        plot(K.k, K.ratio_gray_true, '-d', 'LineWidth',1.6, 'DisplayName','Gray/True');
        xlabel('step k'); ylabel('ratio'); 
        title(sprintf('%s — Per-step size ratio @ %s=%g', hdr, var_axis, x(ix)));
        yline(1,'--');
        legend('Location','best');
    else
        text(0.5,0.5,'ratio\_gray\_true not available','HorizontalAlignment','center'); axis off
    end

    save_plot(f4, plots_dir, sprintf('perstep_row%03d_%s_%s', row_pick, var_axis, rcsi_lbl), ...
        'Formats', {'png','pdf'}, 'Resolution', 200);
end

function v = getcol(T, name)
    % Case-insensitive column lookup that returns a NaN column if missing
    vars = string(T.Properties.VariableNames);
    idx  = find(strcmpi(vars, string(name)), 1);
    if ~isempty(idx)
        v = T.(vars(idx));
    else
        v = nan(height(T),1);
    end
end


%% ---------- Artifact-derived per-step plots ----------
artifact_dir = fullfile(results_dir, 'artifacts');
if exist(artifact_dir,'dir')
    % Use the same row as earlier (last x by default). If you want another row,
    % set row_pick manually here.
    if ~exist('row_pick','var')
        SUMMARY.row = (1:height(SUMMARY))';
        [~,ix] = max(x);
        row_pick = SUMMARY.row(ix);
    end

    art_path = find_artifact_by_row(artifact_dir, row_pick);
    if ~isempty(art_path)
        % Base name for saving
        savebase = fullfile(plots_dir, sprintf('row%04d_artifacts_%s_%s', row_pick, var_axis, rcsi_lbl));

        % (A) Per-step fidelity & sizes from artifact
        plot_perstep_fidelity_sizes(art_path, ...
            'Dims', cfg.io.plot_dims, ...
            'Reduce', getfielddef(cfg.lowmem,'zonotopeOrder_cap', 60), ...
            'SaveBase', savebase, ...
            'Show', true);

        % (B) Overlays of reachable sets at selected steps
        nkv = unique(coerce_numeric(SUMMARY.n_k)); 
        if numel(nkv)==1, nkv = nkv(1); else, nkv = cfg.shared.n_k_val; end
        Klist = unique([1, round(nkv/2), nkv]);  % 1, mid, last
        overlay_reach_sets(art_path, ...
            'Dims', cfg.io.plot_dims, ...
            'Klist', Klist, ...
            'Reduce', getfielddef(cfg.lowmem,'zonotopeOrder_cap', 60), ...
            'SaveBase', savebase, ...
            'Show', true);
    else
        warning('No artifact .mat found for row %d in %s. Skipping artifact-derived plots.', row_pick, artifact_dir);
    end
else
    warning('Artifact directory not found: %s. Set cfg.io.save_artifacts=true if needed.', artifact_dir);
end

% --------- helpers (local to this script) ----------
function p = find_artifact_by_row(artifact_dir, rownum)
    % Try common patterns: row_####.mat, artifact_row_####.mat, etc.
    p = "";
    pat = sprintf('row_%04d', rownum);
    L = dir(fullfile(artifact_dir, '*.mat'));
    % Prefer exact row tag, otherwise fall back to first .mat
    for k = 1:numel(L)
        if contains(L(k).name, pat), p = string(fullfile(L(k).folder, L(k).name)); break; end
    end
    if p=="" && ~isempty(L)
        % Fallback: first .mat (last resort)
        p = string(fullfile(L(1).folder, L(1).name));
    end
end

% ---------- ROC (if safety enabled) ----------
if isfield(cfg,'metrics') && isfield(cfg.metrics,'safety') && cfg.metrics.safety.enable
    % Load any artifact for picked row, reusing previous 'row_pick' logic
    art_path = find_artifact_by_row(fullfile(results_dir,'artifacts'), row_pick);
    if strlength(art_path)>0
        A = load(art_path);
        if isfield(A,'metrics') && isfield(A.metrics,'safety')
            SFT = A.metrics.safety;
            figure('Color','w','Name','ROC – Safety certification');
            hold on; grid on;
            [FPRd,TPRd] = deal(SFT.roc_ddra.FPR, SFT.roc_ddra.TPR);
            [FPRg,TPRg] = deal(SFT.roc_gray.FPR, SFT.roc_gray.TPR);
            plot(FPRd, TPRd, '-o', 'DisplayName', sprintf('DDRA (AUC=%.3f)', SFT.roc_ddra.AUC),'LineWidth',1.6);
            plot(FPRg, TPRg, '-s', 'DisplayName', sprintf('Gray (AUC=%.3f)', SFT.roc_gray.AUC),'LineWidth',1.6);
            plot([0 1],[0 1],'k--','DisplayName','chance');
            xlabel('FPR'); ylabel('TPR'); legend('Location','southeast');
            title(sprintf('%s — ROC over \\tau (spec tightening)', hdr));

            save_plot(gcf, plots_dir, sprintf('roc_%s_%s', var_axis, rcsi_lbl), ...
                'Formats', {'png','pdf'}, 'Resolution', 200);
        end
    end
end
