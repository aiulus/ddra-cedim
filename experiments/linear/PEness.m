%% HOW TO USE — Persistency-of-Excitation Sweep (DDRA vs RCSI/Gray)
% What this does
%   Fixes dimension D and sweeps input shape/mode and PE order L to study effects on:
%     (i) fidelity (containment %)  (ii) conservatism (interval-size proxy)  (iii) runtime
%
% Equal-setting evaluation protocol 
%   • Shared datasets: The exact (x0,u,y) sequences are used for both DDRA and Gray.
%   • Unified noise policy: use_noise = shared.noise_for_gray && shared.noise_for_ddra
%       - When false, both methods operate with W=0 (comparable).
%   • Ridge guard: cfg.ddra.allow_ridge=false → rank-deficient cases are skipped (flagged).
%   • Metrics on the same points:
%       Containment: point-in-interval-hull of OUTPUT sets on VAL (identical for both).
%       Size proxy: aggregated interval width of OUTPUT sets across VAL (identical).
%   • Set reduction policy: common Girard reduction with the same order cap.
%
%% Key knobs
%   • Fix D via sweep_grid.D_list = 2 (for example).
%   • PE sweep: sweep_grid.pe_list = array of structs such as:
%       struct('mode','randn','order',L,'strength',1,'deterministic',true)
%       struct('mode','sinWave','order',L,...)  % for sinusoidal families
%   • Data budgets: cfg.shared.n_m, n_s, n_k (+ *_val for validation)
%
% Memory / IO toggles (recommended ON for large PE sweeps)
%   cfg.lowmem.gray_check_contain = false;   % skip heavy Gray checks (we use the light checker)
%   cfg.lowmem.store_ddra_sets    = false;   % streaming path for DDRA (lower memory)
%   cfg.lowmem.append_csv         = true;    % stream to CSV
%   cfg.lowmem.zonotopeOrder_cap  = 50;      % keep sets compact
%
%% Outputs
%   CSV + plots in experiments/results/{data,plots}/<save_tag>_sweeps


rng(1,'twister');

cfg = struct();
cfg.io = struct('save_tag','kMSD_pe_sweep');

% Shared
cfg.shared = struct();
cfg.shared.dyn   = "k-Mass-SD";
cfg.shared.type  = "standard";
cfg.shared.p_extr = 0.3;
cfg.shared.options_reach = struct('zonotopeOrder',100,'tensorOrder',2,'errorOrder',1,'tensorOrderOutput',2,'verbose',false);
cfg.shared.cs_base = struct('robustnessMargin',1e-9,'verbose',false,'cost',"interval",'constraints',"half");

% Data budgets
cfg.shared.n_m = 3; 
cfg.shared.n_s = 2;  
cfg.shared.n_k = 10;
cfg.shared.n_m_val = 2; 
cfg.shared.n_s_val = cfg.shared.n_s; 
cfg.shared.n_k_val = cfg.shared.n_k;

% DDRA noise
cfg.ddra = struct('eta_w',1,'alpha_w',0.01);
% --- DDRA ridge guard (defaults) ---
cfg.ddra.allow_ridge   = false;   % if false and rank-deficient -> skip point
cfg.ddra.lambda        = 1e-8;    % ridge lambda when allowed
cfg.ddra.ridge_gamma   = 1.0;     % scale for added uncertainty
cfg.ddra.ridge_policy  = "MAB";   % "MAB" (add generator to M_AB) or "W" (inflate W)

% Gray
cfg.gray = struct('methodsGray', ["graySeq"]);

rcsi_lbl = rcsi_label_from_cfg(cfg);                
cfg.io.save_tag = sprintf('%s_%s', cfg.io.save_tag, rcsi_lbl);  

% --- shared noise policy to keep DDRA/Gray comparable when not studying noise
cfg.shared.noise_for_gray = false;   % if false => Gray/RCSI runs with W = 0
cfg.shared.noise_for_ddra = true;    % DDRA uses W unless this is set to false
cfg.shared.use_noise = false;   % or true

% Sweep grid: two shapes × multiple orders
sweep_grid = struct();
sweep_grid.D_list        = 2;                    % fix dimension to isolate PE
sweep_grid.alpha_w_list  = cfg.ddra.alpha_w;
sweep_grid.n_m_list      = cfg.shared.n_m;
sweep_grid.n_s_list      = cfg.shared.n_s;
sweep_grid.n_k_list      = cfg.shared.n_k;
PE_orders = [1 2 3 4 5 6 7 8 9 10];
sweep_grid.pe_list = [ ...
    arrayfun(@(L) struct('mode','randn','order',L,'strength',1,'deterministic',true), PE_orders, 'uni',0), ...
    arrayfun(@(L) struct('mode','sinWave','order',L,'strength',1,'deterministic',true), PE_orders, 'uni',0) ...
];

% New: Memory efficiency toggles
cfg.lowmem = struct();
% Honor an optional cap on zonotopeOrder (keeps sets smaller during sweeps)
if isfield(cfg,'lowmem') && isfield(cfg.lowmem,'zonotopeOrder_cap')
    cfg.shared.options_reach.zonotopeOrder = min( ...
        cfg.shared.options_reach.zonotopeOrder, cfg.lowmem.zonotopeOrder_cap);
end
cfg.lowmem.gray_check_contain = false;   % don’t do expensive Gray containment
cfg.lowmem.store_ddra_sets    = false;   % don’t keep DDRA sets; compute metrics on the fly
cfg.lowmem.append_csv         = true;    % stream CSV row-by-row; don’t keep a giant table
cfg.lowmem.zonotopeOrder_cap  = 50;      % optional: lower order to shrink sets in memory

% === NEW: get plots_dir once ===
[plots_dir, ~] = init_io(cfg);

SUMMARY = run_sweeps(cfg, sweep_grid);

% --- Plots: fidelity/conservatism vs PE order, per shape
f1 = figure; tiledlayout(1,2); 
isRandn = strcmp(SUMMARY.pe_mode,'randn');
isSin   = strcmp(SUMMARY.pe_mode,'sinWave');

nexttile; hold on; grid on; title('Fidelity vs PE order');
plot(SUMMARY.pe_order(isRandn), SUMMARY.cval_ddra(isRandn), '-o','DisplayName','DDRA randn');
plot(SUMMARY.pe_order(isRandn), SUMMARY.cval_gray(isRandn), '-s','DisplayName','Gray randn');
plot(SUMMARY.pe_order(isSin),   SUMMARY.cval_ddra(isSin),   '-.^','DisplayName','DDRA sin');
plot(SUMMARY.pe_order(isSin),   SUMMARY.cval_gray(isSin),   '-v','DisplayName','Gray sin');
xlabel('PE order (L)'); ylabel('Containment on validation (%)'); legend('Location','best');

nexttile; hold on; grid on; title('Conservatism proxy vs PE order');
plot(SUMMARY.pe_order(isRandn), SUMMARY.sizeI_ddra(isRandn), '-o','DisplayName','DDRA randn');
plot(SUMMARY.pe_order(isRandn), SUMMARY.sizeI_gray(isRandn), '-s','DisplayName','Gray randn');
plot(SUMMARY.pe_order(isSin),   SUMMARY.sizeI_ddra(isSin),   '-.^','DisplayName','DDRA sin');
plot(SUMMARY.pe_order(isSin),   SUMMARY.sizeI_gray(isSin),   '-v','DisplayName','Gray sin');
xlabel('PE order (L)'); ylabel('Aggregated interval size');

% === NEW: save the first figure ===
save_plot(f1, plots_dir, 'pe_fidelity_conservatism', 'Formats', {'png','pdf'}, 'Resolution', 200);

%% TODO: This bit is irrelevant for the analysis / remove
% --- Runtime profile vs PE order
%f2 = figure; tiledlayout(1,2);
%nexttile; hold on; grid on; title('DDRA runtime vs L');
%plot(SUMMARY.pe_order, SUMMARY.t_ddra_learn,'-o','DisplayName','learn');
%plot(SUMMARY.pe_order, SUMMARY.t_ddra_check,'-s','DisplayName','check');
%plot(SUMMARY.pe_order, SUMMARY.t_ddra_infer,'-^','DisplayName','infer');
%xlabel('PE order (L)'); ylabel('s'); legend('Location','best');

%nexttile; hold on; grid on; title('Gray runtime vs L');
%plot(SUMMARY.pe_order, SUMMARY.t_gray_learn,'-o','DisplayName','learn');
%plot(SUMMARY.pe_order, SUMMARY.t_gray_val,  '-s','DisplayName','validate');
%plot(SUMMARY.pe_order, SUMMARY.t_gray_infer,'-^','DisplayName','infer');
%xlabel('PE order (L)'); ylabel('s'); legend('Location','best');

% === NEW: save the second figure ===
%save_plot(f2, plots_dir, 'pe_runtime_profiles', 'Formats', {'png','pdf'}, 'Resolution', 200);

disp('PE sweep done.');

close all force
