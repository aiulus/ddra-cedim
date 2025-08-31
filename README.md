# Data-Driven Reachability Analysis: 
# A Contrastive Analysis of Direct and Indirect Methods

## Repository Structure

```plaintext
.
├── dddra/
├── examples/             
├── experiments/            
├── rcsi/               
├── utils/                   # Helper functions
│   ├── io/                  # I/O functions
├── .gitignore
└── README.md                

```


# DDRA vs RCSI/Gray — Experiments Quick Guide

This repo compares Data-Driven Reachability Analysis and Reachset-Conformant System Identification under controlled sweeps. The pipelines share the **same datasets**, use a **single noise policy**, and evaluate with **identical metrics** on the **same validation points** for like-for-like results.

## What’s here

* `experiments/linear/sample_size.m` — Sweep number of input trajectories **n\_m** at fixed system dimension.
* `experiments/linear/PEness.m` — Persistency-of-excitation (PE) sweep: input shape + order **L**.
* `experiments/linear/scalability.m` — State-dimension sweep on the k-MSD chain.
* `experiments/run_sweeps.m` — Core runner (CSV streaming, memory toggles, shared VAL).

Results land in:

```
experiments/
  results/
    data/<save_tag>_sweeps/summary.csv
    plots/<save_tag>_sweeps/*.png|pdf
```

## Requirements

* MATLAB (R2022b+ recommended)
* CORA 2024+ (interval membership + `testCase` APIs).
  Heads-up: CORA ≥ 2024 deprecates `isprop` checks on contSets and some options.

## Quick start (copy/paste)

Minimal block you can adapt in any script:

```matlab
rng(1,'twister');
cfg = struct();
cfg.io = struct('save_tag','my_experiment');

% System & shared options
cfg.shared = struct();
cfg.shared.dyn = "k-Mass-SD";     % dynamics family
cfg.shared.type = "standard";     % same uncertainty preset across D
cfg.shared.p_extr = 0.3;
cfg.shared.options_reach = struct('zonotopeOrder',100,'tensorOrder',2, ...
    'errorOrder',1,'tensorOrderOutput',2,'verbose',false);
cfg.shared.cs_base = struct('robustnessMargin',1e-9,'verbose',false, ...
    'cost',"interval",'constraints',"half");

% Data budgets (train + val)
cfg.shared.n_m = 3; cfg.shared.n_s = 5; cfg.shared.n_k = 10;
cfg.shared.n_m_val = 2; cfg.shared.n_s_val = cfg.shared.n_s; cfg.shared.n_k_val = cfg.shared.n_k;

% DDRA options (incl. ridge guard)
cfg.ddra = struct('eta_w',1,'alpha_w',0.01, ...
    'allow_ridge',false,'lambda',1e-8,'ridge_gamma',1.0,'ridge_policy',"MAB");

% Gray method(s)
cfg.gray = struct('methodsGray', ["grayLS"]);  % or ["graySeq"]

% Unified noise policy (set both true or both false for comparability)
cfg.shared.noise_for_gray = false;
cfg.shared.noise_for_ddra = true;  % -> use_noise=false overall
cfg.shared.use_noise = false;      % legacy field; not required

% Memory/IO toggles
cfg.lowmem = struct('gray_check_contain',true,'store_ddra_sets',true, ...
                    'append_csv',true,'zonotopeOrder_cap',50);

% Sweep grid (example: sample size)
sweep_grid = struct();
sweep_grid.D_list        = [2];
sweep_grid.alpha_w_list  = cfg.ddra.alpha_w;
sweep_grid.n_m_list      = [2 4 8 16 32];
sweep_grid.n_s_list      = cfg.shared.n_s;
sweep_grid.n_k_list      = cfg.shared.n_k;
sweep_grid.pe_list       = {struct('mode','randn','order',2,'deterministic',true,'strength',1)};

SUMMARY = run_sweeps(cfg, sweep_grid);
```

Then plot from `SUMMARY` (examples already in each script).

## What the CSV contains

Per sweep point:

* **Fidelity**: `cval_gray`, `cval_ddra` (containment % on VAL).
* **Conservatism**: `sizeI_gray`, `sizeI_ddra` (aggregated output interval width).
* **Runtime**: `t_*` fields: learn / check (validation) / infer, plus totals (see `ensure_time_totals`).
* **Diagnostics**: `rankZ`, `condZ`, PE labels (`pe_mode`, `pe_order`), `use_noise`.

## Script-specific notes

* **`sample_size.m`**: vary `sweep_grid.n_m_list`; keep `pe_list` fixed for clarity.
* **`PEness.m`**: set `sweep_grid.pe_list` to arrays of structs; combine shapes (`randn`, `sinWave`) and orders `L`.
* **`scalability.m`**: vary `sweep_grid.D_list`; keep budgets fixed; use `randn` inputs.

## Call Hierarchy Details
Main path
```
└─ run_sweeps(cfg, grid)
   ├─ init_io(cfg) → plots_dir, results_dir
   ├─ init_sweep_axes(cfg, grid) → axes, baseC
   └─ FOR each (D, α_w, n_m, n_s, n_k, pe) in axes:
       ├─ (seeding) rng(row_seed,'twister');     % row_seed = 10000+row_index
       ├─ build_true_system(C) → sys_cora (linearSysDT), sys_ddra, R0, U
       │    • R0: state-space uncertainty zonotope around nominal x0
       │    • sys_cora: the object passed to CORA (RCSI/Gray side)
       │
       ├─ R0_cora = zonotope(0, R0.G)            % zero-centered copy for CORA
       ├─ pe_eff = finalize_pe_order(pe, sys_cora, C)
       │    • Policy: “explicit” in PEness; “minimal” elsewhere (nx+lag+1) ∧ feasible
       │
       ├─ [DDRA data gen / learn]  (see DDRA path)
       ├─ [RCSI/Gray identify + evaluate] (see Gray path)
       ├─ (optional) plot_reach_all_onepanel(...)   % unified comparison panel (online/both)
       ├─ Pack row → pack_row(...) → CSV
       └─ (optional) save artifacts row_####.mat     % sys_gray/sys_ddra/VAL/W_eff/M_AB/meta

```
DDRA path
```
DDRA: data generation & gating
└─ ddra_generate_data(C, sys_ddra, dt, R0, U, pe_eff)
   → Xminus, Uminus, Xplus, W, Zinfo, DATASET
     • genPEInput(...) inside (your latest), uses:
       – mode (randn/sinWave), L = pe_eff.order, n_u, n_k, dt
       – opts.strength/deterministic/contInput
       – hankel_stats() for PE diagnostics
     • DATASET holds x0_blocks{b}, u_blocks{b} (blocks of length n_k)

└─ testSuite_fromDDRA(sys_cora, R0, DATASET, n_k)
   → TS_train  % cell of testCase (y,u,x0) with measured sequences

└─ check_PE_order(TS_train, Lgate)  % gate row if PE not satisfied
   – In PEness: Lgate = requested L (explicit policy)
   – In others: Lgate = pe_eff.order (minimal sufficient)

DDRA: learning
└─ ddra_learn_Mab(Xminus, Uminus, Xplus, W, Zinfo, sys_ddra, C)
   → M_AB, ridgeInfo, W_eff
     • ridge guard & optional inflation of uncertainty
     • W_eff is **state-space** additive uncertainty (after E)

DDRA: validation & size on VAL
└─ VAL = VAL_from_TS(TS_val, DATASET_val)   % exact (x0,u,y)
└─ if store_ddra_sets:
   ├─ ddra_infer(sys_ddra, R0, U, W_used, M_AB, C, VAL) → {X_k}
   ├─ map_to_output(X_k, C_true) → {Y_k}   % or identity if C=I
   ├─ sizeI_ddra = mean_k sum|generators(Y_k)|
   └─ cval_ddra = containsY_on_VAL({Y_k}, VAL, tol)

   else (streaming):
   └─ ddra_infer_size_streaming(...) → (~, cval_ddra, widths_k) → sizeI_ddra
```
RCSI/Gray-Box path
```
RCSI (Gray) identify (shared data, unified options)
└─ optTS = ts_options_from_pe(C, pe, sys_cora)   % TS generation options
└─ gray_identify(sys_cora, R0_cora, U, C, pe, ...
       'overrideW', W_for_gray, 'options_testS', optTS,
       'externalTS_train', TS_train, 'externalTS_val', TS_val)
   → configs (cell), each with:
     • configs{i}.sys     : linearSysDT (identified Gray model)
     • configs{i}.params  : struct (R0,U, possibly W)
     • configs{i}.options : options + cs block (for conformance solver)

   NOTE: You already zero-centered R0 for CORA (R0_cora), which avoids
         double-counting initial offsets in validateReach().

Where W_for_gray is set (this was the bug you hit):
└─ In run_sweeps before calling gray_*:
   if use_noise
       W_for_gray = zonotope(zeros(sys_cora.nrOfDisturbances,1));
   else
       W_for_gray = zonotope(zeros(sys_cora.nrOfDisturbances,1));
   end
   % Never set W_for_gray = W_eff (wrong space!).

RCSI containment (fidelity) on VAL
└─ gray_containment(configs, sys_cora, R0_cora, U, C, pe_eff, ...)
   ├─ accum_suite(TS)
   │  └─ For each testCase tc:
   │     • If tc has multiple samples (size(u,3)>1) → split_tc_samples(tc)
   │       (forces validateReach to work on y, not y_a zero-input fallback)
   │     • [~, ev] = validateReach(tc_single, configs, CC=true, PS)
   │       (inside CORA, for each config:
   │        - strip options.cs before reach()
   │        - build params: R0+initialState, u', tFinal, pad inputs)
   │       ev.num_out accumulates mis-containments per config)
   ├─ ctrain = 100*(1 - num_out_train(i_gray)/num_all_train)
   └─ cval   = 100*(1 - num_out_val(i_gray)/num_all_val)

RCSI size proxy (conservatism) on VAL
└─ gray_infer_size_on_VAL(sys= configs{idxGray}.sys, TS_val, C, params_in, 'overrideW', W_for_gray)
   For each testCase tc in TS_val:
   ├─ params.R0 = params_in.R0 + tc.initialState
   ├─ params.u  = tc.u' (pad if needed), params.tFinal = sys.dt*(nk-1)
   ├─ if sys.nrOfDisturbances>0:
   │     params.W = W_for_gray    % **disturbance space** zonotope
   ├─ options = C.shared.options_reach (strip 'cs' if present)
   ├─ R = reach(sys, params, options); R = R.timePoint.set;
   └─ Map to outputs & accumulate size:
       Y_k = linearMap(R{k}, sys.C)   % set first, then matrix
       sizeI_gray += sum(abs(generators(Y_k)))
   → sizeI_gray = mean over all timepoints across all tcs

```

## References


---

