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


## References


---

