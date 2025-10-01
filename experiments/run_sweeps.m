function SUMMARY = run_sweeps(cfg, grid)
% RUN_SWEEPS  Batch experiments for DDRA vs Gray under controlled sweeps.
% Produces CSV summary in results_dir and plots in plots_dir (via callers).
    if ~getfielddef(cfg,'allow_parallel',false)
        try
            ps = parallel.Settings; ps.Pool.AutoCreate = false;
            p = gcp('nocreate'); if ~isempty(p), delete(p); end
        catch 
        end
    end

    [plots_dir, results_dir] = init_io(cfg);

    % -------- Sweep axes & base config -----------------------------------
    [axes, baseC] = init_sweep_axes(cfg, grid);

    % DDRA variant selector: "std" (default) or "meas"
    variant = string(getfielddef(getfielddef(cfg,'ddra',struct()), 'variant', "std"));
    
    % Precompute for in-memory mode only (IO helper will use this)
    Ntot = numel(axes.D) * numel(axes.alpha_w) * numel(axes.n_m) * ...
           numel(axes.n_s) * numel(axes.n_k) * numel(axes.pe);
    
    % Centralized I/O state (handles CSV init, per-step header, etc.)
    IO = sweepio_begin(cfg, plots_dir, results_dir, Ntot);
    t0_all = tic;
    NALL = Ntot;
    

    % -------- Low-memory / IO toggles (with safe defaults) ---------------
    LM = getfielddef(cfg, 'lowmem', struct());
    LM.gray_check_contain = getfielddef(LM, 'gray_check_contain', true);
    LM.store_ddra_sets    = getfielddef(LM, 'store_ddra_sets', true);
    LM.append_csv         = getfielddef(LM, 'append_csv', false);
    %LM.store_ddra_sets = false; 


    % -------- Sweep axes & base config -----------------------------------

    if ~isfield(cfg,'metrics') || ~isfield(cfg.metrics,'tol'), cfg.metrics.tol = 1e-6; end
    % ---- PE policy: explicit for "pe_sweep", minimal otherwise unless forced ----
    is_pe_sweep_tag = contains(lower(string(getfielddef(cfg.io,'save_tag',""))), 'pe_sweep');
    has_forced = any(cellfun(@(p) isfield(p,'force_order') && p.force_order, axes.pe)) ...
          || any(cellfun(@(p) isfield(p,'order_gate')   && ~isempty(p.order_gate), axes.pe));
    baseC.shared.pe_policy = getfielddef(baseC.shared,'pe_policy', ...
        ternary(is_pe_sweep_tag || has_forced, "explicit", "minimal"));


    function out = ternary(cond, a, b), if cond, out=a; else, out=b; end, end

    % Force explicit PE whenever there are multiple PE orders on the grid
    is_multi_pe = numel(axes.pe) > 1;
    if is_multi_pe, baseC.shared.pe_policy = "explicit"; end
    % ------------------------ Main sweeps --------------------------------
    for D = axes.D
      for alpha_w = axes.alpha_w
        for n_m = axes.n_m
          for n_s = axes.n_s
            for n_k = axes.n_k
              for ip = 1:numel(axes.pe)
                pe = axes.pe{ip};
                row_index = IO.rowi + 1;
                
                C = baseC;   % <-- make C first, then compute seed from C & pe
                tag = getfielddef(getfielddef(cfg,'io',struct()),'save_tag','');
                
                row_seed = stable_seed( ...
                    D, ...
                    alpha_w, ...
                    pe, ...                                
                    getfielddef(C.shared,'dyn',""), ...
                    getfielddef(C.shared,'type',""), ...
                    tag);
                
                rng(row_seed,'twister');


                % --- instantiate config for this run ---
                C = baseC;
                C.shared.seed = row_seed;
                C.shared.n_m = n_m;  C.shared.n_s = n_s;  C.shared.n_k = n_k;
                C.shared.n_m_val = max(2, min(n_m, getfielddef(baseC.shared, 'n_m_val', 2)));
                C.shared.n_s_val = n_s;  C.shared.n_k_val = n_k;
                C.ddra.alpha_w   = alpha_w;
            
                % Make nominal PE input reproducible and L-dependent without freezing RNG
                pe.deterministic = getfielddef(pe,'deterministic', true);
                pe.seed_base = uint32(mod( double(row_seed) + 131*double(getfielddef(pe,'order',0)) + 977*double(ip), 2^31-1 ));

                if C.shared.dyn == "k-Mass-SD"; C.shared.dyn_p = D; end
    
                % --- Build true systems ---
                [sys_cora, sys_ddra, R0, U] = build_true_system(C);
                % Zero-centered copy for CORA (deviation only)
                c0 = center(R0);
                G0 = R0.G;               
                R0_cora = zonotope(zeros(size(c0)), G0);   % center = 0, keep same generators

                use_noise   = resolve_use_noise(C.shared);
                % Respect policy: PEness -> explicit L; others -> minimal sufficient
                C.shared.pe_policy = baseC.shared.pe_policy;   % carry into row config

                %% Patch
                % Harmonize PE strength across scripts if strength_rel is provided
                if ~isfield(pe,'strength') || isempty(pe.strength)
                    if isfield(pe,'strength_rel') && ~isempty(pe.strength_rel)
                        G = generators(U); hw = sum(abs(G),2);                    % per-channel half-width
                        pe.strength = mean(pe.strength_rel .* (hw + (hw==0)));    % scalar strength for genPEInput
                    end
                end
                
                if ~isfield(pe,'order') || isempty(pe.order)
                    pe.order = max(1, min(sys_cora.nrOfDims+1, C.shared.n_k-1));
                end

                respect = (string(C.shared.pe_policy) == "explicit");
                pe_eff = pe_normalize(pe, U, sys_cora, C.shared.n_k, ...
                                      'respect_explicit_order', respect);

                % ================= DDRA =================
                t0 = tic;
                [Xminus, Uminus, Xplus, W, Zinfo, DATASET] = ddra_generate_data(C, sys_ddra, sys_cora.dt, R0, U, pe_eff, sys_cora);
                Tlearn = toc(t0);
    
                % Build TRAIN and VAL suites deterministically from generator
                TS_train = testSuite_fromDDRA(sys_cora, R0, DATASET, C.shared.n_k, C.shared.n_m, C.shared.n_s);

                if string(C.shared.pe_policy) == "explicit"
                    % Gate on the requested order when sweeping PE explicitly
                    Lgate = getfielddef(pe,'order_gate', getfielddef(pe,'order', pe_eff.order));
                else
                    Lgate = pe_eff.order; % Gate on the effective minimal sufficient order
                end

                % DEBUG STATEMENT
                fprintf('L requested=%g  effective=%g  policy=%s\n', ...
                    getfielddef(pe,'order',NaN), Lgate, string(C.shared.pe_policy));
               
                okPE = check_PE_order(TS_train, Lgate);

                
                if ~okPE
                    timers = struct('Tlearn', exist('Tlearn','var')*Tlearn + ~exist('Tlearn','var')*NaN, ...
                                    'Tcheck', NaN, 'Tinfer', NaN, 'Tlearn_g', NaN, ...
                                    'Tvalidate_g', NaN, 'Tinfer_g', NaN);
                    row = build_skip_row(C, D, alpha_w, pe, Zinfo, timers, use_noise, "PE_not_satisfied", {});
                    IO  = sweepio_write_row(IO, row);
                    clear TS_train DATASET
                    continue; % keep grid aligned
                end
    
                C_val = C;
                C_val.shared.n_m = C.shared.n_m_val;
                C_val.shared.n_s = C.shared.n_s_val;
                C_val.shared.n_k = C.shared.n_k_val;
                rng(row_seed+1,'twister');[~,~,~,~,~, DATASET_val] = ...
                        ddra_generate_data(C_val, sys_ddra, sys_cora.dt, R0, U, pe_eff, sys_cora);
                TS_val   = testSuite_fromDDRA(sys_cora, R0, DATASET_val, C_val.shared.n_k, C_val.shared.n_m, C_val.shared.n_s);
                fprintf('VAL suite: %d cases, n_k=%d (expects %d=m*s)\n', ...
                    numel(TS_val), size(TS_val{1}.y,1), C_val.shared.n_m*C_val.shared.n_s);

                % Learn M_AB
                % t1 = tic;
                % try
                %     % New signature with ridge control & W_eff (possibly inflated for ridge)
                %     [M_AB, ridgeInfo, W_eff] = ddra_learn_Mab(Xminus, Uminus, Xplus, W, Zinfo, sys_ddra, C);
                % catch
                %     % Old signature fallback (keeps older ddra_learn_Mab working)
                %     M_AB = ddra_learn_Mab(Xminus, Uminus, Xplus, W, Zinfo, sys_ddra);
                %     ridgeInfo = struct('used', false, 'skipped', false, ...
                %                        'lambda', 0, 'kappa', NaN, 'policy', "none");
                %     W_eff = W;  % no inflation in legacy path
                % end
                % Tcheck = toc(t1);
                
                % TODO: bring time measurement to just around ddra_learn_*
                t1 = tic;
                switch variant
                  case "std"
                    try
                      [M_AB, ridgeInfo, W_eff] = ddra_learn_Mab(Xminus, Uminus, Xplus, W, Zinfo, sys_ddra, C);
                    catch
                      M_AB = ddra_learn_Mab(Xminus, Uminus, Xplus, W, Zinfo, sys_ddra);
                      ridgeInfo = struct('used', false, 'skipped', false, 'lambda', 0, 'kappa', NaN, 'policy', "none");
                      W_eff = W;
                    end
                    learned_kind = "std";
                
                    case "meas"
                    [ABc, AV_one, V_meas, info_meas] = ddra_learn_measnoise(Xminus, Uminus, Xplus, W, sys_ddra, C);
                    % For a clean schema, emulate outputs:
                    M_AB = [];   % not used in this branch
                    ridgeInfo = struct('used', false, 'skipped', false, 'lambda', 0, 'kappa', NaN, 'policy', "measnoise");
                    W_eff = W;   % base W stays (can still be mapped/inflated later)
                    learned_kind = "meas";
                
                  otherwise
                    error('Unknown cfg.ddra.variant = %s', variant);
                end
                Tcheck = toc(t1);

                try
                    row.ab_ngen = size(M_AB.G,3); % #gens in matZonotope (if available)
                catch
                    row.ab_ngen = NaN;
                end
                
                % If rank-deficient and allow_ridge = false, skip this sweep point
                if isfield(ridgeInfo,'skipped') && ridgeInfo.skipped
                    timers = struct('Tlearn', exist('Tlearn','var')*Tlearn + ~exist('Tlearn','var')*NaN, ...
                                    'Tcheck', exist('Tcheck','var')*Tcheck + ~exist('Tcheck','var')*NaN, ...
                                    'Tinfer', NaN, 'Tlearn_g', NaN, 'Tvalidate_g', NaN, 'Tinfer_g', NaN);
                    extras = {'ddra_ridge',false, 'ddra_lambda',0, 'ddra_kappa',NaN, 'ddra_ridge_policy',"skip"};
                    row = build_skip_row(C, D, alpha_w, pe, Zinfo, timers, use_noise, "skip", extras);
                    IO  = sweepio_write_row(IO, row);
                    clear Xminus Xplus TS_train TS_val DATASET DATASET_val
                    continue;
                end

                
                clear Xminus Xplus  % free data blocks early

    
                % ---- unified disturbance policy ----
                if resolve_use_noise(C.shared)
                    % baseline from learner
                    W_used = W_eff;
                    % enforce alpha_W scale if the learner didn't
                    g_now  = zono_gnorm(W_used);
                    if ~isnan(g_now) && g_now > 0
                        W_used = zono_scale(W_used, C.ddra.alpha_w / g_now);
                    else
                        % if no usable scale, default to raw alpha_W
                        W_used = zono_scale(W_used, C.ddra.alpha_w);
                    end
                else
                    W_used = zonotope(zeros(size(center(W_eff),1),1));
                end
                
                % (diagnostics into the row/CSV so you can confirm scaling worked)
                row.w_eff_gmean  = zono_gnorm(W_eff);
                row.w_used_gmean = zono_gnorm(W_used);

                

                % ---- Build W_for_gray only if disturbance channels exist ----
                W_for_gray = build_W_for_gray(sys_cora, U, W_used);
                
            
                % ================= GRAY =================
                optTS = ts_options_from_pe(C, pe, sys_cora);
                t3 = tic;
                configs = gray_identify(sys_cora, R0_cora, U, C, pe, ...
                    'overrideW', W_for_gray, ...
                    'options_testS', optTS, ...
                    'externalTS_train', TS_train, ...
                    'externalTS_val',   TS_val);  
                Tlearn_g = toc(t3);
                
                idxGray = pick_gray_config(configs, C);
                  
                % ---- Set W_pred only if identified Gray sys has disturbance channels ----
                % W_pred = build_W_pred(configs{idxGray}.sys, U, W_used);
                W_pred = build_W_pred(configs{idxGray}.sys, C, W_used);  % map state-space W_used -> Gray disturbance space
                if isempty(W_pred)
                    W_pred = zonotope(zeros(configs{idxGray}.sys.nrOfDisturbances,1)); % safe fallback
                end
                configs{idxGray}.params.W = W_pred;


                VAL = pack_VAL_from_TS(TS_val);

                fprintf('VAL(flat) m*s=%d  nk=%d\n', numel(VAL.y), size(VAL.y{1},1));

                %% PATCH
                % --- hard invariants for evaluation parity ---
                must( strcmp(val_hash(VAL), val_hash(pack_VAL_from_TS(TS_val))), 'VAL mismatch (pack parity)');
                must( abs(configs{idxGray}.sys.dt - sys_cora.dt) < 1e-12, 'dt mismatch between Gray and true');
                
                % same U and x0 grids used on both sides
                must( numel(VAL.u) == C_val.shared.n_m*C_val.shared.n_s, 'VAL blocks count mismatch');
                for b = 1:numel(VAL.u)
                    must( isequal(size(VAL.u{b}), [C_val.shared.n_k, size(sys_cora.B,2)]), 'VAL.u block shape mismatch');
                    must( isequal(size(VAL.y{b}), [C_val.shared.n_k, size(sys_cora.C,1)]), 'VAL.y block shape mismatch');
                end
                
                pmode = lower(string(getfielddef(cfg.io,'plot_mode','offline')));
                % === Unified artifact save (after configs, VAL, M_AB, W_eff exist) ===
                sys_true_dt = normalize_to_linearSysDT(sys_ddra, sys_cora.dt);

                artifact = struct();
                artifact.sys_gray = configs{idxGray}.sys;
                artifact.sys_ddra = sys_true_dt;
                artifact.VAL      = VAL;
                artifact.W_eff    = W_used;
                artifact.meta     = struct('row', row_index, 'dt', sys_true_dt.dt, ...
                  'D', D, 'alpha_w', alpha_w, 'n_m', C.shared.n_m, 'n_s', C.shared.n_s, ...
                  'n_k', C.shared.n_k, 'pe', pe, 'save_tag', getfielddef(cfg.io,'save_tag',''), 'schema', 2);
                
                if variant == "std"
                    artifact.M_AB = M_AB;
                else
                    artifact.AB_meas_center = ABc;
                    artifact.AV_oneterm     = AV_one;
                    artifact.V_meas         = V_meas;
                end


                artifact.VAL.R0   = R0;
                artifact.VAL.U    = U;
                
                sweepio_save_artifact(IO, row_index, artifact);


                want_online = any(pmode == ["online","both"]) && getfielddef(cfg.io,'make_reach_plot',true);
                if want_online && ismember(row_index, getfielddef(cfg.io,'plot_rows',[1]))
                    dims   = getfielddef(cfg.io,'plot_dims',[1 2]);
                
                    sys_true_dt = normalize_to_linearSysDT(sys_ddra, sys_cora.dt);
                
                    W_plot = pick_W_for_plot(W_used, use_noise);

                    reach_dir = fullfile(plots_dir,'reach');
                    if ~exist(reach_dir,'dir'), mkdir(reach_dir); end
                    savename_png = fullfile(reach_dir, sprintf('row_%04d_y%d%d.png', row_index, dims(1), dims(2)));
                    
                    extraArgs = {};
                    if variant == "std"
                        extraArgs = {'MAB', M_AB};
                    else
                        extraArgs = {'Variant','meas', 'ABc', ABc, 'AV', AV_one, 'Vmeas', V_meas};
                    end
                    
                    plot_reach_all_onepanel( ...
                        sys_true_dt, ...
                        configs{idxGray}.sys, ...
                        VAL, ...
                        'W', pick_W_for_plot(W_used, use_noise), ...
                        'Dims', dims, ...
                        'ShowSamples', false, ...
                        'Save', savename_png, ...
                        'ReachOptions', getfielddef(C.shared,'options_reach', struct('zonotopeOrder',100,'reductionTechnique','girard')), ...
                        extraArgs{:});


                end

                % --- SAFETY EVAL (if enabled) ---
                SFT = getfielddef(getfielddef(cfg,'metrics',struct()), 'safety', struct('enable',false));
                if getfielddef(SFT,'enable',false)
                    Hs   = SFT.H; hs = SFT.h;
                    kset = getfielddef(SFT,'k_set', 1:C_val.shared.n_k_val);
                    tauv = getfielddef(SFT,'tau_sweep', linspace(0,0.2,21));
                    rord = getfielddef(SFT,'reach_reduce', 80);

                    [safe_true, safe_alg_tau] = eval_safety_labels( ...
                        sys_true_dt, configs{idxGray}.sys, VAL, R0, M_AB, W_used, W_pred, ...
                        Hs, hs, kset, tauv, rord);
                
                    ROC = roc_from_decisions(safe_true, safe_alg_tau);
                
                    % (logging stays the same)
                    i0 = 1;
                    row.sfty_tp_ddra = ROC.ddra.TP(i0); row.sfty_fp_ddra = ROC.ddra.FP(i0);
                    row.sfty_tn_ddra = ROC.ddra.TN(i0); row.sfty_fn_ddra = ROC.ddra.FN(i0);
                    row.sfty_tp_gray = ROC.gray.TP(i0); row.sfty_fp_gray = ROC.gray.FP(i0);
                    row.sfty_tn_gray = ROC.gray.TN(i0); row.sfty_fn_gray = ROC.gray.FN(i0);
                    row.roc_auc_ddra = ROC.ddra.AUC;
                    row.roc_auc_gray = ROC.gray.AUC;
                
                    if exist('artifact','var')
                        artifact.metrics.safety = struct( ...
                            'tau', tauv, ...
                            'safe_true', safe_true, ...
                            'dec_ddra', safe_alg_tau.ddra, ...
                            'dec_gray', safe_alg_tau.gray, ...
                            'roc_ddra', ROC.ddra, ...
                            'roc_gray', ROC.gray);
                        sweepio_save_artifact(IO, row_index, artifact);
                    end
                end


    
                % ---- Gray validation on EXACT same points (contains_interval metric)
                % CORA's validateReach (internally splits multi-sample testCases)
                do_plot = any(lower(string(getfielddef(cfg.io,'plot_mode','offline'))) == ["online","both"]);
                
                if do_plot
                    PS = struct();
                    PS.k_plot   = 1:C.shared.n_k_val;                     % <- plot all steps
                    PS.dims     = getfielddef(cfg.io,'plot_dims',[1 2]);  % <- 2D slice
                    PS.plot_Yp  = false;                                  % <- like the example
                    ps_args = {'plot_settings', PS};
                else
                    ps_args = {'plot_settings', []};  % truly suppress plotting
                end
                
                [ctrain_gray, cval_gray, Tvalidate_g, VAL] = gray_containment( ...
                                configs, sys_cora, R0_cora, U, C, pe_eff, ...
                                'externalTS_train', TS_train, ...
                                'externalTS_val',   TS_val,   ...
                                'check_contain',    true,     ...
                                ps_args{:});

    
                % ---- DDRA inference (stored or streaming) on SAME validation points
                if variant == "std"
                    % START PATCH
                    if LM.store_ddra_sets
                        t2 = tic;
                    
                        % same reduction policy as streaming
                        ordCap = getfielddef(getfielddef(C,'lowmem',struct()), 'zonotopeOrder_cap', 100);
                        Kred   = max(1, min(100, round(ordCap)));
                    
                        nkv        = C.shared.n_k_val;
                        m          = size(sys_true_dt.D,2);
                        Bv         = numel(VAL.x0);
                        wid_sums   = zeros(nkv,1);
                        wid_counts = zeros(nkv,1);
                    
                        % containment counters (exactly like streaming)
                        num_in = 0; num_all = 0; tol = getfielddef(cfg.metrics,'tol',1e-6);
                    
                        % PRE-update loop mirroring streaming
                        for b = 1:Bv
                            Xk = reduce(VAL.R0 + VAL.x0{b}, 'girard', Kred);  % PRE-update X_0
                            Uk = VAL.u{b};                                    % (n_k × m)
                            Yb = VAL.y{b};                                    % (n_k × ny)
                    
                            for ps = 1:nkv
                                Xk   = reduce(Xk, 'girard', Kred);
                                uk   = Uk(ps,:).';
                                U_pt = zonotope(uk);
                    
                                % PRE-update output at time k INCLUDING feedthrough
                                Yset = sys_true_dt.C * Xk + sys_true_dt.D * uk;
                    
                                % widths (SUM across outputs; canonicalize later)
                                Iv   = interval(Yset);
                                lo   = try_get(@() infimum(Iv), @() Iv.inf);
                                hi   = try_get(@() supremum(Iv), @() Iv.sup);
                                wvec = max(hi - lo, 0);
                                wid_sums(ps)   = wid_sums(ps)   + sum(double(wvec), 'omitnan');
                                wid_counts(ps) = wid_counts(ps) + 1;
                    
                                % containment in interval hull of Yset
                                y_meas = Yb(ps,:).';
                                if contains_interval(y_meas, Yset, tol), num_in = num_in + 1; end
                                num_all = num_all + 1;
                    
                                % POST-update propagation
                                Xk = M_AB * cartProd(Xk, U_pt) + W_used;
                            end
                        end
                    
                        Tinfer = toc(t2);
                    
                        % average across blocks; canonical report = MEAN across outputs
                        wid_ddra_k_raw = wid_sums ./ max(1, wid_counts);
                        [sizeI_ddra, wid_ddra_k] = normalize_widths(wid_ddra_k_raw, size(sys_true_dt.C,1), "mean");
                        sizeI_ddra_k = wid_ddra_k;
                    
                        cval_ddra = (num_all>0) * (100 * num_in / max(1,num_all));
                    else
                        t2 = tic;
                        [sizeI_ddra, cval_ddra, wid_ddra_k] = ddra_infer_size_streaming( ...
                            sys_ddra, R0, [], W_used, M_AB, C, VAL);
                        Tinfer = toc(t2);
                        [sizeI_ddra, wid_ddra_k] = normalize_widths(wid_ddra_k, size(sys_true_dt.C,1), "mean");
                        sizeI_ddra_k = wid_ddra_k;
                    end
                else
                    % measurement-noise-aware
                    t2 = tic;
                    [sizeI_ddra, cval_ddra, wid_ddra_k] = ddra_infer_size_streaming_meas( ...
                        sys_true_dt, R0, [], W_eff, ABc, AV_one, V_meas, C, VAL);
                    Tinfer = toc(t2);
                    sizeI_ddra_k = wid_ddra_k;
                end


                % ---- Gray size metric (interval proxy), consistent with DDRA
                t4 = tic;
                [sizeI_gray, wid_gray_k] = gray_infer_size_on_VAL(configs{idxGray}.sys, TS_val, C, ...
                    configs{idxGray}.params, 'overrideW', W_pred, 'width_agg', 'sum');
                Tinfer_g   = toc(t4);
                [sizeI_gray, wid_gray_k] = normalize_widths(wid_gray_k, size(configs{idxGray}.sys.C,1), "mean");

                cov_ddra_k = []; cov_gray_k = []; fv_ddra = []; fv_gray = [];

                if getfielddef(cfg.metrics,'enhanced', true)
                    % GRAY per-step metrics from VAL only
                    try
                        %% ACHTUNG - PLACEHOLDER - HASN'T BEEN IMPLEMENTED YET
                        [~, cov_gray_k, fv_gray, ~] = perstep_gray_metrics(configs{idxGray}.sys, VAL, W_pred, struct('zonotopeOrder',60));
                    catch
                        cov_gray_k = []; fv_gray = []; 
                    end
                
                    % DDRA per-step metrics from VAL only (uses M_AB; fast)
                    try
                        %% ACHTUNG - PLACEHOLDER - HASN'T BEEN IMPLEMENTED YET
                        [~, cov_ddra_k, fv_ddra, ~] = perstep_ddra_metrics(sys_true_dt, VAL, M_AB, W_used, 60);
                    catch
                        cov_ddra_k = []; fv_ddra = [];
                    end
                end

                % --- Shape-aware / directional metrics (optional; guarded, computed once) ---
                DIR = [];   % local struct to carry metrics forward; if empty -> not computed
                
                dir_cfg = getfielddef(getfielddef(cfg,'metrics',struct()), 'directional', struct());
                want_dir = getfielddef(dir_cfg,'enable', false);
                
                if want_dir
                    % Build reach options (never pass 'cs' here)
                    optReach = getfielddef(C.shared,'options_reach', struct('zonotopeOrder',60,'reductionTechnique','girard'));
                    if isfield(optReach,'cs'), optReach = rmfield(optReach,'cs'); end
                
                    nkv = C.shared.n_k_val;
                    ny  = size(sys_true_dt.C,1);
                    Nd  = getfielddef(dir_cfg,'Nd',64);
                    dseed = getfielddef(dir_cfg,'seed',12345);
                    Ddirs = sample_unit_dirs(ny, Nd, dseed);
                
                    % Accumulators (scalar + per-step)
                    dir_eps_all = []; dir_del_all = [];
                    hd_sym_k = zeros(nkv,1); mw_gray_k = zeros(nkv,1); mw_ddra_k = zeros(nkv,1);
                
                    % One representative PRE→POST chain for DDRA mean width
                    ordCap = getfielddef(getfielddef(C,'lowmem',struct()), 'zonotopeOrder_cap', 100);
                    Kred   = max(1, min(100, round(ordCap)));
                
                    for b = 1:numel(VAL.x0)
                        % Precompute VAL reaches once per block
                        params_true = struct('R0', R0, 'U', U, 'u', VAL.u{b}', 'tFinal', sys_true_dt.dt*(nkv-1));
                        params_gray = struct('R0', configs{idxGray}.params.R0 + VAL.x0{b}, ...
                                             'u',  VAL.u{b}', 'tFinal', configs{idxGray}.sys.dt*(nkv-1));
                        if isprop(configs{idxGray}.sys,'nrOfDisturbances') && configs{idxGray}.sys.nrOfDisturbances>0
                            params_gray.W = W_pred; % already in disturbance space
                        end
                
                        Rt = reach(sys_true_dt, params_true, optReach).timePoint.set;
                        Rg = reach(configs{idxGray}.sys, params_gray, optReach).timePoint.set;
                
                        % Representative chain for DDRA widths
                        Xk_rep = reduce(R0 + VAL.x0{b}, 'girard', Kred);
                
                        for k = 1:nkv
                            % True/Gray output sets as zonotopes
                            Yt_k = toZono(linearMap(Rt{k}, sys_true_dt.C));
                            Yg_k = toZono(linearMap(Rg{k}, configs{idxGray}.sys.C));
                
                            % Directional supports (outer metrics)
                            s_t = support_zono_vec(Ddirs, center(Yt_k), generators(Yt_k));
                            s_g = support_zono_vec(Ddirs, center(Yg_k), generators(Yg_k));
                
                            eps_kd = s_g ./ max(s_t, eps);        % ratio
                            del_kd = max(s_g - s_t, 0);           % outer gap
                
                            dir_eps_all = [dir_eps_all; eps_kd(:)];
                            dir_del_all = [dir_del_all; del_kd(:)];
                
                            % Symmetric Hausdorff on sampled directions
                            hd_sym_k(k) = max(abs(s_g - s_t));
                
                            % Gray mean width (via support sums)
                            mw_gray_k(k) = mean( s_g + support_zono_vec(-Ddirs, center(Yg_k), generators(Yg_k)) );
                
                            % DDRA PRE-update output width along representative chain
                            Xk_rep = reduce(Xk_rep, 'girard', Kred);
                            uk     = VAL.u{b}(k,:).';
                            Yd_k   = sys_true_dt.C*Xk_rep + sys_true_dt.D*uk;
                            Yd_k   = toZono(Yd_k);
                            mw_ddra_k(k) = mean( ...
                                support_zono_vec(Ddirs,  center(Yd_k), generators(Yd_k)) + ...
                                support_zono_vec(-Ddirs, center(Yd_k), generators(Yd_k)) );
                
                            % Advance PRE→POST for next step
                            if variant == "std"
                                Xk_rep = M_AB * cartProd(Xk_rep, zonotope(uk)) + W_used;
                            else
                                % meas-noise-aware propagation: add V_meas to the regressor for k>1
                                X_for_AB = Xk_rep;
                                if k > 1 && exist('V_meas','var') && ~isempty(generators(V_meas))
                                    X_for_AB = Xk_rep + V_meas;
                                end
                                Xk_rep = ABc * cartProd(X_for_AB, zonotope(uk)) + AV_one + W_used;
                            end

                        end
                    end
                
                    % Scalar summaries
                    DIR.eps_med      = median(dir_eps_all, 'omitnan');
                    DIR.eps_p90      = prctile(dir_eps_all, 90);
                    DIR.hout_med     = median(dir_del_all, 'omitnan');
                    DIR.hout_p90     = prctile(dir_del_all, 90);
                    DIR.haus_sym_med = median(hd_sym_k, 'omitnan');
                    DIR.haus_sym_p90 = prctile(hd_sym_k, 90);
                    DIR.mw_gray_mean = mean(mw_gray_k, 'omitnan');
                    DIR.mw_ddra_mean = mean(mw_ddra_k, 'omitnan');
                
                    % Keep per-step arrays for artifact (if you save it later)
                    DIR.hd_sym_k   = hd_sym_k;
                    DIR.mw_gray_k  = mw_gray_k;
                    DIR.mw_ddra_k  = mw_ddra_k;
                end


                
                % --- Pack row ---
                % ensure percentages
                if (ctrain_gray <= 1), ctrain_gray = 100*ctrain_gray; end
                if (cval_gray   <= 1), cval_gray   = 100*cval_gray;   end
                if (cval_ddra   <= 1), cval_ddra   = 100*cval_ddra;   end
    
                row = pack_row(C, D, alpha_w, pe, ...
                    ctrain_gray, cval_gray, cval_ddra, sizeI_ddra, sizeI_gray, ...
                    Zinfo.rankZ, Zinfo.condZ, ...
                    Tlearn, Tcheck, Tinfer, Tlearn_g, Tvalidate_g, Tinfer_g);

                % Summaries (AUC_k and FV percentiles)
                row.cov_auc_gray = mean(cov_gray_k, 'omitnan');
                row.cov_auc_ddra = mean(cov_ddra_k, 'omitnan');
                if ~isempty(fv_gray), row.fv_gray_med = median(fv_gray(~isinf(fv_gray))); else, row.fv_gray_med = NaN; end
                if ~isempty(fv_ddra), row.fv_ddra_med = median(fv_ddra(~isinf(fv_ddra))); else, row.fv_ddra_med = NaN; end
                
                % Predeclare directional summaries so schema is stable (will be overwritten if computed)
                row.dir_eps_med = NaN;
                row.dir_eps_p90 = NaN;
                row.hout_med    = NaN;
                row.hout_p90    = NaN;

                % Overwrite directional metrics if computed earlier
                if ~isempty(DIR)
                    row.dir_eps_med  = DIR.eps_med;
                    row.dir_eps_p90  = DIR.eps_p90;
                    row.hout_med     = DIR.hout_med;
                    row.hout_p90     = DIR.hout_p90;
                    row.haus_sym_med = DIR.haus_sym_med;
                    row.haus_sym_p90 = DIR.haus_sym_p90;
                    row.mw_gray_mean = DIR.mw_gray_mean;
                    row.mw_ddra_mean = DIR.mw_ddra_mean;
                else
                    % Make sure fields exist even if not enabled (schema stability)
                    row.haus_sym_med = NaN;
                    row.haus_sym_p90 = NaN;
                    row.mw_gray_mean = NaN;
                    row.mw_ddra_mean = NaN;
                end



                % --- Optional: save plotting artifact for this row ---
                if isfield(cfg,'io') && isfield(cfg.io,'save_artifacts') && cfg.io.save_artifacts                
                    % pack only what plotting needs
                    artifact = struct();
                    artifact.sys_gray = configs{idxGray}.sys;     % linearSysDT
                    artifact.sys_ddra = sys_true_dt;              % linearSysDT
                    artifact.VAL      = VAL;  artifact.VAL.R0 = R0;  artifact.VAL.U  = U;
                    artifact.W_eff    = W_used;                    % actually used in inference
                    artifact.M_AB     = M_AB;
                    artifact.meta     = struct('row', row_index, 'dt', sys_true_dt.dt, ...
                     'D', D, 'alpha_w', alpha_w, 'n_m', C.shared.n_m, 'n_s', C.shared.n_s, ...
                     'n_k', C.shared.n_k, 'pe', pe, 'save_tag', cfg.io.save_tag, 'schema', 2);

                    % artifact = struct();
                    % artifact.sys_gray = configs{idxGray}.sys;
                    % artifact.sys_ddra = sys_true_dt;
                    % artifact.VAL      = VAL;
                    % artifact.W_eff    = W_used;
                    % artifact.meta     = struct('row', row_index, 'dt', sys_true_dt.dt, ...
                    %   'D', D, 'alpha_w', alpha_w, 'n_m', C.shared.n_m, 'n_s', C.shared.n_s, ...
                    %   'n_k', C.shared.n_k, 'pe', pe, 'save_tag', getfielddef(cfg.io,'save_tag',''), 'schema', 2);
                    % 
                    % if variant == "standard"
                    %     artifact.M_AB = M_AB;
                    % else
                    %     artifact.AB_meas_center = ABc;
                    %     artifact.AV_oneterm     = AV_one;
                    %     artifact.V_meas         = V_meas;
                    % end

                    artifact.metrics.sizeI_gray_k = wid_gray_k;

                    if exist('sizeI_ddra_k','var'), artifact.metrics.sizeI_ddra_k = sizeI_ddra_k; end

                    try
                        % Per-step volume ratio: Gray vs True (averaged over blocks)
                        nkv = C.shared.n_k_val;
                        ratio_k = zeros(nkv,1); cnt_k = zeros(nkv,1);
                    
                        % Build reach options without cs
                        optReach = getfielddef(C.shared,'options_reach', struct());
                        if isfield(optReach,'cs'), optReach = rmfield(optReach,'cs'); end
                    
                        for b = 1:numel(VAL.x0)
                            params_true = struct('R0', VAL.R0, 'U', U, 'u', VAL.u{b}', 'tFinal', sys_true_dt.dt*(nkv-1));
                            params_gray = struct('R0', configs{idxGray}.params.R0 + VAL.x0{b}, ...
                                                 'u',  VAL.u{b}', 'tFinal', configs{idxGray}.sys.dt*(nkv-1));
                            if configs{idxGray}.sys.nrOfDisturbances>0
                                params_gray.W = W_pred;  % already in disturbance space
                            end
                    
                            Rt = reach(sys_true_dt, params_true, optReach).timePoint.set;
                            Rg = reach(configs{idxGray}.sys, params_gray, optReach).timePoint.set;
                    
                            for ps = 1:nkv
                                Yt = linearMap(Rt{ps}, sys_true_dt.C);  yt = norm(Yt);
                                Yg = linearMap(Rg{ps}, configs{idxGray}.sys.C); yg = norm(Yg);
                                if yt > 0
                                    ratio_k(ps) = ratio_k(ps) + yg/yt;  cnt_k(ps) = cnt_k(ps) + 1;
                                end
                            end
                        end

                        if ~isempty(DIR)
                            artifact.metrics.dir_eps_med = DIR.eps_med;
                            artifact.metrics.dir_eps_p90 = DIR.eps_p90;
                            artifact.metrics.hout_med    = DIR.hout_med;
                            artifact.metrics.hout_p90    = DIR.hout_p90;
                    
                            artifact.metrics.haus_sym_k  = DIR.hd_sym_k;
                            artifact.metrics.mw_gray_k   = DIR.mw_gray_k;
                            artifact.metrics.mw_ddra_k   = DIR.mw_ddra_k;

                            artifact.metrics.ratio_gray_vs_true_k = ratio_k ./ max(1, cnt_k);
                    
                            % Also ensure row is overwritten here (if not already)
                            row.dir_eps_med  = DIR.eps_med;
                            row.dir_eps_p90  = DIR.eps_p90;
                            row.hout_med     = DIR.hout_med;
                            row.hout_p90     = DIR.hout_p90;
                            row.haus_sym_med = DIR.haus_sym_med;
                            row.haus_sym_p90 = DIR.haus_sym_p90;
                            row.mw_gray_mean = DIR.mw_gray_mean;
                            row.mw_ddra_mean = DIR.mw_ddra_mean;
                        end
  
                    catch ME
                        warning(ME.message);
                    end

                    % Make sure vectors exist before streaming
                    if ~exist('wid_ddra_k','var'), wid_ddra_k = []; end
                    if ~exist('wid_gray_k','var'), wid_gray_k = []; end
                    if ~exist('cov_ddra_k','var'), cov_ddra_k = []; end
                    if ~exist('cov_gray_k','var'), cov_gray_k = []; end
                    if ~exist('ratio_k','var') || isempty(ratio_k)
                        % fallback: try pull from artifact, else NaNs
                        ratio_k = [];
                        if exist('artifact','var') && isfield(artifact,'metrics') && isfield(artifact.metrics,'ratio_gray_vs_true_k')
                            ratio_k = artifact.metrics.ratio_gray_vs_true_k;
                        end
                    end
                   
                    sweepio_save_artifact(IO, row_index, artifact);

                    clear configs  

                end

                % Ensure vectors exist
                if ~exist('wid_ddra_k','var'), wid_ddra_k = []; end
                if ~exist('wid_gray_k','var'), wid_gray_k = []; end
                if ~exist('cov_ddra_k','var'), cov_ddra_k = []; end
                if ~exist('cov_gray_k','var'), cov_gray_k = []; end
                
                % Try to pull ratio_k from artifact if it was computed; else leave empty
                ratio_k = [];
                if exist('artifact','var') && isfield(artifact,'metrics') && ...
                        isfield(artifact.metrics,'ratio_gray_vs_true_k')
                    ratio_k = artifact.metrics.ratio_gray_vs_true_k;
                end
                
                sweepio_stream_perstep(IO, row_index, wid_ddra_k, wid_gray_k, ...
                       ratio_k, cov_ddra_k, cov_gray_k);

                
                %% end patch
                row.use_noise = use_noise;
                row.ddra_ridge       = ridgeInfo.used;
                row.ddra_lambda      = ridgeInfo.lambda;
                row.ddra_kappa       = ridgeInfo.kappa;
                row.ddra_ridge_policy= char(ridgeInfo.policy);

                % --- Initialize schema & write/accumulate ----
                IO = sweepio_write_row(IO, row);
                fprintf('[%s] row %d/%d (%.1f%%) | elapsed %s | ERA %s | avg/row %.2fs | tag=%s\n', ...
                    datestr(now,'HH:MM:SS'), IO.rowi, NALL, 100*IO.rowi/NALL, ...
                    char(duration(0,0,toc(t0_all),'Format','hh:mm:ss')), ...
                    char(duration(0,0,(NALL-IO.rowi)*(toc(t0_all)/max(IO.rowi,1)),'Format','hh:mm:ss')), ...
                    toc(t0_all)/max(IO.rowi,1), char(getfielddef(cfg.io,'save_tag','')) );

                clear M_AB Zinfo sizeI_ddra sizeI_gray TS_train TS_val DATASET DATASET_val VAL
              end
            end
          end
        end
      end
    end

    % ------------------------- Save & return ------------------------------
    SUMMARY = sweepio_finalize(IO);
   
end

% ---------- Helpers ----------
function s = zono_gnorm(Z)
% mean 2-norm of generators (a simple, monotone width proxy)
try
    G = generators(Z);
    s = mean(vecnorm(double(G),2,1), 'omitnan');
    if ~isfinite(s), s = NaN; end
catch
    s = NaN;
end
end

function Z2 = zono_scale(Z, fac)
% scale both center and generators by fac
try
    c = double(center(Z)); G = double(generators(Z));
    Z2 = zonotope(c*fac, G*fac);
catch
    Z2 = Z; % safe fallback
end
end


