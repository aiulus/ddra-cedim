function SUMMARY = run_sweeps(cfg, grid)
% RUN_SWEEPS  Batch experiments for DDRA vs Gray under controlled sweeps.
% Produces CSV summary in results_dir and plots in plots_dir (via callers).
    if ~getfielddef(cfg,'allow_parallel',false)
        try
            ps = parallel.Settings; ps.Pool.AutoCreate = false;
            p = gcp('nocreate'); if ~isempty(p), delete(p); end
        catch, end
    end

    [plots_dir, results_dir] = init_io(cfg); 
    csv_path = fullfile(results_dir, 'summary.csv');
    csv_perstep = fullfile(results_dir, 'summary_perstep.csv');
    if exist(csv_perstep,'file'), delete(csv_perstep); end


    % -------- Low-memory / IO toggles (with safe defaults) ---------------
    LM = getfielddef(cfg, 'lowmem', struct());
    LM.gray_check_contain = getfielddef(LM, 'gray_check_contain', true);
    LM.store_ddra_sets    = getfielddef(LM, 'store_ddra_sets', true);
    LM.append_csv         = getfielddef(LM, 'append_csv', false);
    %LM.store_ddra_sets = false; 

    % If streaming to CSV, start clean
    if LM.append_csv && exist(csv_path, 'file')
        delete(csv_path);
    end

    % -------- Sweep axes & base config -----------------------------------
    [axes, baseC] = init_sweep_axes(cfg, grid);
    t0_all = tic;
    NALL = numel(axes.D) * numel(axes.alpha_w) * numel(axes.n_m) * numel(axes.n_s) * numel(axes.n_k) * numel(axes.pe);

    if ~isfield(cfg,'metrics') || ~isfield(cfg.metrics,'tol'), cfg.metrics.tol = 1e-6; end
    % ---- PE policy: explicit for "pe_sweep", minimal otherwise unless forced ----
    is_pe_sweep_tag = contains(lower(string(getfielddef(cfg.io,'save_tag',""))), 'pe_sweep');
    has_forced = any(cellfun(@(p) isfield(p,'force_order') && p.force_order, axes.pe));
    baseC.shared.pe_policy = getfielddef(baseC.shared,'pe_policy', ...
        ternary(is_pe_sweep_tag || has_forced, "explicit", "minimal"));

    function out = ternary(cond, a, b), if cond, out=a; else, out=b; end, end


    % Precompute for in-memory mode only
    if ~LM.append_csv
        Ntot = numel(axes.D) * numel(axes.alpha_w) * numel(axes.n_m) * ...
               numel(axes.n_s) * numel(axes.n_k) * numel(axes.pe);
        cells = [];  % allocate after first row when the header is known
    end
    hdr  = {};
    rowi = 0;
    % Force explicit PE whenever there are multiple PE orders
    % on the grid
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
                row_index = rowi + 1;
                
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

                %% New
                pe_eff = pe_normalize(pe, U, sys_cora, C.shared.n_k);

                % ================= DDRA =================
                t0 = tic;
                % DATASET.x0_blocks{b}, DATASET.u_blocks{b} with block length n_k
                [Xminus, Uminus, Xplus, W, Zinfo, DATASET] = ddra_generate_data(C, sys_ddra, sys_cora.dt, R0, U, pe_eff, sys_cora);
                Tlearn = toc(t0);
    
                % Build TRAIN and VAL suites deterministically from generator
                TS_train = testSuite_fromDDRA(sys_cora, R0, DATASET, C.shared.n_k, C.shared.n_m, C.shared.n_s);

                if string(C.shared.pe_policy) == "explicit"
                    % Gate on the requested order when sweeping PE explicitly
                    Lgate = getfielddef(pe,'order_gate', getfielddef(pe,'order', pe_eff.order));
                else
                    % Gate on the effective minimal sufficient order
                    Lgate = pe_eff.order;
                end

                % DEBUG STATEMENT
                fprintf('L requested=%g  effective=%g  policy=%s\n', ...
                    getfielddef(pe,'order',NaN), Lgate, string(C.shared.pe_policy));
               
                okPE = check_PE_order(TS_train, Lgate);

                
                if ~okPE
                    % mark and skip this row to keep the grid aligned
                    emit_skip_row("PE_not_satisfied");
                    clear TS_train DATASET
                    continue; % skip to next sweep point
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
                t1 = tic;
                try
                    % New signature with ridge control & W_eff (possibly inflated for ridge)
                    [M_AB, ridgeInfo, W_eff] = ddra_learn_Mab(Xminus, Uminus, Xplus, W, Zinfo, sys_ddra, C);
                catch
                    % Old signature fallback (keeps older ddra_learn_Mab working)
                    M_AB = ddra_learn_Mab(Xminus, Uminus, Xplus, W, Zinfo, sys_ddra);
                    ridgeInfo = struct('used', false, 'skipped', false, ...
                                       'lambda', 0, 'kappa', NaN, 'policy', "none");
                    W_eff = W;  % no inflation in legacy path
                end
                Tcheck = toc(t1);

                %% Debug sequence
                boxM = intervalMatrix(M_AB).int;
                ab_span = boxM.sup - boxM.inf;
                row.ab_hw_mean = mean(ab_span(:))/2;     % mean half-width of all AB entries
                row.ab_hw_max  = max(ab_span(:))/2;      % max half-width
                fprintf("ab_hw_mean: %.8f, ab_hw_max: %.8f\n", row.ab_hw_mean, row.ab_hw_max)
                try
                    row.ab_ngen = size(M_AB.G,3);        % #gens in matZonotope (if available)
                catch
                    row.ab_ngen = NaN;
                end


                % Fairness: if ridge would inflate MAB and we're not mirroring it into Gray, skip row for both
                if isfield(ridgeInfo,'used') && ridgeInfo.used && ...
                   isfield(ridgeInfo,'policy') && string(ridgeInfo.policy)=="inflate_MAB"
                    emit_skip_row("ridge_inflate_MAB_unmirrored", ...
                        'ddra_ridge',true, ...
                        'ddra_lambda',ridgeInfo.lambda, ...
                        'ddra_kappa',ridgeInfo.kappa, ...
                        'ddra_ridge_policy',char(ridgeInfo.policy));
                    clear Xminus Xplus TS_train TS_val DATASET DATASET_val
                    continue;

                end

                
                % If rank-deficient and allow_ridge = false, skip this sweep point
                if isfield(ridgeInfo,'skipped') && ridgeInfo.skipped
                    emit_skip_row("skip", ...
                        'ddra_ridge',false, ...
                        'ddra_lambda',0, ...
                        'ddra_kappa',NaN, ...
                        'ddra_ridge_policy',"skip");
                    clear Xminus Xplus TS_train TS_val DATASET DATASET_val
                    continue;
                end
                
                clear Xminus Xplus  % free data blocks early

    
                % ---- unified disturbance policy ----
                % W_eff: from ddra_learn_Mab (may include ridge inflation)
                % ---- unified disturbance policy (unchanged) ----
                if resolve_use_noise(C.shared)
                    W_used = W_eff;    % state-dim nx
                else
                    W_used = zonotope(zeros(size(center(W_eff),1),1));
                end
                

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
                W_pred = build_W_pred(configs{idxGray}.sys, U, W_used);
                if ~isfield(configs{idxGray},'params') || isempty(configs{idxGray}.params)
                    configs{idxGray}.params = struct();
                end
                if ~isempty(W_pred)
                    configs{idxGray}.params.W = W_pred;
                end

                %fprintf('Reporting from line 257@run_sweeps.m. \n ||W_eff||_1=%.3g  ||E*Wd||_1=%.3g\n', ...
                %   sum(abs(generators(W_used)),'all'), ...
                %   sum(abs(sys_cora.E*generators(coerceWToSys(sys_cora, W_used))),'all'));

                % ---- Debug prints ----
                valE = 0;
                if isprop(sys_cora,'nrOfDisturbances') && sys_cora.nrOfDisturbances > 0
                    try
                        Wc = coerceWToSys(sys_cora, W_used);
                        if isa(Wc,'zonotope')
                            valE = sum(abs(sys_cora.E * generators(Wc)), 'all');
                        end
                    catch
                        % leave valE = 0
                    end
                end
                try
                    valW = sum(abs(generators(W_used)),'all');
                catch
                    valW = 0;
                end
                %fprintf('Reporting from run_sweeps: ||W_eff||_1=%.3g  ||E*Wd||_1=%.3g\n', valW, valE);


                % Build VAL record (x0,u) from TS_val and pass to both sides
                % VAL = VAL_from_TS(TS_val, DATASET_val);
                VAL = pack_VAL_from_TS(TS_val);

                fprintf('VAL(flat) m*s=%d  nk=%d  hash=%s\n', ...
                     numel(VAL.y), size(VAL.y{1},1), val_hash(VAL));
    
                % DEBUG SEQUENCE - START
                % Quick internal parity check (non-fatal; hits breakpoint if mismatch)
                % --- Debug: VAL parity check (non-fatal)
                try
                    Vcanon = pack_VAL_from_TS(TS_val);
                    okVAL  = strcmp(val_hash(VAL), val_hash(Vcanon));
                    if ~okVAL
                        warning('VAL != canonical(TS_val) — check pack_VAL_from_TS / TS construction.');
                    end
                catch ME
                    % If val_hash is missing, just warn (don't stop experiments)
                    warning(ME.identifier, '%s',  ME.message);
                end

                % DEBUG SEQUENCE - END

                pmode = lower(string(getfielddef(cfg.io,'plot_mode','offline')));
                % === Unified artifact save (after configs, VAL, M_AB, W_eff exist) ===
                row_index = rowi + 1;
                artdir = fullfile(results_dir, 'artifacts');
                if ~exist(artdir,'dir'), mkdir(artdir); end
                
                % normalize “true” to linearSysDT for consistent on-disk schema
                sys_true_dt = normalize_to_linearSysDT(sys_ddra, sys_cora.dt);
                
                artifact = struct();
                artifact.sys_gray = configs{idxGray}.sys;     % RCSI/Gray (linearSysDT)
                artifact.sys_ddra = sys_true_dt;        % “true” model (linearSysDT)
                artifact.VAL      = VAL;                % exact (x0,u,y) used on VAL
                artifact.W_eff    = W_eff;              % effective W used by DDRA
                artifact.M_AB     = M_AB;               % matrix zonotope for DDRA
                artifact.meta     = struct( ...
                    'row', row_index, 'dt', sys_true_dt.dt, ...
                    'D', D, 'alpha_w', alpha_w, ...
                    'n_m', C.shared.n_m, 'n_s', C.shared.n_s, 'n_k', C.shared.n_k, ...
                    'pe', pe, 'save_tag', cfg.io.save_tag, 'schema', 2 );
                
                % include R0/U so the plotter can exactly reproduce reach
                artifact.VAL.R0 = R0; 
                artifact.VAL.U  = U;

                want_online = any(pmode == ["online","both"]) && getfielddef(cfg.io,'make_reach_plot',true);
                if want_online && ismember(row_index, getfielddef(cfg.io,'plot_rows',[1]))
                    dims   = getfielddef(cfg.io,'plot_dims',[1 2]);
                
                    sys_true_dt = normalize_to_linearSysDT(sys_ddra, sys_cora.dt);
                
                    W_plot = pick_W_for_plot(W_used, use_noise);

                    reach_dir = fullfile(plots_dir,'reach');
                    if ~exist(reach_dir,'dir'), mkdir(reach_dir); end
                    
                    savename_png = fullfile(reach_dir, sprintf('row_%04d_y%d%d.png', row_index, dims(1), dims(2)));
                    % (If save_plot also writes PDF, it can derive .pdf from this base. If not,
                    % it will at least save a valid .png and not warn.)

                    plot_reach_all_onepanel( ...
                        sys_true_dt, ...
                        configs{idxGray}.sys, ...
                        VAL, ...
                        'MAB', M_AB, ...
                        'W', W_plot, ...
                        'Dims', dims, ...
                        'ShowSamples', false, ...
                        'Save', savename_png, ...
                        'ReachOptions', getfielddef(C.shared,'options_reach', ...
                                            struct('zonotopeOrder',100,'reductionTechnique','girard')));


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

                %% Debug prints
                % 1) How many disturbance channels did the Gray model end up with?
                %fprintf('Gray: nrOfDisturbances = %d, size(E) = [%d %d]\n', ...
                %    configs{idxGray}.sys.nrOfDisturbances, size(configs{idxGray}.sys.E));
                
                % 2) What is the dimension of W_for_gray you pass in?
                %if isa(W_for_gray,'zonotope')
                %    fprintf('W_for_gray dim = %d, #gens = %d\n', size(center(W_for_gray),1), size(generators(W_for_gray),2));
                %else
                %    disp('W_for_gray is empty or not a zonotope');
                %end
                
                % 3) How big is W_true vs W_for_gray?
                %if isa(W_eff,'zonotope')
                %    r_true = sum(abs(generators(W_eff)), 'all');
                %    fprintf('||W_eff.generators||_1 = %.3g\n', r_true);
                %end
                %if isa(W_for_gray,'zonotope')
                %    r_gray = sum(abs(generators(W_for_gray)), 'all');
                %    fprintf('||W_for_gray.generators||_1 = %.3g\n', r_gray);
                %end
                %% End - Debug prints
                
                [ctrain_gray, cval_gray, Tvalidate_g, VAL] = gray_containment( ...
                                configs, sys_cora, R0_cora, U, C, pe_eff, ...
                                'externalTS_train', TS_train, ...
                                'externalTS_val',   TS_val,   ...
                                'check_contain',    true,     ...
                                ps_args{:});

    
                % ---- DDRA inference (stored or streaming) on SAME validation points
                % ---- DDRA inference (stored or streaming) on SAME validation points
                if LM.store_ddra_sets
                    t2 = tic;
                    [Xsets_ddra, ~] = ddra_infer(sys_ddra, R0, U, W_used, M_AB, C, VAL);
                    Tinfer = toc(t2);
                
                    % Map sets to output space correctly
                    Ysets_ddra = cellfun(@(X) linearMap(X, sys_true_dt.C), Xsets_ddra, 'uni', 0);
                
                    % Per-step interval widths (mean over outputs), then canonicalize to "mean"
                    wid_ddra_k_raw = zeros(numel(Ysets_ddra),1);
                    for k = 1:numel(Ysets_ddra)
                        Ik = interval(Ysets_ddra{k});
                        try, lo = infimum(Ik); hi = supremum(Ik);
                        catch, lo = Ik.inf;    hi = Ik.sup;
                        end
                        wvec = max(hi - lo, 0);
                        wid_ddra_k_raw(k) = mean(double(wvec), 'omitnan');
                    end
                    [ sizeI_ddra, wid_ddra_k ] = normalize_widths( ...
                        wid_ddra_k_raw, size(sys_true_dt.C,1), "mean");
                    sizeI_ddra_k = wid_ddra_k;   % keep per-step array for artifacts/CSV
                
                    % Containment must be in output-space too
                    cval_ddra = containsY_on_VAL(Ysets_ddra, VAL, cfg.metrics.tol);
                
                    clear Xsets_ddra
                else
                    t2 = tic;
                    [sizeI_ddra, cval_ddra, wid_ddra_k] = ddra_infer_size_streaming( ...
                        sys_ddra, R0, [], W_used, M_AB, C, VAL);
                    Tinfer = toc(t2);
                    [ sizeI_ddra, wid_ddra_k ] = normalize_widths( ...
                        wid_ddra_k, size(sys_true_dt.C,1), getfielddef(cfg.metrics,'width_agg',"mean"));
                    sizeI_ddra_k = wid_ddra_k;
                end


                % ---- Gray size metric (interval proxy), consistent with DDRA
                t4 = tic;
                [sizeI_gray, wid_gray_k] = gray_infer_size_on_VAL(configs{idxGray}.sys, TS_val, C, ...
                    configs{idxGray}.params, 'overrideW', W_pred, 'width_agg', 'sum');
                Tinfer_g   = toc(t4);
                [ sizeI_gray, wid_gray_k ] = normalize_widths(wid_gray_k, size(configs{idxGray}.sys.C,1), getfielddef(cfg.metrics,'width_agg',"mean"));

                cov_ddra_k = []; cov_gray_k = []; fv_ddra = []; fv_gray = [];

                if getfielddef(cfg.metrics,'enhanced', true)
                    % GRAY per-step metrics from VAL only
                    try
                        [~, cov_gray_k, fv_gray, ~] = perstep_gray_metrics(configs{idxGray}.sys, VAL, W_pred, struct('zonotopeOrder',60));
                    catch
                        cov_gray_k = []; fv_gray = []; 
                    end
                
                    % DDRA per-step metrics from VAL only (uses M_AB; fast)
                    try
                        [~, cov_ddra_k, fv_ddra, ~] = perstep_ddra_metrics(sys_true_dt, VAL, M_AB, W_used, 60);
                    catch
                        cov_ddra_k = []; fv_ddra = [];
                    end
                end

                
                % --- Pack row ---
                % ensure percentages
                if (ctrain_gray <= 1), ctrain_gray = 100*ctrain_gray; end
                if (cval_gray   <= 1), cval_gray   = 100*cval_gray;   end
                if (cval_ddra   <= 1), cval_ddra   = 100*cval_ddra;   end

                rowi = rowi + 1;      
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


                % --- Optional: save plotting artifact for this row ---
                if isfield(cfg,'io') && isfield(cfg.io,'save_artifacts') && cfg.io.save_artifacts
                    % ensure folder
                    artdir = fullfile(results_dir, 'artifacts');
                    if ~exist(artdir,'dir'), mkdir(artdir); end
                
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
                    
                            for k = 1:nkv
                                Yt = linearMap(Rt{k}, sys_true_dt.C);  yt = norm(Yt);
                                Yg = linearMap(Rg{k}, configs{idxGray}.sys.C); yg = norm(Yg);
                                if yt > 0
                                    ratio_k(k) = ratio_k(k) + yg/yt;  cnt_k(k) = cnt_k(k) + 1;
                                end
                            end
                        end

                        % --- Directional support / Hausdorff-outer (shape-aware) ---
                        dir_cfg = getfielddef(cfg,'metrics',struct());
                        dir_cfg = getfielddef(dir_cfg,'directional', struct('Nd',64,'seed',12345));
                        Nd   = getfielddef(dir_cfg,'Nd',64);
                        dseed= getfielddef(dir_cfg,'seed',12345);
                        
                        ny = size(sys_true_dt.C,1);
                        Ddirs = sample_unit_dirs(ny, Nd, dseed); % (ny x Nd), unit directions
                        
                        dir_eps_all = [];   % accumulate eps over blocks x steps x dirs
                        dir_del_all = [];   % accumulate delta over blocks x steps x dirs
                        Hout_k = zeros(nkv,1); cnt_kH = zeros(nkv,1);
                        
                        for b = 1:numel(VAL.x0)
                            Ublk = VAL.u{b};
                            for k = 1:nkv
                                uk = get_uk(Ublk, k, size(sys_true_dt.D,2)); % (m x 1)
                        
                                % True Y_k and Gray Y_k as zonotopes (outer-approx if needed)
                                Yt_k = toZono(linearMap(Rt{k}, sys_true_dt.C));
                                Yg_k = toZono(linearMap(Rg{k}, configs{idxGray}.sys.C));
                        
                                % Support values along all directions (vectorized)
                                st = support_zono_vec(Ddirs, center(Yt_k), generators(Yt_k)) + Ddirs'*(sys_true_dt.D*uk);
                                sg = support_zono_vec(Ddirs, center(Yg_k), generators(Yg_k)) + Ddirs'*(configs{idxGray}.sys.D*uk);
                        
                                eps_kd = sg ./ max(st, eps);
                                del_kd = max(sg - st, 0);
                        
                                dir_eps_all = [dir_eps_all; eps_kd(:)]; 
                                dir_del_all = [dir_del_all; del_kd(:)]; 

                                Hout_k(k) = Hout_k(k) + max(del_kd);
                                cnt_kH(k) = cnt_kH(k) + 1;
                            end
                        end
                        
                        artifact.metrics.dir_eps_med = median(dir_eps_all, 'omitnan');
                        artifact.metrics.dir_eps_p90 = prctile(dir_eps_all, 90);
                        artifact.metrics.hout_med    = median(dir_del_all, 'omitnan');
                        artifact.metrics.hout_p90    = prctile(dir_del_all, 90);
                        artifact.metrics.hout_k      = Hout_k ./ max(1,cnt_kH);

                        % Overwrite row-level directional summaries now that we have them
                        row.dir_eps_med = artifact.metrics.dir_eps_med;
                        row.dir_eps_p90 = artifact.metrics.dir_eps_p90;
                        row.hout_med    = artifact.metrics.hout_med;
                        row.hout_p90    = artifact.metrics.hout_p90;


                        artifact.metrics.ratio_gray_vs_true_k = ratio_k ./ max(1, cnt_k);
                        
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
                   
                    save(fullfile(artdir, sprintf('row_%04d.mat', row_index)), '-struct','artifact','-v7.3');

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
                
                stream_perstep(csv_perstep, row_index, wid_ddra_k, wid_gray_k, ...
                               ratio_k, cov_ddra_k, cov_gray_k);


                
                %% end patch
                row.use_noise = use_noise;
                row.ddra_ridge       = ridgeInfo.used;
                row.ddra_lambda      = ridgeInfo.lambda;
                row.ddra_kappa       = ridgeInfo.kappa;
                row.ddra_ridge_policy= char(ridgeInfo.policy);

    
                % --- Initialize schema & write/accumulate ----
                write_row(row, csv_path, LM);
                fprintf('[%s] row %d/%d (%.1f%%) | elapsed %s | ERA %s | avg/row %.2fs | tag=%s\n', ...
                        datestr(now,'HH:MM:SS'), rowi, NALL, 100*rowi/NALL, ...
                        char(duration(0,0,toc(t0_all),'Format','hh:mm:ss')), ...
                        char(duration(0,0,(NALL-rowi)*(toc(t0_all)/max(rowi,1)),'Format','hh:mm:ss')), ...
                        toc(t0_all)/max(rowi,1), char(getfielddef(cfg.io,'save_tag','')));


                clear M_AB Zinfo sizeI_ddra sizeI_gray TS_train TS_val DATASET DATASET_val VAL
              end
            end
          end
        end
      end
    end

    % ------------------------- Save & return ------------------------------
    if LM.append_csv
        fprintf('Sweeps done. Rows: %d\nCSV -> %s\n', rowi, csv_path);
        SUMMARY = readtable(csv_path);
    else
        writecell([hdr; cells(1:rowi,:)], csv_path);
        fprintf('Sweeps done. Rows: %d\nCSV -> %s\n', rowi, csv_path);
        SUMMARY = cell2table(cells(1:rowi,:), 'VariableNames', hdr);
    end

    % ---------- nested helpers (drop-in) ------------------------------------

    function write_row_init_if_needed(row, csv_path, LM)
        % Initializes header and in-memory cells (if needed) using caller scope.
        % Relies on outer-scope variables: hdr, cells, Ntot, rowi
        persistent initialized
        if isempty(initialized) || rowi == 1
            hdr = fieldnames(orderfields(row))';
            if LM.append_csv
                fid = fopen(csv_path, 'w'); fprintf(fid, '%s\n', strjoin(hdr, ',')); fclose(fid);
            else
                cells = cell(Ntot, numel(hdr));
            end
            initialized = true;
        end
    end

    function write_row(row, csv_path, LM)
        % Writes one row to CSV or in-memory cells (same logic you already use)
        % Relies on outer-scope: hdr, cells, rowi
        write_row_init_if_needed(row, csv_path, LM);
        if LM.append_csv
            append_row_csv(csv_path, hdr, row);
        else
            for j = 1:numel(hdr)
                v = row.(hdr{j}); if isstring(v), v = char(v); end
                cells{rowi, j} = v;
            end
        end
    end
    
    function emit_skip_row(reason, varargin)
        % Common skip pathway to avoid copy-paste
        % varargin lets you pass name/value that end up on row fields
        % Relies on outer-scope vars already in run_sweeps: C, D, alpha_w, pe,
        % Zinfo, Tlearn, Tcheck, use_noise, csv_path, LM, rowi
        if exist('Tlearn','var'),        Tlearn_val = Tlearn;        else, Tlearn_val = NaN; end
        if exist('Tcheck','var'),        Tcheck_val = Tcheck;        else, Tcheck_val = NaN; end
        if exist('Tinfer','var'),        Tinfer_val = Tinfer;        else, Tinfer_val = NaN; end
        if exist('Tlearn_g','var'),      Tlearn_g_val = Tlearn_g;    else, Tlearn_g_val = NaN; end
        if exist('Tvalidate_g','var'),   Tvalidate_g_val = Tvalidate_g; else, Tvalidate_g_val = NaN; end
        if exist('Tinfer_g','var'),      Tinfer_g_val = Tinfer_g;    else, Tinfer_g_val = NaN; end
    
        rowi = rowi + 1;
        row = pack_row(C, D, alpha_w, pe, ...
            NaN, NaN, NaN, NaN, NaN, ...
            Zinfo.rankZ, Zinfo.condZ, ...
            Tlearn_val, Tcheck_val, Tinfer_val, Tlearn_g_val, Tvalidate_g_val, Tinfer_g_val);
        row.skipped     = true;
        row.skip_reason = string(reason);
        row.use_noise   = use_noise;

        % Ensure schema stability for columns that appear in non-skip rows
        row.cov_auc_gray = NaN;
        row.cov_auc_ddra = NaN;
        row.fv_gray_med  = NaN;
        row.fv_ddra_med  = NaN;
        
        row.dir_eps_med  = NaN;
        row.dir_eps_p90  = NaN;
        row.hout_med     = NaN;
        row.hout_p90     = NaN;
        
        % Ridge defaults (match non-skip rows)
        row.ddra_ridge        = false;
        row.ddra_lambda       = 0;
        row.ddra_kappa        = NaN;
        row.ddra_ridge_policy = "none";

    
        % attach any extra fields (e.g., ridge info)
        for k = 1:2:numel(varargin)
            row.(varargin{k}) = varargin{k+1};
        end
    
        write_row(row, csv_path, LM);
    end
    
    function idxGray = pick_gray_config(configs, C)
        % Same logic you use in multiple places, consolidated
        want = "graySeq";
        try
            if isfield(C,'gray') && isfield(C.gray,'methodsGray') && ~isempty(C.gray.methodsGray)
                want = string(C.gray.methodsGray(1));
            end
        catch
        end
        idxGray = find(cellfun(@(c) isfield(c,'name') && want==string(c.name), configs), 1, 'first');
        if isempty(idxGray), idxGray = min(2, numel(configs)); end
    end
    
    function Wfg = build_W_for_gray(sys_cora, ~, W_used)
        % Map state-noise W_used -> disturbance space so E*Wfg \supseteq W_used.
        if isprop(sys_cora,'nrOfDisturbances') && sys_cora.nrOfDisturbances > 0
            Wfg = normalizeWForGray(sys_cora, W_used);
        else
            Wfg = [];
        end
    end
    
    function Wpred = build_W_pred(gray_sys, ~, W_used)
        if isprop(gray_sys,'nrOfDisturbances') && gray_sys.nrOfDisturbances > 0
            Wpred = normalizeWForGray(gray_sys, W_used);
        else
            Wpred = [];
        end
    end

    function ensure_perstep_header(csv_perstep)
        if ~exist(csv_perstep,'file') || dir(csv_perstep).bytes==0
            fid_h = fopen(csv_perstep,'w');
            fprintf(fid_h,'row,k,wid_ddra,wid_gray,ratio_gray_true,cov_ddra,cov_gray\n');
            fclose(fid_h);
        end
    end
    
    function stream_perstep(csv_perstep, row_index, wid_ddra_k, wid_gray_k, ratio_k, cov_ddra_k, cov_gray_k)
        ensure_perstep_header(csv_perstep);
        nkv = max([numel(wid_gray_k), numel(wid_ddra_k), numel(ratio_k), numel(cov_ddra_k), numel(cov_gray_k)]);
        if isempty(wid_ddra_k), wid_ddra_k = nan(nkv,1); end
        if isempty(wid_gray_k), wid_gray_k = nan(nkv,1); end
        if isempty(ratio_k),    ratio_k    = nan(nkv,1); end
        if isempty(cov_ddra_k), cov_ddra_k = nan(nkv,1); end
        if isempty(cov_gray_k), cov_gray_k = nan(nkv,1); end
    
        fid_ps = fopen(csv_perstep,'a');
        for kk = 1:nkv
            fprintf(fid_ps, '%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g\n', ...
                row_index, kk, wid_ddra_k(min(kk,end)), wid_gray_k(min(kk,end)), ...
                ratio_k(min(kk,end)), cov_ddra_k(min(kk,end)), cov_gray_k(min(kk,end)));
        end
        fclose(fid_ps);
    end

    function Wplot = pick_W_for_plot(W_used, use_noise)
        if use_noise, Wplot = W_used;
        else,        Wplot = zonotope(zeros(size(center(W_used),1),1));
        end
    end
    
    function logf(varargin)
        % toggle with cfg.io.verbose=true/false if you want
        try, vb = logical(getfielddef(cfg,'io',struct()).verbose);
        catch, vb = true;
        end
        if vb, fprintf(varargin{:}); end
    end

    function D = sample_unit_dirs(p, Nd, seed)
        if nargin<3 || isempty(seed), seed = 12345; end
        rng(seed,'twister');
        X = randn(p, Nd);
        n = sqrt(sum(X.^2,1)) + eps;
        D = X ./ n;
    end
    
    function s = support_zono_vec(Ddirs, c, G)
    % Ddirs: (p x Nd), c: (p x 1), G: (p x g)
        if isempty(G), s = Ddirs' * c; return; end
        s = Ddirs' * c + sum(abs(Ddirs' * G), 2);
    end
    
    function uk = get_uk(Ublk, k, m)
        if m==0, uk = zeros(0,1); return; end
        Ublk = squeeze(Ublk);
        if size(Ublk,1) == m
            uk = Ublk(:, min(k, size(Ublk,2)));
        elseif size(Ublk,2) == m
            uk = Ublk(min(k, size(Ublk,1)), :)';
        else
            error('get_uk: U has incompatible size for m=%d.', m);
        end
    end

    function Z = toZono(S)
    % Convert various CORA set types to a (conservative) zonotope.
        if isa(S,'zonotope')
            Z = S;
        elseif isa(S,'polyZonotope') || isa(S,'conZonotope')
            Z = zonotope(S);   % CORA outer-approx
        else
            try
                Z = zonotope(S);
            catch
                error('toZono: unsupported set type %s', class(S));
            end
        end
    end


end

function s = stable_seed(D, alpha_w, pe, dyn, typ, tag)
    o     = double(getfielddef(pe,'order',0));
    smode = double(sum(char(string(getfielddef(pe,'mode','')))));
    dyn   = char(string(dyn));
    typ   = char(string(typ));
    tag   = char(string(tag));
    s64   = uint64(1469598103934665603);
    K     = uint64([D, round(alpha_w*1e6), o, smode, sum(double(dyn)), sum(double(typ)), sum(double(tag))]);
    for v = K
        s64 = bitxor(s64, uint64(v));
        s64 = s64 * 1099511628211;
    end
    s = uint32(mod(s64, uint64(2^31-1))); if s==0, s=1; end
end


function h = val_hash(v)
%VAL_HASH Deterministic hash of structs/cells/arrays for parity checks.
% Returns lowercase hex SHA-256 string.

    md = java.security.MessageDigest.getInstance('SHA-256');
    feed(md, v);
    h = lower(reshape(dec2hex(typecast(md.digest(),'uint8'))',1,[]));
end

function feed(md, v)
    if isnumeric(v) || islogical(v)
        % shape + class + raw bytes
        updateStr(md, class(v));
        updateIntVec(md, int64(size(v)));
        if ~isempty(v)
            b = typecast(v(:), 'uint8');  % column-major
            md.update(b);
        end
    elseif ischar(v)
        updateStr(md, 'char'); md.update(uint8(v(:)));
    elseif isstring(v)
        updateStr(md, 'string');
        for s = v(:).'
            feed(md, char(s));
        end
    elseif iscell(v)
        updateStr(md, 'cell'); updateIntVec(md, int64(size(v)));
        for i = 1:numel(v), feed(md, v{i}); end
    elseif isstruct(v)
        updateStr(md, 'struct'); updateIntVec(md, int64(size(v)));
        f = sort(fieldnames(v));
        for k = 1:numel(f)
            fn = f{k}; md.update(uint8(fn));
            for i = 1:numel(v)
                feed(md, v(i).(fn));
            end
        end
    elseif isa(v,'zonotope')
        % If you ever hash sets: reduce to center/generators
        try
            feed(md, center(v)); feed(md, generators(v));
        catch
            updateStr(md, class(v));
        end
    else
        % Fallback to class+num2str (last resort)
        updateStr(md, ['<', class(v), '>']);
        try
            feed(md, num2str(v));
        catch
        end
    end
end

function updateStr(md, s), md.update(uint8(s)); end
function updateIntVec(md, v), md.update(typecast(v(:),'uint8')); end
