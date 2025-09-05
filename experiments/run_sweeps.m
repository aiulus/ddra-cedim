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

    % -------- Low-memory / IO toggles (with safe defaults) ---------------
    LM = getfielddef(cfg, 'lowmem', struct());
    LM.gray_check_contain = getfielddef(LM, 'gray_check_contain', true);
    LM.store_ddra_sets    = getfielddef(LM, 'store_ddra_sets', true);
    LM.append_csv         = getfielddef(LM, 'append_csv', false);

    % If streaming to CSV, start clean
    if LM.append_csv && exist(csv_path, 'file')
        delete(csv_path);
    end

    % -------- Sweep axes & base config -----------------------------------
    [axes, baseC] = init_sweep_axes(cfg, grid);
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

    % ------------------------ Main sweeps --------------------------------
    for D = axes.D
      for alpha_w = axes.alpha_w
        for n_m = axes.n_m
          for n_s = axes.n_s
            for n_k = axes.n_k
              for ip = 1:numel(axes.pe)
                pe = axes.pe{ip};
                % --- per-row deterministic seed (distinct datasets per sweep row) ---
                row_index     = rowi + 1;             % the row we're about to produce
                row_seed      = 10000 + row_index;    % any fixed offset is fine
                C.shared.seed = row_seed;
                rng(row_seed,'twister');
            
                % Make nominal PE input reproducible and L-dependent without freezing RNG
                pe.deterministic = getfielddef(pe,'deterministic', true);
                pe.seed_base = uint32(mod(100000 + 131*row_index + 977*ip + 1009*n_k + 7*n_s + 3*n_m, 2^31-1));

    
                % --- instantiate config for this run ---
                C = baseC;
                C.shared.n_m = n_m;  C.shared.n_s = n_s;  C.shared.n_k = n_k;
                C.shared.n_m_val = max(2, min(n_m, getfielddef(baseC.shared, 'n_m_val', 2)));
                C.shared.n_s_val = n_s;  C.shared.n_k_val = n_k;
                C.ddra.alpha_w   = alpha_w;
                % --- per-row deterministic seed (distinct datasets per sweep row) ---
                row_index     = rowi + 1;             % the row we're about to produce
                row_seed      = 10000 + row_index;    % any fixed offset is fine
                C.shared.seed = row_seed;
                rng(row_seed,'twister');
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
                pe_eff = finalize_pe_order(pe, sys_cora, C);   % pass C, not cfg

                % ================= DDRA =================
                t0 = tic;
                % NOTE: ddra_generate_data must return DATASET as 6th out:
                % DATASET.x0_blocks{b}, DATASET.u_blocks{b} with block length n_k
                [Xminus, Uminus, Xplus, W, Zinfo, DATASET] = ddra_generate_data(C, sys_ddra, sys_cora.dt, R0, U, pe_eff);
                Tlearn = toc(t0);
    
                % Build TRAIN and VAL suites deterministically from generator
                TS_train = testSuite_fromDDRA(sys_cora, R0, DATASET, C.shared.n_k);

                if string(C.shared.pe_policy) == "explicit"
                    % Gate on the requested order when sweeping PE explicitly
                    Lgate = getfielddef(pe,'order_gate', getfielddef(pe,'order', pe_eff.order));
                else
                    % Gate on the effective minimal sufficient order
                    Lgate = pe_eff.order;
                end
                okPE = check_PE_order(TS_train, Lgate);

                if ~okPE
                    % mark and skip this row to keep the grid aligned
                    rowi = rowi + 1;
                    row = pack_row(C, D, alpha_w, pe, ...
                        NaN, NaN, NaN, NaN, NaN, ...
                        Zinfo.rankZ, Zinfo.condZ, ...
                        Tlearn, NaN, NaN, NaN, NaN, NaN);
                    row.skipped      = true;
                    row.skip_reason  = "PE_not_satisfied";
                    row.use_noise    = resolve_use_noise(C.shared);
                
                    % --- write row ---
                    if rowi == 1
                        hdr = fieldnames(orderfields(row))';
                        if LM.append_csv
                            fid = fopen(csv_path, 'w'); fprintf(fid, '%s\n', strjoin(hdr, ',')); fclose(fid);
                        else
                            cells = cell(Ntot, numel(hdr));
                        end
                    end
                    if LM.append_csv
                        append_row_csv(csv_path, hdr, row);
                    else
                        for j = 1:numel(hdr)
                            v = row.(hdr{j}); if isstring(v), v = char(v); end
                            cells{rowi, j} = v;
                        end
                    end
                
                    clear TS_train DATASET  
                    continue;               % skip to next sweep point
                end
    
                C_val = C;
                C_val.shared.n_m = C.shared.n_m_val;
                C_val.shared.n_s = C.shared.n_s_val;
                C_val.shared.n_k = C.shared.n_k_val;
                rng(row_seed+1,'twister');
                [~,~,~,~,~, DATASET_val] = ddra_generate_data(C_val, sys_ddra, sys_cora.dt, R0, U, pe_eff);
                TS_val   = testSuite_fromDDRA(sys_cora, R0, DATASET_val, C_val.shared.n_k);
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
                
                % If rank-deficient and allow_ridge = false, skip this sweep point
                if isfield(ridgeInfo,'skipped') && ridgeInfo.skipped
                    rowi = rowi + 1;
                    row = pack_row(C, D, alpha_w, pe, ...
                        NaN, NaN, NaN, NaN, NaN, ...
                        Zinfo.rankZ, Zinfo.condZ, ...
                        Tlearn, Tcheck, NaN, NaN, NaN, NaN);
                    row.skipped        = true;
                    row.ddra_ridge     = false;
                    row.ddra_lambda    = 0;
                    row.ddra_kappa     = NaN;
                    row.ddra_ridge_policy = "skip";
                    row.use_noise      = use_noise;
                
                    % --- write row ---
                    if rowi == 1
                        hdr = fieldnames(orderfields(row))';
                        if LM.append_csv
                            fid = fopen(csv_path, 'w'); fprintf(fid, '%s\n', strjoin(hdr, ',')); fclose(fid);
                        else
                            cells = cell(Ntot, numel(hdr));
                        end
                    end
                    if LM.append_csv
                        append_row_csv(csv_path, hdr, row);
                    else
                        for j = 1:numel(hdr)
                            v = row.(hdr{j}); if isstring(v), v = char(v); end
                            cells{rowi, j} = v;
                        end
                    end
                
                    % Skip inference (and optionally Gray) to keep grid aligned
                    clear Xminus Uminus Xplus TS_train TS_val DATASET DATASET_val
                    continue;
                end
                
                clear Xminus Uminus Xplus  % free data blocks early

    
                % ---- unified disturbance policy ----
                % W_eff: from ddra_learn_Mab (may include ridge inflation)
                if resolve_use_noise(C.shared)
                    W_used   = W_eff;           % both branches use effective W
                else
                    W_used   = zonotope(zeros(size(center(W_eff),1),1)); % zero disturbance
                end
                W_for_gray = W_used;            % keep Gray and DDRA aligned at inference
                %% Temporary fix: Overwrites the process noise on gray-box RCSI
                W_for_gray = zonotope(zeros(sys_cora.nrOfDisturbances,1));

                % ================= GRAY =================
                optTS = ts_options_from_pe(C, pe, sys_cora);
                t3 = tic;
                configs = gray_identify(sys_cora, R0_cora, U, C, pe, ...
                    'overrideW', W_for_gray, ...
                    'options_testS', optTS, ...
                    'externalTS_train', TS_train, ...
                    'externalTS_val',   TS_val);  
                Tlearn_g = toc(t3);
                
                want = "graySeq";
                idxGray = find(arrayfun(@(c) isfield(c{1},'method') && want==string(c{1}.method), configs), 1, 'first');
                if isempty(idxGray), idxGray = 1; end

                % Build VAL record (x0,u) from TS_val and pass to both sides
                VAL = VAL_from_TS(TS_val, DATASET_val);
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
                
                if getfielddef(cfg.io,'save_artifacts',true)
                    save(fullfile(artdir, sprintf('row_%04d.mat', row_index)), ...
                         '-struct', 'artifact', '-v7.3');
                end

                want_online = any(pmode == ["online","both"]) && getfielddef(cfg.io,'make_reach_plot',true);
                if want_online && ismember(row_index, getfielddef(cfg.io,'plot_rows',[1]))
                    dims   = getfielddef(cfg.io,'plot_dims',[1 2]);
                    everyK = getfielddef(cfg.io,'plot_every_k', 1);
                
                    sys_true_dt = normalize_to_linearSysDT(sys_ddra, sys_cora.dt);
                
                    if use_noise, W_plot = W_used; else, W_plot = zonotope(zeros(size(center(W_used),1),1)); end
                
                    reach_dir = fullfile(plots_dir,'reach'); if ~exist(reach_dir,'dir'), mkdir(reach_dir); end
                    savename = fullfile(reach_dir, sprintf('row_%04d_y%d%d', row_index, dims(1), dims(2)));
                    optReach  = getfielddef(C.shared,'options_reach', struct());
                    orderCap  = getfielddef(optReach, 'zonotopeOrder', 100);
                     
                    plot_reach_all_onepanel( ...
                        sys_true_dt, ...
                        configs{idxGray}.sys, ...
                        VAL, ...
                        'MAB', M_AB, ...
                        'W', W, ...
                        'Dims', [1 2], ...
                        'ShowSamples', false, ...
                        'Save', cfg.io.save_tag);

                end

    
                % ---- Gray validation on EXACT same points (contains_interval metric)
                % CORA's validateReach (internally splits multi-sample testCases)
                do_plot = any(lower(string(getfielddef(cfg.io,'plot_mode','offline'))) == ["online","both"]);
                ps_args = {};
                if do_plot
                    PS = struct();
                    % Optional: thin the workload
                    % PS.k_plot = 1:cfg.io.plot_every_k:C.shared.n_k_val; 
                    % PS.dims   = getfielddef(cfg.io,'plot_dims',[1 2]);
                    ps_args = {'plot_settings', PS};
                else
                    % IMPORTANT: pass [] to suppress plotting
                    ps_args = {'plot_settings', []};
                end
                
                [ctrain_gray, cval_gray, Tvalidate_g] = gray_containment( ...
                    configs, sys_cora, R0_cora, U, C, pe_eff, ...
                    'externalTS_train', TS_train, ...
                    'externalTS_val',   TS_val,   ...
                    'check_contain',    true,     ...
                    ps_args{:});

              
    
                % ---- DDRA inference (stored or streaming) on SAME validation points
                if LM.store_ddra_sets
                    t2 = tic;
                    % add VAL as an optional arg; legacy code can ignore it
                    [Xsets_ddra, ~] = ddra_infer(sys_ddra, R0, U, W_used, M_AB, C, VAL);
                    Tinfer = toc(t2);

                    Ysets_ddra =  Xsets_ddra;

                    % If you model output noise V as a zonotope, inflate here:
                    % Ysets_ddra = cellfun(@(Y) Y + V, Ysets_ddra, 'uni', 0);
                    
                    wid_ddra = cellfun(@(Z) sum(abs(generators(Z)), 'all'), Ysets_ddra);
                    sizeI_ddra = mean(wid_ddra(:));
                    % Containment should also be done in output-space:
                    %cval_ddra = contains_on_VAL_linear(sys_true_dt, W_used, Ysets_ddra, VAL);  
                    cval_ddra = containsY_on_VAL(Ysets_ddra, VAL, 1e-6);

                    clear Xsets_ddra
                else
                    t2 = tic;
                    [~, cval_ddra, wid_ddra] = ddra_infer_size_streaming(sys_ddra, R0, U, W_used, M_AB, C, VAL);
                    Tinfer = toc(t2);
                    sizeI_ddra = mean(wid_ddra(:)); 
                end

    
                % ---- Gray size metric (interval proxy), consistent with DDRA
                %% Start patch
                % ---- Gray size metric (interval proxy), consistent with DDRA
                t4 = tic;
                sizeI_gray = gray_infer_size_on_VAL(configs{idxGray}.sys, TS_val, C, ...
                    configs{idxGray}.params, 'overrideW', W_for_gray);

                Tinfer_g   = toc(t4);
                
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
                    
                    save(fullfile(artdir, sprintf('row_%04d.mat', row_index)), '-struct','artifact','-v7.3');

                end
                
                clear configs  % free Gray objects early
                
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


                %% end patch
                row.use_noise = use_noise;
                row.ddra_ridge       = ridgeInfo.used;
                row.ddra_lambda      = ridgeInfo.lambda;
                row.ddra_kappa       = ridgeInfo.kappa;
                row.ddra_ridge_policy= char(ridgeInfo.policy);

    
                % --- Initialize schema & write/accumulate ----
                if rowi == 1
                    hdr = fieldnames(orderfields(row))';
                    if LM.append_csv
                        fid = fopen(csv_path, 'w'); fprintf(fid, '%s\n', strjoin(hdr, ',')); fclose(fid);
                    else
                        cells = cell(Ntot, numel(hdr));
                    end
                end
    
                if LM.append_csv
                    append_row_csv(csv_path, hdr, row);
                else
                    for j = 1:numel(hdr)
                        v = row.(hdr{j}); if isstring(v), v = char(v); end
                        cells{rowi, j} = v;
                    end
                end
    
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

end

% ------------------- Helpers-------------------

function use_noise = resolve_use_noise(S)
    % Single switch: default = true unless explicitly disabled
    if isfield(S,'noise_for_ddra') || isfield(S,'noise_for_gray')
        % require both ON to be true; any false -> off
        g = ~isfield(S,'noise_for_gray') || logical(S.noise_for_gray);
        d = ~isfield(S,'noise_for_ddra') || logical(S.noise_for_ddra);
        use_noise = g && d;
    else
        use_noise = true;
    end
end


