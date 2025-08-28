function SUMMARY = run_sweeps(cfg, grid)
% RUN_SWEEPS  Batch experiments for DDRA vs Gray under controlled sweeps.
% Produces CSV summary in results_dir and plots in plots_dir (via callers).

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
    
                % --- instantiate config for this run ---
                C = baseC;
                C.shared.n_m = n_m;  C.shared.n_s = n_s;  C.shared.n_k = n_k;
                C.shared.n_m_val = max(2, min(n_m, getfielddef(baseC.shared, 'n_m_val', 2)));
                C.shared.n_s_val = n_s;  C.shared.n_k_val = n_k;
                C.ddra.alpha_w   = alpha_w;
                C.shared.seed    = 1; rng(C.shared.seed, 'twister');
                if C.shared.dyn == "k-Mass-SD"; C.shared.dyn_p = D; end
    
                % --- Build true systems ---
                [sys_cora, sys_ddra, R0, U] = build_true_system(C);
    
                % ================= DDRA =================
                t0 = tic;
                % NOTE: ddra_generate_data must return DATASET as 6th out:
                % DATASET.x0_blocks{b}, DATASET.u_blocks{b} with block length n_k
                [Xminus, Uminus, Xplus, W, Zinfo, DATASET] = ddra_generate_data(sys_ddra, R0, U, C, pe);
                Tlearn = toc(t0);
    
                % Build TRAIN and VAL suites deterministically from generator
                TS_train = local_TS_from_dataset(sys_cora, R0, DATASET, C.shared.n_k);
    
                C_val = C;
                C_val.shared.n_m = C.shared.n_m_val;
                C_val.shared.n_s = C.shared.n_s_val;
                C_val.shared.n_k = C.shared.n_k_val;
                [~,~,~,~,~, DATASET_val] = ddra_generate_data(sys_ddra, R0, U, C_val, pe);
                TS_val   = local_TS_from_dataset(sys_cora, R0, DATASET_val, C_val.shared.n_k);
    
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
                
                    % --- write row (unchanged from your code) ---
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

    
                clear Xminus Uminus Xplus  % free data blocks early
    
                % ---- unified disturbance policy ----
                use_noise   = resolve_use_noise(C.shared);
                if use_noise
                    W_for_gray = W;
                else
                    W_for_gray = zonotope(zeros(size(center(W),1),1));
                end

                % ================= GRAY =================
                optTS = ts_options_from_pe(C, pe, sys_cora);
                t3 = tic;
                configs = gray_identify(sys_cora, R0, U, C, pe, ...
                    'overrideW', W_for_gray, ...
                    'options_testS', optTS, ...
                    'externalTS_train', TS_train);
                Tlearn_g = toc(t3);
    
                % Build VAL record (x0,u) from TS_val and pass to both sides
                VAL = local_VAL_from_TS(TS_val);
    
                % ---- Gray validation on EXACT same points (contains_interval metric)
                [ctrain_gray, cval_gray, Tvalidate_g] = gray_containment_on_VAL( ...
                    configs, TS_train, TS_val);  % thin wrapper calling reach() with cs stripped
    
                % ---- DDRA inference (stored or streaming) on SAME validation points
                if LM.store_ddra_sets
                    t2 = tic;
                    [Xsets_ddra, sizeI_ddra] = ddra_infer(sys_ddra, R0, U, W_eff, M_AB, C);
                    Tinfer = toc(t2);
                    % Recompute containment with EXACT VAL points instead of MC
                    cval_ddra = contains_on_VAL_linear(sys_ddra, W_eff, Xsets_ddra, VAL);
                    clear Xsets_ddra
                else
                    t2 = tic;
                    [sizeI_ddra, cval_ddra] = ddra_infer_size_streaming(sys_ddra, R0, U, W_eff, M_AB, C, VAL);
                    Tinfer = toc(t2);
                end
    
                % ---- Gray size metric (interval proxy), consistent with DDRA
                t4 = tic;
                sizeI_gray = gray_infer_size_on_VAL(configs{2}.sys, TS_val, C);
                Tinfer_g   = toc(t4);
    
                clear configs  % free Gray objects early
    
                % --- Pack row ---
                rowi = rowi + 1;
                row = pack_row(C, D, alpha_w, pe, ...
                    ctrain_gray, cval_gray, cval_ddra, sizeI_ddra, sizeI_gray, ...
                    Zinfo.rankZ, Zinfo.condZ, ...
                    Tlearn, Tcheck, Tinfer, Tlearn_g, Tvalidate_g, Tinfer_g);
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

function TS = local_TS_from_dataset(sys_cora, R0, DATASET, n_k)
    % Build a cell array of CORA testCase-like structs with fields: y,u,initialState
    % using the exact (x0,u) blocks recorded in DATASET.
    % Assumes you already have a simulator that produced y to match your identification code.
    B = numel(DATASET.x0_blocks);
    TS = cell(1,B);
    for b = 1:B
        obj = struct();
        obj.initialState = DATASET.x0_blocks{b};              % (dim_x × 1)
        obj.u            = reshape(DATASET.u_blocks{b}, [], n_k, 1); % (m × n_k × 1)
        % y is only needed by validate-like plotters; Gray/size code can simulate it
        obj.y            = nan(n_k, sys_cora.nrOfOutputs, 1);
        TS{b} = obj;
    end
end

function VAL = local_VAL_from_TS(TS_val)
    % Normalize validation points for “exact same points” containment checks.
    B = numel(TS_val);
    VAL = struct('x0',{cell(1,B)}, 'u',{cell(1,B)});
    for b = 1:B
        VAL.x0{b} = TS_val{b}.initialState;
        VAL.u{b}  = squeeze(TS_val{b}.u);  % (m × n_k)
        if isvector(VAL.u{b}), VAL.u{b} = VAL.u{b}(:)'; end
    end
end
