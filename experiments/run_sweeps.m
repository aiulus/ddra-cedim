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
                use_noise   = resolve_use_noise(C.shared);
    
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
                    'externalTS_train', TS_train, ...
                    'externalTS_val',   TS_val);  
                Tlearn_g = toc(t3);
    
                % Build VAL record (x0,u) from TS_val and pass to both sides
                VAL = local_VAL_from_TS(TS_val, DATASET_val);
    
                % ---- Gray validation on EXACT same points (contains_interval metric)
                [ctrain_gray, cval_gray, Tvalidate_g] = gray_containment_on_VAL(configs, TS_train, TS_val, 1e-6);
              
    
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
                %% Start patch
                % ---- Gray size metric (interval proxy), consistent with DDRA
                t4 = tic;
                sizeI_gray = gray_infer_size_on_VAL(configs{2}.sys, TS_val, C, configs{2}.params);
                Tinfer_g   = toc(t4);
                
                % --- Optional: save plotting artifact for this row ---
                if isfield(cfg,'io') && isfield(cfg.io,'save_artifacts') && cfg.io.save_artifacts
                    % ensure folder
                    artdir = fullfile(results_dir, 'artifacts');
                    if ~exist(artdir,'dir'), mkdir(artdir); end
                
                    % pack only what plotting needs
                    artifact = struct();
                    if numel(configs) >= 2
                        artifact.sys_gray = configs{2}.sys;    % identified RCSI/Gray model (CORA linearSysDT)
                    else
                        artifact.sys_gray = configs{1}.sys;    % fallback
                    end
                    artifact.sys_ddra = sys_ddra;              % true/nominal model (A,B,C,D,dt struct)
                    artifact.VAL      = VAL;                   % exact (x0,u,y[,w]) used in validation
                    artifact.W_eff    = W_eff;                 % disturbance used by DDRA
                    artifact.meta = struct( ...
                        'rowi', rowi+1, 'D', D, 'alpha_w', alpha_w, ...
                        'n_m', C.shared.n_m, 'n_s', C.shared.n_s, 'n_k', C.shared.n_k, ...
                        'pe',  pe, 'save_tag', cfg.io.save_tag);
                
                    % include R0 for nicer initial-set plotting (if not already present)
                    if ~isfield(artifact.VAL,'R0') || isempty(artifact.VAL.R0)
                        artifact.VAL.R0 = R0;
                    end
                
                    artfile = fullfile(artdir, sprintf('row_%04d.mat', rowi+1));
                    save(artfile, '-struct', 'artifact');
                end
                
                clear configs  % free Gray objects early
                
                % --- Pack row ---
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

function TS = local_TS_from_dataset(sys_cora, ~, DATASET, n_k_override)
% Build CORA testCase objects from DATASET.
% Accepts BOTH new (x0_list/U_blocks/Y_blocks) and old (x0_blocks/u_blocks) layouts.

    % --- detect layout & define accessors ---
    if isfield(DATASET,'x0_list') && isfield(DATASET,'U_blocks')
        M       = size(DATASET.x0_list, 2);
        get_x0  = @(b) DATASET.x0_list(:, b);              % (n_x×1)
        get_U   = @(b) DATASET.U_blocks(:, :, b);          % (n_u×n_k)
        have_Y  = isfield(DATASET,'Y_blocks');
        if have_Y, get_Y = @(b) DATASET.Y_blocks(:, :, b); % (n_y×n_k)
        end
        n_k_ds  = size(DATASET.U_blocks, 2);
    elseif isfield(DATASET,'x0_blocks') && isfield(DATASET,'u_blocks')
        M       = numel(DATASET.x0_blocks);
        get_x0  = @(b) DATASET.x0_blocks{b};               % (n_x×1)
        get_U   = @(b) DATASET.u_blocks{b};                % (n_u×n_k)
        have_Y  = false;
        n_k_ds  = size(get_U(1), 2);
    else
        error('local_TS_from_dataset: DATASET missing x0/U fields.');
    end

    % --- pick horizon ---
    if nargin >= 4 && ~isempty(n_k_override)
        n_k = n_k_override;
    else
        n_k = n_k_ds;
    end

    % --- construct testCase objects ---
    ny  = sys_cora.nrOfOutputs;
    TS  = cell(1, M);
    for b = 1:M
        x0   = get_x0(b);         % (n_x×1)
        Ublk = get_U(b);          % (n_u×n_k)

        % match requested horizon by crop/pad (rare)
        if size(Ublk,2) > n_k, Ublk = Ublk(:,1:n_k); end
        if size(Ublk,2) < n_k, Ublk = [Ublk, zeros(size(Ublk,1), n_k-size(Ublk,2))]; end

        % CORA testCase expects time-major u: (a×p×s) = (n_k × n_u × 1)
        u_tc = reshape(Ublk.', n_k, size(Ublk,1), 1);

        % Measured outputs y: (a×q×s) = (n_k × n_y × 1)
        if have_Y
            %Yblk = get_Y(b);            % (n_y×n_k)
            %y_tc = reshape(Yblk.', n_k, ny, 1);
            y_tc = permute(DATASET.Y_blocks(:,:,b), [2 1 3]);
        else
            % Fallback: if DATASET has no Y, hand in NaNs (RCSI can still reach/compare)
            %y_tc = nan(n_k, ny, 1);
            y_tc = nan(n_k, ny, 1);
        end

        % Create the actual CORA testCase object
        TS{b} = testCase(y_tc, u_tc, x0, sys_cora.dt, class(sys_cora));
    end
end

function VAL = local_VAL_from_TS(TS_val, DATASET_val)
    % Build synchronized validation record for EXACT (x0,u,y[,w]) points.

    B = numel(TS_val);
    VAL = struct('x0',{cell(1,B)}, 'u',{cell(1,B)}, 'y',{cell(1,B)}, 'w',{cell(1,B)});

    have_DS = (nargin >= 2) && ~isempty(DATASET_val);

    for b = 1:B
        % x0
        VAL.x0{b} = TS_val{b}.initialState;

        % u : (m × n_k)
        Ui = squeeze(TS_val{b}.u);         % (n_k × m) or (n_k × m × 1)
        if ndims(TS_val{b}.u) == 3, Ui = squeeze(TS_val{b}.u); end
        if size(Ui,1) < size(Ui,2), Ui = Ui.'; end    % ensure (m × n_k)
        VAL.u{b} = Ui;

        % y : prefer DATASET_val if present; otherwise use the TS field
        yi = [];
        if have_DS
            if isfield(DATASET_val,'Y_blocks')           % new layout: (ny × n_k × M)
                yi = DATASET_val.Y_blocks(:,:,b).';      % (n_k × ny)
            elseif isfield(DATASET_val,'y_blocks')       % rare legacy naming
                yi = DATASET_val.y_blocks(:,:,b).';
            end
        end
        if isempty(yi) && isfield(TS_val{b},'y') && ~isempty(TS_val{b}.y)
            ti = TS_val{b}.y;                            % (n_k × ny × 1 or s)
            yi = squeeze(ti);                            % (n_k × ny) if s==1
            if isvector(yi), yi = yi(:); end
        end
        VAL.y{b} = yi;  

        if have_DS && isfield(DATASET_val,'W_blocks')
            VAL.w{b} = DATASET_val.W_blocks(:,:,b);      % (nx × n_k)
        end
    end
end
