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
                C.ddra.alpha_w   = alpha_w;               % match noise scale
                C.shared.seed    = 1; rng(C.shared.seed, 'twister');
                if C.shared.dyn == "k-Mass-SD"; C.shared.dyn_p = D; end

                % --- Build true systems ---
                [sys_cora, sys_ddra, R0, U] = build_true_system(C);

                % ================= DDRA =================
                t0 = tic;
                [Xminus, Uminus, Xplus, W, Zinfo] = ddra_generate_data(sys_ddra, R0, U, C, pe);
                Tlearn = toc(t0);

                t1 = tic;
                M_AB = ddra_learn_Mab(Xminus, Uminus, Xplus, W, Zinfo, sys_ddra);
                Tcheck = toc(t1);

                if LM.store_ddra_sets
                    t2 = tic;
                    [Xsets_ddra, sizeI_ddra] = ddra_infer(sys_ddra, R0, U, W, M_AB, C);
                    Tinfer = toc(t2);
                    cval_ddra = mc_containment_X(sys_ddra, R0, U, W, Xsets_ddra, C);
                    clear Xsets_ddra
                else
                    t2 = tic;
                    [sizeI_ddra, cval_ddra] = ddra_infer_size_streaming(sys_ddra, R0, U, W, M_AB, C);
                    Tinfer = toc(t2);
                end
                clear Xminus Uminus Xplus W  % free data blocks early

                % ================= GRAY =================

                % ================= Start: New Patch =================
                % Force W=0 for Gray/RCSI if requested (keeps fair comparison when not studying noise)
                W_for_gray = W;
                if isfield(C.shared,'noise_for_gray') && ~C.shared.noise_for_gray
                    W_for_gray = zonotope(zeros(size(center(W),1),1)); % zero disturbance
                end

                t3 = tic;
                configs = gray_identify(sys_cora, R0, U, C, pe, 'overrideW', W_for_gray);
                Tlearn_g = toc(t3);
                % ================= End: New Patch =================

                %t3 = tic;
                %configs  = gray_identify(sys_cora, R0, U, C, pe);
                %Tlearn_g = toc(t3);

                [ctrain_gray, cval_gray, Tvalidate_g] = gray_containment( ...
                    configs, sys_cora, R0, U, C, pe, 'check_contain', LM.gray_check_contain);

                t4 = tic;
                sizeI_gray = gray_infer_size(configs{2}.sys, R0, U, C);
                Tinfer_g   = toc(t4);

                clear configs  % free Gray objects early

                % --- Pack row ---
                rowi = rowi + 1;
                row = pack_row(C, D, alpha_w, pe, ...
                    ctrain_gray, cval_gray, cval_ddra, sizeI_ddra, sizeI_gray, ...
                    Zinfo.rankZ, Zinfo.condZ, ...
                    Tlearn, Tcheck, Tinfer, Tlearn_g, Tvalidate_g, Tinfer_g);

                % add PE diagnostics if present
                row.pe_hankel_rank = getfielddef(Zinfo,'hankel_rank', NaN);
                row.pe_hankel_full = getfielddef(Zinfo,'hankel_full', false);
                row.pe_rank_frac   = getfielddef(Zinfo,'rank_frac',   NaN);

                % --- Initialize schema on first row & write/accumulate ----
                if rowi == 1
                    hdr = fieldnames(orderfields(row))';
                    if LM.append_csv
                        % header once
                        fid = fopen(csv_path, 'w');
                        fprintf(fid, '%s\n', strjoin(hdr, ','));
                        fclose(fid);
                    else
                        cells = cell(Ntot, numel(hdr));  % fully preallocate
                    end
                end

                if LM.append_csv
                    append_row_csv(csv_path, hdr, row);   % stream to disk
                else
                    for j = 1:numel(hdr)
                        v = row.(hdr{j});
                        if isstring(v), v = char(v); end
                        cells{rowi, j} = v;
                    end
                end

                % free leftover biggies (paranoia)
                clear M_AB Zinfo sizeI_ddra sizeI_gray
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
