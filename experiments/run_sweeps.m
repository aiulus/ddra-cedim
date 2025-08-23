function SUMMARY = run_sweeps(cfg, grid)
    % RUN_SWEEPS  Batch experiments for DDRA vs Gray under controlled sweeps.
    % Produces CSV summary + overview plots in results_dir/plots_dir.
    
    [plots_dir, results_dir] = init_io(cfg);
    [axes, baseC] = init_sweep_axes(cfg, grid);
    
    rows = struct([]);   
    rowi = 0;
    
    for D = axes.D
      for alpha_w = axes.alpha_w
        for n_m = axes.n_m
          for n_s = axes.n_s
            for n_k = axes.n_k
              for ip = 1:numel(axes.pe)
                pe = axes.pe{ip};
    
                % --- instantiate config for this run ---
                C = baseC; C.shared.n_m = n_m; C.shared.n_s = n_s; C.shared.n_k = n_k;
                C.shared.n_m_val = max(2, min(n_m, getfielddef(baseC.shared,'n_m_val',2)));
                C.shared.n_s_val = n_s; C.shared.n_k_val = n_k;
                C.ddra.alpha_w = alpha_w;  % match noise scale
                C.shared.seed  = 1; rng(C.shared.seed,'twister');
                if C.shared.dyn == "k-Mass-SD"; C.shared.dyn_p = D; end
    
                % --- Build true systems ---
                [sys_cora, sys_ddra, R0, U] = build_true_system(C);
    
                % ================= DDRA =================
                t0 = tic; [Xminus, Uminus, Xplus, W, Zinfo] = ddra_generate_data(sys_ddra, R0, U, C, pe); Tlearn = toc(t0);
                t1 = tic; M_AB = ddra_learn_Mab(Xminus, Uminus, Xplus, W, Zinfo, sys_ddra); Tcheck = toc(t1);
                t2 = tic; [Xsets_ddra, sizeI_ddra] = ddra_infer(sys_ddra, R0, U, W, M_AB, C); Tinfer = toc(t2);
                cval_ddra = mc_containment_X(sys_ddra, R0, U, W, Xsets_ddra, C); % fidelity
    
                % ================= GRAY =================
                t3 = tic; configs = gray_identify(sys_cora, R0, U, C, pe); Tlearn_g = toc(t3);
                [ctrain_gray, cval_gray, Tvalidate_g] = gray_containment(configs, sys_cora, R0, U, C, pe);
                t4 = tic; sizeI_gray = gray_infer_size(configs{2}.sys, R0, U, C); Tinfer_g = toc(t4);
    
                % --- Log a row ---
                rowi = rowi + 1;
                rows(rowi) = pack_row(C, D, alpha_w, pe, ...
                    ctrain_gray, cval_gray, cval_ddra, sizeI_ddra, sizeI_gray, ...
                    Zinfo.rankZ, Zinfo.condZ, ...
                    Tlearn, Tcheck, Tinfer, Tlearn_g, Tvalidate_g, Tinfer_g); %#ok<AGROW>
    
                % Overview text artifact
                save_overview_text(plots_dir, sprintf('D%d_nm%d_ns%d_nk%d_pe%s',D,n_m,n_s,n_k,pe.mode), ...
                    rows(rowi));
    
              end
            end
          end
        end
      end
    end
    
    % --- Save CSV ---
    SUMMARY = struct2table(rows);
    writetable(SUMMARY, fullfile(results_dir, 'summary.csv'));
    fprintf('Sweeps done. Rows: %dCSV -> %s', height(SUMMARY), fullfile(results_dir,'summary.csv'));
end