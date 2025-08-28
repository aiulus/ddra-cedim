function SUMMARY = run_sweeps_square_black_vs_ddraLip(cfg, grid)
% Sweep over n_m for CORA's "Square" using:
%   - Black-box RCSI (conform: GP/CGP)  and
%   - DDRA-Lipschitz (Alg. 6)
%
% Writes rows to CSV and returns SUMMARY table.

    [~, results_dir] = init_io(cfg);
    csv_path = fullfile(results_dir, 'summary.csv');

    LM = getfielddef(cfg,'lowmem',struct());
    append_csv = getfielddef(LM,'append_csv',true);

    % Sweep axes
    axes = struct();
    axes.n_m = as_list(getfielddef(grid,'n_m_list', getfielddef(cfg,'shared',struct()).n_m));
    axes.n_s = as_list(getfielddef(grid,'n_s_list', getfielddef(cfg,'shared',struct()).n_s));
    axes.n_k = as_list(getfielddef(grid,'n_k_list', getfielddef(cfg,'shared',struct()).n_k));
    axes.pe  = getfielddef(grid,'pe_list', {struct('mode','randn')});

    % Header (includes black + ddra columns)
    hdr = {'n_m','n_s','n_k', ...
           'cval_black','sizeI_black', ...
           't_black_learn','t_black_val','t_black_infer', ...
           'cval_ddra','sizeI_ddra', ...
           't_ddra_build','t_ddra_val','t_ddra_reach_avg'};

    if append_csv
        fid = fopen(csv_path,'w'); fprintf(fid,'%s\n',strjoin(hdr,',')); fclose(fid);
    else
        Ntot = numel(axes.n_m)*numel(axes.n_s)*numel(axes.n_k)*numel(axes.pe);
        cells = cell(Ntot, numel(hdr));
    end

    rowi = 0;

    % ------------------------ Main sweeps ------------------------
    for n_m = axes.n_m
      for n_s = axes.n_s
        for n_k = axes.n_k
          for ip = 1:numel(axes.pe)
            pe = axes.pe{ip};

            % -------- Build system & suites --------
            shared = getfielddef(cfg,'shared',struct());
            dyn    = getfielddef(shared,'dyn',"Square");
            type   = getfielddef(shared,'type',"standard");
            p_extr = getfielddef(shared,'p_extr',0.3);

            [sys, R0, U, ~] = custom_loadDynamics(dyn, type);  %# <- required
            options_reach = getfielddef(shared,'options_reach',struct());

            % True params and suites
            params_true = struct();
            params_true.R0     = R0;
            params_true.U      = U;
            params_true.tFinal = sys.dt * (n_k - 1);

            ts_opts = struct('p_extr', p_extr);
            params_true.testSuite       = createTestSuite(sys, params_true, n_k, n_m, n_s, ts_opts);

            n_m_train = getfielddef(cfg, {'black','train','n_m'}, 100);
            n_s_train = getfielddef(cfg, {'black','train','n_s'}, 10);
            n_k_train = getfielddef(cfg, {'black','train','n_k'}, n_k);
            params_true.testSuite_train = createTestSuite(sys, params_true, n_k_train, n_m_train, n_s_train, ts_opts);

            n_m_val = getfielddef(shared,'n_m_val', max(2,min(n_m,5)));
            n_s_val = getfielddef(shared,'n_s_val', n_s);
            n_k_val = getfielddef(shared,'n_k_val', n_k);
            params_true.testSuite_val = createTestSuite(sys, params_true, n_k_val, n_m_val, n_s_val, ts_opts);

            % ======== BLACK-BOX RCSI (conform) ========
            options = options_reach;
            options.cs = getfielddef(shared,'cs_base',struct());
            p = getfielddef(cfg, {'black','approx','p'}, getfielddef(sys,'n_p',1));
            options.approx = struct('p', p);

            apx = getfielddef(cfg,'black',struct()); apx = getfielddef(apx,'approx',struct());
            map_if(options, apx, 'gp_parallel');       map_if(options, apx, 'gp_pop_size');
            map_if(options, apx, 'gp_num_gen');        map_if(options, apx, 'gp_func_names');
            map_if(options, apx, 'gp_max_genes');      map_if(options, apx, 'gp_max_depth');
            map_if(options, apx, 'cgp_num_gen');       map_if(options, apx, 'cgp_pop_size_base');
            map_if(options, apx, 'save_res');

            % Initial guesses (centered sets)
            cR0 = center(R0);
            cU  = center(U);  cU = cU(:);
            params_id_init        = params_true;          % keep testSuite_* fields
            params_id_init.R0     = zonotope([cR0]);      % 0-radius
            params_id_init.U      = zonotope([cU, eye(numel(cU)) ones(numel(cU),1)]);
            params_id_init.n_p    = p;

            methods = getfielddef(cfg, {'black','methodsBlack'}, ["blackCGP"]);
            t_black_learn = NaN; learned_idx = 2;
            configs = cell(numel(methods)+1,1);

            % configs{1}: true
            configs{1}.sys     = sys;
            configs{1}.params  = rmfield(params_true,'testSuite');
            configs{1}.options = options_reach;
            configs{1}.name    = "true";

            for iM = 1:numel(methods)
                method = string(methods(iM));
                fprintf("Identification (BLACK %s): n_m=%d, n_s=%d, n_k=%d\n", method, n_m, n_s, n_k);
                t_id = tic;
                [configs{iM+1}.params, results] = conform(sys, params_id_init, options, method);
                if iM == 1, t_black_learn = toc(t_id); else, ~toc(t_id); end
                configs{iM+1}.sys     = results.sys;
                configs{iM+1}.options = options_reach;
                configs{iM+1}.name    = method;
            end

            % Validation (BLACK) — containment on validation suite
            t_val = tic;
            num_in = 0; num_out = 0;
            for mval = 1:length(params_true.testSuite_val)
                [~, eval_val] = validateReach(params_true.testSuite_val{mval}, configs, 1);
                num_in  = num_in  + eval_val.num_in(learned_idx);
                num_out = num_out + eval_val.num_out(learned_idx);
            end
            denom = max(1, num_in + num_out);
            cval_black   = 100 * (num_in / denom);
            t_black_val  = toc(t_val);

            % Conservatism proxy (BLACK) — size of intervals for learned model
            t_inf = tic;
            params_val = params_true;
            params_val.tFinal = sys.dt * (n_k_val - 1);
            Rlearn = reach(configs{learned_idx}.sys, params_val, options_reach);
            sizeI_black  = agg_interval_size(Rlearn);
            t_black_infer = toc(t_inf);

            % ======== DDRA-Lipschitz (Alg. 6) ========
            % Build data D from training suite
            p = getfielddef(cfg, {'black','approx','p'}, 1);   % Square: p=1
            D = build_D_from_suite(params_true.testSuite_train, p);
            nx = size(D.Xminus,1);

            % Disturbance W (possibly off for Black)
            if getfielddef(cfg.shared,'noise_for_ddra', true)
                eta_w  = getfielddef(cfg.ddra,'eta_w',1);
                alphaW = getfielddef(cfg.ddra,'alpha_w',0);
                W = zonotope(zeros(nx,1), alphaW*ones(nx, max(1,eta_w)));
            else
                W = zonotope(zeros(nx,1));
            end

            % Estimate Lipschitz & cover radius
            t_ddra_build_tic = tic;
            [Lvec, delta] = estimate_L_and_delta_from_data(D);
            Cnl = struct();
            Cnl.shared = struct('n_k', n_k, 'n_k_val', n_k_val);
            Cnl.lowmem = struct('zonotopeOrder_cap', getfielddef(cfg.lowmem,'zonotopeOrder_cap',100));
            Cnl.nlip   = struct('ZepsFlag', true, 'Lvec', Lvec(:), 'gamma', delta*ones(nx,1));
            t_ddra_build = toc(t_ddra_build_tic);

            % Validation with DDRA-Lip
            t_ddra_val_tic = tic;
            num_in_d = 0; num_out_d = 0; sum_sizeI_ddra = 0; t_reach_sum = 0;
            for mval = 1:length(params_true.testSuite_val)
                TS = params_true.testSuite_val{mval};
                U_nom = params_true.testSuite_val{mval}.u;   % [nu x n_k_val]
                Uk = cell(1, n_k_val);
                for k = 1:n_k_val
                    u_k   = U_nom(:,k);
                    u_km1 = U_nom(:, max(1, k-1));  % pad with u_1 at k=1
                    u_stack = [u_k; u_km1];
                    Uk{k} = zonotope(u_stack, zeros(numel(u_stack),0)); % singleton set
                end

                t_case = tic;
                [Xsets, sizeI_case] = ddra_reach_lipschitz(params_true.R0, Uk, W, D, Cnl);
                t_reach_sum = t_reach_sum + toc(t_case);

                [nin, nout] = contain_points_in_sets(TS.y, Xsets);
                num_in_d  = num_in_d  + nin;
                num_out_d = num_out_d + nout;
                sum_sizeI_ddra = sum_sizeI_ddra + sizeI_case;
            end
            t_ddra_val       = toc(t_ddra_val_tic);
            t_ddra_reach_avg = t_reach_sum / max(1,length(params_true.testSuite_val));
            cval_ddra        = 100 * (num_in_d / max(1, num_in_d + num_out_d));
            sizeI_ddra       = sum_sizeI_ddra / max(1,length(params_true.testSuite_val));

            % -------- Row & write --------
            rowi = rowi + 1;
            row = {n_m, n_s, n_k, ...
                   cval_black, sizeI_black, ...
                   t_black_learn, t_black_val, t_black_infer, ...
                   cval_ddra, sizeI_ddra, ...
                   t_ddra_build, t_ddra_val, t_ddra_reach_avg};

            if append_csv
                fid = fopen(csv_path,'a');
                fprintf(fid, '%s\n', strjoin(cellfun(@num2str_cell, row,'uni',0), ',')); fclose(fid);
            else
                cells(rowi,:) = row;
            end

          end
        end
      end
    end

    if ~append_csv
        writecell([hdr; cells(1:rowi,:)], csv_path);
    end
    SUMMARY = readtable(csv_path);
    fprintf('Sweeps done. Rows: %d\nCSV -> %s\n', height(SUMMARY), csv_path);
end