function SUMMARY = run_sweeps_lip_rcsi(cfg, grid)
% RUN_SWEEPS_LIP_RCSI
%   Nonlinear/Lipschitz sweep over n_m for kLipMSD using CORA's black-box RCSI (blackCGP).
%   Writes summary.csv and returns a table with (D, n_m, n_s, n_k, cval_black, sizeI_black, times).

    [~, results_dir] = init_io(cfg);
    csv_path = fullfile(results_dir, 'summary.csv');

    LM = getfielddef(cfg,'lowmem',struct());
    append_csv = getfielddef(LM,'append_csv',true);


    % Axes (no PE axis here)
    axes = struct();
    axes.D   = as_list(getfielddef(grid,'D_list',   4));
    axes.n_m = as_list(getfielddef(grid,'n_m_list',[2 4 8 16]));
    axes.n_s = as_list(getfielddef(grid,'n_s_list',getfielddef(cfg,'shared',struct()).n_s));
    axes.n_k = as_list(getfielddef(grid,'n_k_list',getfielddef(cfg,'shared',struct()).n_k));

    % Header
    hdr = {'D','n_m','n_s','n_k', ...
           'cval_black','sizeI_black', ...
           't_black_learn','t_black_val','t_black_infer'};

    if append_csv
        fid = fopen(csv_path,'w'); fprintf(fid,'%s\n',strjoin(hdr,',')); fclose(fid);
    else
        Ntot = numel(axes.D)*numel(axes.n_m)*numel(axes.n_s)*numel(axes.n_k);
        cells = cell(Ntot, numel(hdr));
    end

    rowi = 0;

    for D = axes.D
      for n_m = axes.n_m
        for n_s = axes.n_s
          for n_k = axes.n_k

            % ---------- Build system ----------
            shared = getfielddef(cfg,'shared',struct());
            type   = getfielddef(shared,'type',"standard");
            p_extr = getfielddef(shared,'p_extr',0.3);

            [sys, R0, U, ~] = custom_loadDynamics('kLipMSD', type, struct('k', D));
            options_reach = getfielddef(shared,'options_reach',struct());

            % Horizon for suites
            params_true = struct();
            params_true.R0     = R0;
            params_true.U      = U;
            params_true.tFinal = sys.dt * (n_k - 1);

            ts_opts = struct('p_extr', p_extr);

            % (a) Identification suite (sanity)
            params_true.testSuite = createTestSuite(sys, params_true, n_k, n_m, n_s, ts_opts);

            % (b) TRAIN suite (REQUIRED by black-box conform)
            n_m_train = getfielddef(cfg, {'black','train','n_m'}, 50);
            n_s_train = getfielddef(cfg, {'black','train','n_s'}, 5);
            n_k_train = getfielddef(cfg, {'black','train','n_k'}, n_k);
            params_true.testSuite_train = createTestSuite(sys, params_true, n_k_train, n_m_train, n_s_train, ts_opts);

            % (c) VALIDATION suite (recommended)
            n_m_val = getfielddef(shared,'n_m_val',2);
            n_s_val = getfielddef(shared,'n_s_val',n_s);
            n_k_val = getfielddef(shared,'n_k_val',n_k);
            params_true_val = params_true; % copy tFinal later
            params_true.testSuite_val = createTestSuite(sys, params_true, n_k_val, n_m_val, n_s_val, ts_opts);

            % ---------- Conformance options ----------
            options = options_reach;
            options.cs = getfielddef(shared,'cs_base',struct());
            
            % single ARX order p (use train horizon - 1)
            p = getfielddef(cfg, {'black','approx','p'}, 1);
            
            options.approx = struct( ...
                'cgp_num_gen',        getfielddef(cfg, {'black','approx','cgp_num_gen'}, 5), ...
                'cgp_pop_size_base',  getfielddef(cfg, {'black','approx','cgp_pop_size_base'}, 10), ...
                'gp_parallel',        false, ...
                'save_res',           false, ...
                'p',                  p ...
            );
            
            % ---------- Initial guesses ----------
            cR0 = center(R0);  cU = center(U); cU = cU(:);
            params_id_init        = params_true;   % KEEP testSuite_* fields
            params_id_init.R0     = zonotope([cR0]);                  % 0-radius
            params_id_init.U      = zonotope([cU, eye(numel(cU))]);
            params_id_init.tFinal = sys.dt * (n_k_train - 1);
            params_id_init.n_p    = p;                                     

            % ---------- Run black-box identification (only blackCGP) ----------
            method  = "blackCGP";
            configs = cell(2,1);
            configs{1}.sys     = sys;
            configs{1}.params  = rmfield(params_true,'testSuite'); % keep *_train/val
            configs{1}.options = options_reach;
            configs{1}.name    = "true";

            fprintf("Identification with method %s  (D=%d, n_m=%d, n_s=%d, n_k=%d)\n", method, D, n_m, n_s, n_k);

            t_id = tic;
            [configs{2}.params, results] = conform(sys, params_id_init, options, method);
            t_black_learn = toc(t_id);

            configs{2}.sys     = results.sys;
            configs{2}.options = options_reach;
            configs{2}.name    = method;

            % ---------- Validation containment on testSuite_val ----------
            t_val = tic;
            num_in = 0; num_out = 0;
            for mval = 1:length(params_true.testSuite_val)
                [~, eval_val] = validateReach(params_true.testSuite_val{mval}, configs, 1);
                % index 2 corresponds to learned model
                num_in  = num_in  + eval_val.num_in(2);
                num_out = num_out + eval_val.num_out(2);
            end
            denom = max(1, num_in + num_out);
            cval_black   = 100 * (num_in / denom);
            t_black_val  = toc(t_val);

            % ---------- Conservatism proxy on validation horizon ----------
            t_inf = tic;
            params_true_val.R0     = R0;
            params_true_val.U      = U;
            params_true_val.tFinal = sys.dt * (n_k_val - 1);
            Rlearn = reach(configs{2}.sys, params_true_val, options_reach);
            sizeI_black  = agg_interval_size(Rlearn);
            t_black_infer = toc(t_inf);

            % ---------- Row & write ----------
            rowi = rowi + 1;
            row = {D, n_m, n_s, n_k, ...
                   cval_black, sizeI_black, ...
                   t_black_learn, t_black_val, t_black_infer};

            if append_csv
                fid = fopen(csv_path,'a');
                fprintf(fid, '%s\n', strjoin(cellfun(@num2str_cell, row,'uni',0), ','));
                fclose(fid);
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

% ----------------- helpers -----------------
function L = as_list(v)
    if isscalar(v), L = v; else, L = v(:)'; end
end
function v = getfielddef(S, path, def)
    % supports nested path: {'a','b','c'}
    if ischar(path) || isstring(path), path = {char(path)}; end
    v = S;
    for i=1:numel(path)
        if isstruct(v) && isfield(v, path{i})
            v = v.(path{i});
        else
            v = def; return;
        end
    end
end
function s = num2str_cell(x)
    if ischar(x) || isstring(x), s = char(x);
    elseif islogical(x), s = char(string(x));
    else, s = num2str(x, '%.10g'); end
end
function S = agg_interval_size(R)
    % Sum interval widths over all time points and all state dims
    S = 0;
    if ~isfield(R,'timePoint') || ~isfield(R.timePoint,'set'), return; end
    for k = 1:numel(R.timePoint.set)
        Ik = interval(R.timePoint.set{k});
        S = S + sum(abs(Ik.sup - Ik.inf), 'all');
    end
end
