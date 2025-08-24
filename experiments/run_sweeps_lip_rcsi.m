function SUMMARY = run_sweeps_lip_rcsi(cfg, grid)
% RUN_SWEEPS_LIP_RCSI
%   Nonlinear/Lipschitz sweep over n_m for kLipMSD using CORA's black-box RCSI.
%   Writes summary.csv and returns a table with (n_m, D, cval_black, sizeI_black, times,...).

    [plots_dir, results_dir] = init_io(cfg);
    csv_path = fullfile(results_dir, 'summary.csv');

    LM = getfielddef(cfg,'lowmem',struct());    
    append_csv = getfielddef(LM,'append_csv',true);

    % Axes
    axes = struct();
    axes.D   = as_list(getfielddef(grid,'D_list',4));
    axes.n_m = as_list(getfielddef(grid,'n_m_list',[2 4 8 16]));
    axes.n_s = as_list(getfielddef(grid,'n_s_list',cfg.shared.n_s));
    axes.n_k = as_list(getfielddef(grid,'n_k_list',cfg.shared.n_k));
    axes.pe  = grid.pe_list;  if ~iscell(axes.pe), axes.pe = {axes.pe}; end

    % header
    hdr = {'D','n_m','n_s','n_k','pe_mode','pe_order', ...
           'cval_black','sizeI_black', ...
           't_black_learn','t_black_val','t_black_infer'};

    if ~append_csv
        cells = cell(numel(axes.D)*numel(axes.n_m)*numel(axes.pe), numel(hdr));
    else
        % fresh file w/ header
        fid = fopen(csv_path,'w'); fprintf(fid,'%s\n',strjoin(hdr,',')); fclose(fid);
    end

    rowi = 0;
    for D = axes.D
      for n_m = axes.n_m
        for n_s = axes.n_s
          for n_k = axes.n_k
            for ip = 1:numel(axes.pe)
                pe = axes.pe{ip};
                % ---- Build system + suites ----
                [sys, R0, U, ~] = custom_loadDynamics('kLipMSD', cfg.shared.type, struct('k',D));
                params_true = struct('R0',R0,'U',U,'tFinal', sys.dt*(n_k-1));
                params_true.testSuite = createTestSuite(sys, params_true, n_k, n_m, n_s, struct('p_extr',cfg.shared.p_extr));
                params_true.tFinal = sys.dt*(cfg.shared.n_k_val-1);
                TS_val = createTestSuite(sys, params_true, cfg.shared.n_k_val, cfg.shared.n_m_val, cfg.shared.n_s_val, struct('p_extr',cfg.shared.p_extr));

                % ---- Identification options (CORA-style) ----
                options_reach = cfg.shared.options_reach;
                options = options_reach;
                options.cs = cfg.shared.cs_base;
                options.cs.constraints = string(getfielddef(options.cs,'constraints',"half"));

                % Black-box options (light defaults; override via cfg.black.approx)
                approx = struct();
                approx.gp_parallel       = false;
                approx.gp_pop_size       = 50;
                approx.gp_num_gen        = 30;
                approx.gp_func_names     = {'times','plus','square'};
                approx.gp_max_genes      = 2;
                approx.gp_max_depth      = 2;
                approx.cgp_num_gen       = 5;
                approx.cgp_pop_size_base = 5;
                approx.save_res          = false;
                approx.p                 = getfielddef(sys,'n_p',1); % ARX/NARX depth if needed
                if isfield(cfg,'black') && isfield(cfg.black,'approx')
                    fn = fieldnames(cfg.black.approx);
                    for ii=1:numel(fn), approx.(fn{ii}) = cfg.black.approx.(fn{ii}); end
                end
                options.approx = approx;

                % ---- Pack configs for conform() ----
                configs = cell(2,1);
                configs{1}.sys = sys;
                configs{1}.params = rmfield(params_true,'testSuite');
                configs{1}.options = options_reach;
                configs{1}.name = "true";

                % Harmonize identification sets expected by conform()
                params_id_init = params_true;
                params_id_init.R0 = R0;                 % state-space, non-ARX
                cU = center(U);                         % keep dimension
                params_id_init.U = zonotope([cU, eye(size(cU,1)), ones(size(cU))]);

                % ---- Learn black-box (first method only for this sweep) ----
                method = string(cfg.black.methodsBlack(1));
                t0 = tic;
                [configs{2}.params, results] = conform(sys, params_id_init, options, method);
                Tlearn = toc(t0);
                configs{2}.sys     = results.sys;
                configs{2}.options = options_reach;
                configs{2}.name    = method;

                % ---- Validation fidelity on TS_val ----
                num_out=zeros(2,1); num_in=zeros(2,1);
                t1 = tic;
                for mvi = 1:length(TS_val)
                    [R, eval_val] = validateReach(TS_val{mvi}, configs, 1); 
                    num_out = num_out + eval_val.num_out;
                    num_in  = num_in  + eval_val.num_in;
                end
                Tval = toc(t1);
                num_all = length(TS_val)*cfg.shared.n_k_val*size(TS_val{1}.y,3);
                cval_black = 100*(1 - num_out(2)/(num_out(2)+num_in(2)));

                % ---- Conservatism proxy (aggregate interval size) ----
                % Use a single model-based reach on the learned sys over the same horizon
                t2 = tic;
                P = configs{2}.params;         % contains R0,U,tFinal
                Rlearn = reach(configs{2}.sys, P, options_reach);
                sizeI_black = agg_interval_size(Rlearn);
                Tinfer = toc(t2);

                % ---- Row & write ----
                rowi = rowi + 1;
                row = {D, n_m, n_s, n_k, string(getfielddef(pe,'mode','randn')), ...
                       getfielddef(pe,'order',1), ...
                       cval_black, sizeI_black, ...
                       Tlearn, Tval, Tinfer};

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
function v = getfielddef(S,f,def)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = def; end
end
function s = num2str_cell(x)
    if ischar(x) || isstring(x), s = char(x);
    elseif islogical(x), s = string(x); s = s{1};
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
