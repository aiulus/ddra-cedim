function [stats, R_last] = ident_and_validate_blackbox(lookup, conf_opts)
% Lightweight version returning validation metrics for the learned model.

    rng(1,'twister');

    [sys, params_true.R0, params_true.U] = custom_loadDynamics(lookup.sys.dyn, lookup.sys.type);
    params_true.tFinal = sys.dt * lookup.n_k - sys.dt;

    params_true.testSuite        = createTestSuite(sys, params_true, lookup.n_k,       lookup.n_m,      lookup.n_s, struct('p_extr',0.3));
    params_true.testSuite_train  = createTestSuite(sys, params_true, lookup.n_k_train, lookup.n_m_train, lookup.n_s_train);
    params_true.testSuite_val    = createTestSuite(sys, params_true, lookup.n_k_val,   lookup.n_m_val,   lookup.n_s_val, struct('p_extr',0.3));

    options = conf_opts.options_reach;
    options.cs = conf_opts.cs; options.cs.constraints = lookup.constraints;
    options.approx = conf_opts.approx;

    configs = cell(2,1);
    configs{1}.sys     = sys;
    configs{1}.params  = rmfield(params_true,'testSuite');
    configs{1}.options = conf_opts.options_reach;
    configs{1}.name    = "true";

    % ARX: 0-dim R0 + interval U around center
    params_id_init = params_true;
    params_id_init.R0 = zonotope(zeros(0,1));
    dim_u = sys.nrOfInputs; cU = center(params_true.U);
    if numel(cU) ~= dim_u, cU = zeros(dim_u,1); end
    params_id_init.U = zonotope([cU(:), eye(dim_u)]);
    if ~isfield(options,'approx') || ~isfield(options.approx,'p') || isempty(options.approx.p) || options.approx.p==0
        options.approx.p = max(1, getfielddef(sys,'n_p',1));
    end

    
    nx_eff = effective_state_dim(sys, options);  % <-- helper below
    nu     = sys.nrOfInputs;
    
    % R0: center at 0 in the correct dimension if missing/mismatched
    if ~isfield(params_id_init,'R0') || isempty(params_id_init.R0) ...
       || size(center(params_id_init.R0),1) ~= nx_eff
        params_id_init.R0 = zonotope(zeros(nx_eff,1));
    end
    
    % U: ensure correct input dimension; keep center if available
    cU = zeros(nu,1);
    if isfield(params_id_init,'U') && ~isempty(params_id_init.U)
        cU0 = center(params_id_init.U);
        if numel(cU0) == nu, cU = cU0(:); end
    end
    params_id_init.U = zonotope([cU, eye(nu), ones(nu,1)]);
    t_learn = tic;
    [configs{2}.params, results] = conform(sys, params_id_init, options, "blackGP");
    t_learn = toc(t_learn);

    configs{2}.sys     = results.sys;
    configs{2}.options = conf_opts.options_reach;
    configs{2}.name    = "blackGP";

    % Validation
    num_out = 0; num_in = 0; size_acc = 0; steps_tot = 0;

    t_check = tic;
    for m=1:length(params_true.testSuite_val)
        [R_last, eval_val] = validateReach(params_true.testSuite_val{m}, configs, 1);
        num_out = num_out + eval_val.num_out; 
        num_in  = num_in  + eval_val.num_in;

        % size proxy for learned model only (index 2)
        for k=1:numel(R_last{2}.timePoint.set)
            S = R_last{2}.timePoint.set{k};
            box = interval(S); width = box.sup - box.inf;
            size_acc = size_acc + sum(abs(width));
            steps_tot = steps_tot + 1;
        end
    end
    t_check = toc(t_check);

    p_contained_val = 100 - (num_out/(num_out+num_in))*100;
    sizeI_val = size_acc / max(1,steps_tot);

    stats = struct('contain_pct', p_contained_val, ...
                   'sizeI', sizeI_val, ...
                   't_total', t_learn + t_check, ...
                   't_learn', t_learn, 't_check', t_check, 't_infer', NaN);
end

function v = getfielddef(S,f,def)
    if isobject(S), has=ismethod(S,'properties'); else, has=isstruct(S); end
    if has && (isprop(S,f) || (isstruct(S)&&isfield(S,f))), v = S.(f); else, v = def; end
end
