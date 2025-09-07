function [ctrain, cval, Tval, VAL] = gray_containment(configs, sys_cora, R0, U, C, pe, varargin)
% GRAY_CONTAINMENT
%   Containment on TRAIN and VAL using CORA's validateReach.
%   Linear systems: combine all test cases into one (as in CORA example).
%   Nonlinear systems: iterate per test case (no per-sample splitting).
%   Honors caller-provided plot_settings ([] → headless).

    % ---- args ----
    p = inputParser;
    addParameter(p,'check_contain', true);
    addParameter(p,'externalTS_train',[]);
    addParameter(p,'externalTS_val',[]);
    addParameter(p,'plot_settings',[]);   % [] suppresses plotting in validateReach
    parse(p, varargin{:});

    CC        = p.Results.check_contain;
    TS_train  = p.Results.externalTS_train;
    TS_val    = p.Results.externalTS_val;
    PS        = p.Results.plot_settings;

    % --- Containment tolerance parity (thread cfg.metrics.tol into configs.*.options.cs)
    tol = 1e-6;
    if isfield(C,'metrics') && isfield(C.metrics,'tol') && ~isempty(C.metrics.tol)
        tol = C.metrics.tol;
    end
    for ii = 1:numel(configs)
        if ~isfield(configs{ii},'options') || ~isstruct(configs{ii}.options)
            configs{ii}.options = struct();
        end
        if ~isfield(configs{ii}.options,'cs') || ~isstruct(configs{ii}.options.cs)
            configs{ii}.options.cs = struct();
        end
        configs{ii}.options.cs.robustnessMargin = tol;
    end


    % ---- build test suites if not provided ----
    optTS = ts_options_from_pe(C, pe, sys_cora);
    optTS.stateSet = R0;

    if isempty(TS_train)
        p_train = struct('R0',R0,'U',U,'tFinal', sys_cora.dt*(C.shared.n_k-1));
        TS_train = createTestSuite(sys_cora, p_train, C.shared.n_k, C.shared.n_m, C.shared.n_s, optTS);
    end
    if isempty(TS_val)
        p_val = struct('R0',R0,'U',U,'tFinal', sys_cora.dt*(C.shared.n_k_val-1));
        TS_val = createTestSuite(sys_cora, p_val, C.shared.n_k_val, C.shared.n_m_val, C.shared.n_s_val, optTS);
    end

    % ---- pick gray config index consistently (robust) ----
    names = strings(numel(configs),1);
    for ii = 1:numel(configs)
        try
            if isstruct(configs{ii}) && isfield(configs{ii},'name') && ~isempty(configs{ii}.name)
                names(ii) = string(configs{ii}.name);
            else
                names(ii) = "";
            end
        catch
            names(ii) = "";
        end
    end
    
    % prefer the first requested method if available; else fall back to "graySeq"; else 2
    want = "graySeq";
    if isstruct(C) && isfield(C,'gray') && isfield(C.gray,'methodsGray') && ~isempty(C.gray.methodsGray)
        want = string(C.gray.methodsGray(1));
        if numel(want) > 1, want = want(1); end
    end
    
    i_gray = find(names == want, 1);
    if isempty(i_gray)
        i_gray = find(names == "graySeq", 1);
        if isempty(i_gray)
            i_gray = min(2, numel(configs));
        end
    end


    % ---- TRAIN ----
    [num_out_tr, num_all_tr] = eval_suite(TS_train, configs, CC, PS);
    ctrain = 100*(1 - num_out_tr(i_gray)/max(1, num_all_tr));

    % ---- VAL (timed) ----
    t0 = tic;
    [num_out_val, num_all_val] = eval_suite(TS_val, configs, CC, PS);
    Tval = toc(t0);
    cval = 100*(1 - num_out_val(i_gray)/max(1, num_all_val));

    % ---- VAL pack (x0/u/y for first sample of each testcase) ----
    VAL = [];
    try
        VAL = pack_VAL_from_TS(TS_val);
    catch
        % leave VAL=[]
    end
end

% ---- Helpers ----
function VAL = pack_VAL_from_TS(TS)
    B = numel(TS);
    VAL = struct('x0',{cell(1,B)},'u',{cell(1,B)},'y',{cell(1,B)});
    for b = 1:B
        VAL.x0{b} = TS{b}.initialState;
        VAL.u{b}  = squeeze(TS{b}.u(:,:,1));
        VAL.y{b}  = squeeze(TS{b}.y(:,:,1));
    end
end

% ---- evaluator (mirrors CORA example usage) ----
function [num_out_vec, num_all] = eval_suite(TS, configs, CC, PS)
    num_out_vec = zeros(numel(configs),1);
    num_all = 0;
    for m = 1:numel(TS)
        if isempty(PS)
            [~, ev] = validateReach(TS{m}, configs, CC);
        else
            [~, ev] = validateReach(TS{m}, configs, CC, PS);
        end
        num_out_vec = num_out_vec + ev.num_out;
        nk = size(TS{m}.y,1); ns = size(TS{m}.y,3);
        num_all = num_all + nk * ns;
    end
end