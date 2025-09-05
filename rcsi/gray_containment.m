function [ctrain, cval, Tval, VAL] = gray_containment(configs, sys_cora, R0, U, C, pe, varargin)
% GRAY_CONTAINMENT
%   Containment on TRAIN and VAL using CORA's validateReach, with
%   per-sample testCases to avoid y_a/zero-input fallback.

    % ---- args ----
    p = inputParser;
    addParameter(p,'check_contain', true);
    addParameter(p,'externalTS_train',[]);
    addParameter(p,'externalTS_val',[]);
    addParameter(p,'plot_settings',[]);     % optional, passed to validateReach
    parse(p, varargin{:});
    CC        = p.Results.check_contain;
    TS_train  = p.Results.externalTS_train;
    TS_val    = p.Results.externalTS_val;
    PS        = p.Results.plot_settings;
    
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

    % ---- pick gray config index consistently ----
    want = "graySeq";
    i_gray = find(cellfun(@(cfg) isstruct(cfg) && isfield(cfg,'name') ...
                          && want==string(cfg.name), configs), 1);
    if isempty(i_gray), i_gray = min(2, numel(configs)); end

    % ---- accumulate helpers ----
    function [num_out_vec, num_all] = accum_suite(TS)
        num_out_vec = zeros(numel(configs),1);
        num_all = 0;
        for m = 1:numel(TS)
            tc = TS{m};
            sysClass = class(sys_cora);
            if size(tc.u,3) > 1
                TCs = split_tc_samples(tc, sysClass);
            else
                TCs = {tc};
            end
            for s = 1:numel(TCs)
                % delegate to CORA (strips cs, handles padding, plots if PS given)
                [~, ev] = validateReach(TCs{s}, configs, CC);        
                close all;
                num_out_vec = num_out_vec + ev.num_out;
                % denominator = #time steps × #samples (here 1)
                nk = size(TCs{s}.y,1);
                ns = size(TCs{s}.y,3);  % 1 for single-sample tc
                num_all = num_all + nk*ns;
            end
        end
    end

    % ---- TRAIN ----
    [num_out_tr, num_all_tr] = accum_suite(TS_train);
    ctrain = 100*(1 - num_out_tr(i_gray)/max(1, num_all_tr));

    % ---- VAL (timed) ----
    t0 = tic;
    [num_out_val, num_all_val] = accum_suite(TS_val);
    Tval = toc(t0);
    cval = 100*(1 - num_out_val(i_gray)/max(1, num_all_val));
    

    % ---- VAL pack (optional, for downstream plotting outside CORA) ----
    VAL = [];
    try
        VAL = pack_VAL_from_TS(TS_val); 
    catch
        % leave VAL = []
    end
end

% ------ Helpers ------

function TS1 = split_tc_samples(tc, sysClass)
    S = size(tc.u,3);
    TS1 = cell(1,S);
    for s = 1:S
        u_s = tc.u(:,:,s);   % (n_k × n_u)
        y_s = tc.y(:,:,s);   % (n_k × n_y)
        TS1{s} = testCase(y_s, u_s, tc.initialState, tc.sampleTime, sysClass);
    end
end

function VAL = pack_VAL_from_TS(TS)
% Minimal VAL pack: x0/u/y for first sample of each testcase
    B = numel(TS);
    VAL = struct('x0',{cell(1,B)},'u',{cell(1,B)},'y',{cell(1,B)});
    for b = 1:B
        VAL.x0{b} = TS{b}.initialState;
        VAL.u{b}  = squeeze(TS{b}.u(:,:,1));   % (n_k × n_u)
        VAL.y{b}  = squeeze(TS{b}.y(:,:,1));   % (n_k × n_y)
    end
end
