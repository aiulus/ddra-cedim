function [ctrain, cval, Tval] = gray_containment(configs, sys_cora, R0, U, C, pe, varargin)
    p = inputParser;
    addParameter(p,'check_contain', true);
    addParameter(p,'externalTS_train',[]);
    addParameter(p,'externalTS_val',[]);
    parse(p, varargin{:});
    TS_train_ext = p.Results.externalTS_train;
    TS_val_ext   = p.Results.externalTS_val;
    CC = p.Results.check_contain;

    optTS = ts_options_from_pe(C, pe, sys_cora);
    optTS.stateSet = R0;

    %params_true = struct('R0',R0,'U',U,'tFinal', sys_cora.dt*(C.shared.n_k-1));
    %params_true.testSuite = createTestSuite(sys_cora, params_true, C.shared.n_k, C.shared.n_m, C.shared.n_s, optTS);
    %params_true.tFinal = sys_cora.dt*(C.shared.n_k_val-1);
    %TS_val = createTestSuite(sys_cora, params_true, C.shared.n_k_val, C.shared.n_m_val, C.shared.n_s_val, optTS);

    % ------------ NEW ------------
    % TRAINING
    if ~isempty(TS_train_ext)
        TS_train = TS_train_ext;
    else
        params_true = struct('R0',R0,'U',U,'tFinal', sys_cora.dt*(C.shared.n_k-1));
        TS_train = createTestSuite(sys_cora, params_true, C.shared.n_k, C.shared.n_m, C.shared.n_s, optTS);
    end
    
    % VALIDATION
    if ~isempty(TS_val_ext)
        TS_val = TS_val_ext;
    else
        params_true = struct('R0',R0,'U',U,'tFinal', sys_cora.dt*(C.shared.n_k_val-1));
        TS_val = createTestSuite(sys_cora, params_true, C.shared.n_k_val, C.shared.n_m_val, C.shared.n_s_val, optTS);
    end

    num_out = zeros(numel(configs),1); t0 = tic;
    for m=1:length(params_true.testSuite)
        [~, eval_id] = validateReach(params_true.testSuite{m}, configs, CC);
        num_out = num_out + eval_id.num_out;
    end
    num_all_id = length(params_true.testSuite)*C.shared.n_k*size(params_true.testSuite{1}.y,3);
    ctrain = 100*(1 - num_out(2)/max(1,num_all_id));

    num_out_val = zeros(numel(configs),1);
    for m=1:length(TS_val)
        [~, eval_val] = validateReach(TS_val{m}, configs, CC);
        num_out_val = num_out_val + eval_val.num_out;
    end
    num_all_val = length(TS_val)*C.shared.n_k_val*size(TS_val{1}.y,3);
    cval = 100*(1 - num_out_val(2)/max(1,num_all_val));
    Tval = toc(t0);
end