function sizeI = gray_infer_size_on_VAL(sys_gray, TS_val, C, params_gray_opt)
%GRAY_INFER_SIZE_ON_VAL  Aggregate interval-size proxy for Gray on VAL points.
% Usage:
%   sizeI = gray_infer_size_on_VAL(configs{2}.sys, TS_val, C, configs{2}.params)
%   sizeI = gray_infer_size_on_VAL(configs{2}.sys, TS_val, C)   % fallback (no R0/U uncertainty)
%
% Notes:
% - Uses the *same validation (x0,u)* sequences as DDRA.
% - Strips options.cs before reach().
% - Size metric = sum_k sum |Iv.sup - Iv.inf|, aggregated across all VAL cases.
%
% Recommended: pass params_gray_opt = configs{2}.params so R0/U match Gray ID.

    % Reach options: remove conformance fields to match plain reach()
    if isfield(C, 'shared') && isfield(C.shared, 'options_reach')
        options = C.shared.options_reach;
    else
        options = struct();
    end
    if isfield(options,'cs'), options = rmfield(options,'cs'); end

    % Base params (prefer identified params if provided)
    has_gray_params = (nargin >= 4) && ~isempty(params_gray_opt);
    if has_gray_params
        params_base = params_gray_opt;  % contains R0, U (and maybe more)
    else
        % Fallback: no uncertainty (point R0/U) – still well-defined
        params_base = struct();
        params_base.R0 = zonotope(zeros(sys_gray.nrOfDims,1));
        params_base.U  = zonotope(zeros(sys_gray.nrOfInputs,1));
    end

    sizeI = 0;

    % Loop over VAL test cases
    for m = 1:numel(TS_val)
        % Build params for this case
        params_m = params_base;

        % Shift R0 center by the test's initial state (keep generators if present)
        if isfield(params_m,'R0') && isa(params_m.R0,'zonotope')
            params_m.R0 = params_m.R0 + TS_val{m}.initialState;
        else
            params_m.R0 = zonotope(TS_val{m}.initialState);
        end

        % TS_val{m}.u is (n_k × n_u × 1) -> need (n_u × n_k)
        u_time_major = TS_val{m}.u;                         % (n_k × n_u × 1)
        params_m.u   = reshape(permute(u_time_major,[2 1 3]), ...
                               sys_gray.nrOfInputs, []);    % (n_u × n_k)

        % Horizon
        n_k = size(params_m.u, 2);
        params_m.tFinal = sys_gray.dt * (n_k - 1);

        % Reach
        Rg = reach(sys_gray, params_m, options);
        sets_k = Rg.timePoint.set;  % cell {1..n_k}

        % Aggregate interval-width sizes
        for k = 1:numel(sets_k)
            Zk = sets_k{k};
            if ~isa(Zk,'contSet'), Zk = zonotope(Zk); end
            Iv = interval(Zk);
            sizeI = sizeI + sum(abs(Iv.sup(:) - Iv.inf(:)));
        end
    end
end
