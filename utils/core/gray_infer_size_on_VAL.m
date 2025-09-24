function [sizeI, wid_k] = gray_infer_size_on_VAL(sys, TS_val, C, params_in, varargin)
%GRAY_INFER_SIZE_ON_VAL  Output-interval size proxy for Gray on VAL points.
% Returns:
%   sizeI : scalar = average across (blocks × time) of the per-step aggregator
%   wid_k : n_k-by-1 = per-step aggregator averaged across blocks
%
% Aggregator over outputs:
%   - 'sum'  : sum of interval widths across outputs (legacy default)
%   - 'mean' : mean of interval widths across outputs (matches DDRA tests)

    % ------------ args / options ------------
    p = inputParser;
    addParameter(p, 'overrideW', []);
    addParameter(p, 'width_agg', []);   % NEW: 'sum' | 'mean'
    parse(p, varargin{:});
    W_override = p.Results.overrideW;

    % decide aggregator (varargin > C.metrics.width_agg > 'sum')
    width_agg = lower(local_getfielddef(C, 'metrics', struct()));
    width_agg = lower(local_getfielddef(width_agg, 'width_agg', ''));
    if ~isempty(p.Results.width_agg), width_agg = lower(string(p.Results.width_agg)); end
    if width_agg == "", width_agg = "sum"; end           % backward-compatible default
    if ~ismember(width_agg, ["sum","mean"])
        error('gray_infer_size_on_VAL: width_agg must be "sum" or "mean".');
    end

    options = local_getfielddef(local_getfielddef(C,'shared',struct()),'options_reach', struct());
    if isfield(options,'cs'); options = rmfield(options,'cs'); end

    % allow nk from config, but tolerate testcases with different nk
    nk_global = local_getfielddef(local_getfielddef(C,'shared',struct()),'n_k_val',0);
    if nk_global<=0
        nk_global = max(cellfun(@(tc) size(tc.y,1), TS_val));
    end
    wid_sums   = zeros(nk_global,1);
    wid_counts = zeros(nk_global,1);

    acc = 0; 
    count = 0;

    % --- zero-center R0 once (guarded)
    nx = size(sys.A,1);
    if isfield(params_in,'R0') && ~isempty(params_in.R0)
        R0c = params_in.R0;
        try
            R0c = zonotope([zeros(size(center(R0c),1),1), generators(R0c)]);
        catch
            % keep as-is if not a zonotope
        end
    else
        % default: exact point at origin in state space
        R0c = zonotope(zeros(nx,1));
    end


    for m = 1:numel(TS_val)
        tc = TS_val{m};
        nk = size(tc.y,1);

        % representative initial state for this testcase (nx×1)
        x0bar = mean(tc.initialState, 3);   % or tc.initialState(:,:,1)

        % --- build params for this testcase
        params = struct();
        params.R0     = R0c + x0bar;        % shift by vector (nx×1)
        params.u      = tc.u.';             % (nu×nk)
        params.tFinal = sys.dt * (nk-1);

        % input padding if model expects more channels
        du = sys.nrOfInputs - size(params.u,1);
        if du > 0
            params.u = [params.u; zeros(du, size(params.u,2))];
        elseif du < 0
            params.u = params.u(1:sys.nrOfInputs, :);
        end

        % disturbance set W
        if sys.nrOfDisturbances > 0
            if ~isempty(W_override) && ~isa(W_override,'emptySet')
                params.W = coerceWToSys(sys, W_override);
            elseif isfield(params_in,'W') && ~isempty(params_in.W)
                params.W = coerceWToSys(sys, params_in.W);
            else
                params.W = zonotope(zeros(sys.nrOfDisturbances,1));
            end
        end

        % --- reach & widths
        R = reach(sys, params, options);
        R = R.timePoint.set;

        for k = 1:nk
            Xk = R{k};
            if isa(Xk,'emptyset') || ~isa(Xk,'contSet'), continue; end

            % map to outputs; if mapping unsupported, use state set (D*u added below)
            try, Yk = linearMap(Xk, sys.C); catch, Yk = Xk; end

            % feedthrough (does not change interval width but keeps Y_k explicit)
            uk = params.u(:,k);
            if isfield(sys,'D') && ~isempty(sys.D)
                Yk = Yk + sys.D * uk;
            end

            Ik = interval(Yk);
            try, lo = infimum(Ik); hi = supremum(Ik); catch, lo = Ik.inf; hi = Ik.sup; end
            wvec = max(hi - lo, 0);            % width per output dimension

            % --- aggregator over outputs (NEW)
            switch width_agg
                case "sum"
                    wk = sum(double(wvec));
                case "mean"
                    wk = mean(double(wvec), 'omitnan');
            end

            wid_sums(k)   = wid_sums(k)   + wk;
            wid_counts(k) = wid_counts(k) + 1;
            acc = acc + wk; count = count + 1;
        end
    end

    sizeI = (count==0) * NaN + (count>0) * (acc / count);
    wid_k = wid_sums ./ max(1, wid_counts);
end

% --- Helpers (local copies to avoid collisions) ---
function v = local_getfielddef(S,f,d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end
