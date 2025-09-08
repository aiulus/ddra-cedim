function [sizeI, wid_k] = gray_infer_size_on_VAL(sys, TS_val, C, params_in, varargin)
%GRAY_INFER_SIZE_ON_VAL  Output-interval size proxy for Gray on VAL points.
% Returns:
%   sizeI : scalar = mean interval width over ALL (blocks × time)
%   wid_k : n_k-by-1 = mean interval width per time step (avg across blocks)

    p = inputParser;
    addParameter(p,'overrideW',[]);
    parse(p, varargin{:});
    W_override = p.Results.overrideW;
    options = getfielddef(C.shared,'options_reach', struct());
    if isfield(options,'cs'); options = rmfield(options,'cs'); end

    % allow nk from config, but tolerate testcases with different nk
    nk_global = getfielddef(getfielddef(C,'shared',struct()),'n_k_val',0);
    if nk_global<=0
        nk_global = max(cellfun(@(tc) size(tc.y,1), TS_val));
    end
    wid_sums   = zeros(nk_global,1);
    wid_counts = zeros(nk_global,1);

    acc = 0; 
    count = 0;

    % --- zero-center R0 once
    R0c = params_in.R0;
    try
        R0c = zonotope([zeros(size(center(R0c),1),1), generators(R0c)]);
    catch
        % if R0c is already zero-centered or not a zonotope, keep as-is
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
            if ~isa(W_override, 'emptySet')
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
            if isa(Xk, 'emptyset') || ~isa(Xk,'contSet'), continue; end
            try, Yk = linearMap(Xk, sys.C); catch, Yk = Xk; end

            Ik = interval(Yk);
            try
                lo = infimum(Ik); hi = supremum(Ik);
            catch
                lo = Ik.inf; hi = Ik.sup;  % very old CORA
            end
            wk = sum(max(hi - lo, 0));

            wid_sums(k)   = wid_sums(k)   + wk;
            wid_counts(k) = wid_counts(k) + 1;

            acc   = acc   + wk;
            count = count + 1;
        end
    end

    sizeI = (count==0) * NaN + (count>0) * (acc / count);
    wid_k = wid_sums ./ max(1, wid_counts);
end

% --- Helpers ---
function v = getfielddef(S,f,d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end
