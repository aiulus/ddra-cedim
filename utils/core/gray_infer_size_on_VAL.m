function [sizeI, wid_k] = gray_infer_size_on_VAL(sys, TS_val, C, params_in, varargin)
%GRAY_INFER_SIZE_ON_VAL  Output-interval size proxy for Gray on VAL points.
% Returns:
%   sizeI : scalar = mean interval width over ALL (blocks × time)
%   wid_k : n_k-by-1 = mean interval width per time step (avg across blocks)

    p = inputParser;
    addParameter(p,'overrideW',[]);     % if provided, use this W (any space)
    parse(p, varargin{:});
    W_override = p.Results.overrideW;

    % base options: remove 'cs' before calling reach()
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

    for m = 1:numel(TS_val)
        tc = TS_val{m};
        nk = size(tc.y,1);

        % --- representative initial state for this testcase (nx×1)
        x0bar = mean(tc.initialState, 3);           % average over samples
        % (alternative: tc.initialState(:,:,1))

        % --- build params for this testcase
        params = struct();
        params.R0     = params_in.R0 + x0bar;       % shift by vector (nx×1)
        params.u      = tc.u.';                     % measured input, (nu×nk)
        params.tFinal = sys.dt * (nk-1);

        % input padding if model expects more channels
        du = sys.nrOfInputs - size(params.u,1);
        if du > 0
            params.u = [params.u; zeros(du, size(params.u,2))];
        elseif du < 0
            params.u = params.u(1:sys.nrOfInputs, :);  % defensive
        end

        % --- disturbance set W (required if nrOfDisturbances>0)
        if sys.nrOfDisturbances > 0
            if ~isempty(W_override)
                params.W = coerceWToSys(sys, W_override);
            elseif isfield(params_in,'W') && ~isempty(params_in.W)
                params.W = coerceWToSys(sys, params_in.W);
            else
                params.W = zonotope(zeros(sys.nrOfDisturbances,1));
            end
        end

        % --- reach in state space
        R = reach(sys, params, options);            % CORA object
        R = R.timePoint.set;                        % cell array of contSets

        % --- per-step sizes and accumulation (avg across blocks only)
        for k = 1:nk
            Xk = R{k};
            if isempty(Xk) || ~isa(Xk,'contSet'), continue; end
            % output map (or identity if C missing)
            try
                Yk = linearMap(Xk, sys.C);
            catch
                Yk = Xk;
            end

            % interval width (sum over outputs)
            Ik = interval(Yk);
            
            % CORA-compatible: width = sup - inf
            try
                lo = infimum(Ik);
                hi = supremum(Ik);
            catch
                % very old CORA builds store fields directly
                if isstruct(Ik) && isfield(Ik,'inf') && isfield(Ik,'sup')
                    lo = Ik.inf;
                    hi = Ik.sup;
                else
                    error('interval() returned an unexpected type; cannot extract bounds.');
                end
            end
            
            wk = sum(max(hi - lo, 0));   % guard tiny negatives from numerics

            wid_sums(k)   = wid_sums(k)   + wk;
            wid_counts(k) = wid_counts(k) + 1;

            % global scalar proxy uses the same widths
            acc   = acc   + wk;
            count = count + 1;
        end
    end

    % outputs
    if count==0
        sizeI = NaN;
    else
        sizeI = acc / count;
    end
    wid_k = wid_sums ./ max(1, wid_counts);
end

% --- Helpers ---
function v = getfielddef(S,f,d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end
