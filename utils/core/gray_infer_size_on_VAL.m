function sizeI = gray_infer_size_on_VAL(sys, TS_val, C, params_in, varargin)
%GRAY_INFER_SIZE_ON_VAL  Output-interval size proxy for Gray on VAL points.
% Uses measured u (per testcase), injects params.W (zero if noise disabled),
% strips optimizer block from options, and maps to output via linearMap.

    p = inputParser;
    addParameter(p,'overrideW',[]);     % if provided, use this W
    parse(p, varargin{:});
    W_override = p.Results.overrideW;

    % base options: remove 'cs' before calling reach()
    options = getfielddef(C.shared,'options_reach', struct());
    if isfield(options,'cs'); options = rmfield(options,'cs'); end

    acc = 0; count = 0;

    for m = 1:numel(TS_val)
        tc = TS_val{m};
        nk = size(tc.y,1);

        % --- build params for this testcase
        params = struct();
        params.R0     = params_in.R0 + tc.initialState;
        params.u      = tc.u';                              % measured input
        params.tFinal = sys.dt * (nk-1);

        % input padding if model expects more channels
        du = sys.nrOfInputs - size(params.u,1);
        if du ~= 0
            params.u = cat(1, params.u, zeros(du, size(params.u,2), size(params.u,3)));
        end

        % --- disturbance set W (required if nrOfDisturbances>0)
        if sys.nrOfDisturbances > 0
            if ~isempty(W_override)
                params.W = W_override;
            elseif isfield(params_in,'W') && ~isempty(params_in.W)
                params.W = params_in.W;
            else
                % zero-disturbance zonotope with correct dimension
                params.W = zonotope(zeros(sys.nrOfDisturbances,1));
            end
        end

        % --- reach in state space
        R = reach(sys, params, options);   % CORA object
        R = R.timePoint.set;               % cell array of contSets

        % --- map to outputs and accumulate interval width
        for k = 1:nk
            Xk = R{k};
            if isempty(Xk) || ~isa(Xk,'contSet'); continue; end
            try
                Yk = linearMap(Xk, sys.C);   % CORA: set first, matrix second
            catch
                % fallback for identity-output or missing C
                Yk = Xk;
            end
            if isa(Yk,'zonotope')
                acc = acc + sum(abs(generators(Yk)), 'all');
                count = count + 1;
            else
                % best-effort fallback
                try
                    Z = zonotope(Yk);
                    acc = acc + sum(abs(generators(Z)), 'all');
                    count = count + 1;
                catch
                    % skip if unhandled set type
                end
            end
        end
    end

    sizeI = acc / max(1, count);
end

% --- tiny util (kept local to avoid cross-file deps) ---
function v = getfielddef(S,f,d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end
