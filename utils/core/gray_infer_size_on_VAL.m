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
                params.W = coerceWToSys(sys, W_override);
            elseif isfield(params_in,'W') && ~isempty(params_in.W)
                params.W = coerceWToSys(sys, params_in.W);
            else
                params.W = zonotope(zeros(sys.nrOfDisturbances,1));
            end
        end

        % --- reach in state space
        R = reach(sys, params, options);   % CORA object
        R = R.timePoint.set;               % cell array of contSets

        % --- map to outputs and accumulate volume
        vols = [];
        for k = 1:nk
            Xk = R{k};
            if isempty(Xk) || ~isa(Xk,'contSet'), continue; end
            % output map (or identity if C is missing)
            try
                Yk = linearMap(Xk, sys.C);
            catch
                Yk = Xk;
            end
            vols(end+1,1) = norm(Yk);
        end
        acc   = acc   + sum(vols);
        count = count + numel(vols);
        if count == 0
            sizeI = NaN;   
        else
            sizeI = acc / count;
        end 
    end
end

% --- Helpers ---
function v = getfielddef(S,f,d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end

function Wd = coerceWToSys(sys_, W_in)
    % Return [] for no disturbances, or a zonotope in R^{nrOfDisturbances}
    nw = sys_.nrOfDisturbances;
    if nw == 0 || isempty(W_in)
        Wd = []; return;
    end
    if ~isa(W_in,'zonotope')
        warning('overrideW is not a zonotope; dropping to zero-disturbance.'); 
        Wd = zonotope(zeros(nw,1)); return;
    end
    d_in = size(center(W_in),1);
    if d_in == nw
        Wd = W_in; return;  % already correct dimension
    end
    % Try to map state-noise to disturbance-noise via left pseudo-inverse of E
    E = getfielddef(sys_,'E',[]);
    if ~isempty(E) && size(E,1) == d_in && size(E,2) == nw
        % Conservative preimage via pinv: Wd := pinv(E) * Ws
        % NOTE: E*Wd == Proj_{col(E)}(Ws); this is an approximation for size metric.
        Wd = linearMap(W_in, pinv(E));
        % If you want an explicit safety margin, add a tiny ε-box here:
        % epsBox = zonotope(zeros(nw,1), 1e-12*eye(nw));
        % Wd = Wd + epsBox;
        return;
    end
    % Fallback: zero disturbance with correct dimension
    warning(['overrideW dim=%d incompatible with sys.nrOfDisturbances=%d; ' ...
             'falling back to zero-disturbance for Gray size metric.'], d_in, nw);
    Wd = zonotope(zeros(nw,1));
end