function sizeI = gray_infer_size_on_VAL(sys, TS_val, C, params_in, varargin)
% Output-interval size proxy for Gray on the exact VAL points.
% Uses measured inputs; always sets params.W if disturbances exist.

    p = inputParser;
    addParameter(p,'overrideW',[]);
    parse(p, varargin{:});
    W_override = p.Results.overrideW;

    % reachability options (strip optimizer block)
    options = getfielddef(C.shared,'options_reach', struct());
    if isfield(options,'cs'); options = rmfield(options,'cs'); end

    acc = 0; count = 0;

    for m = 1:numel(TS_val)
        tc = TS_val{m};

        % split multi-sample testCase into single-sample cases
        if size(tc.u,3) > 1
            TCs = split_tc_samples(tc);
        else
            TCs = {tc};
        end

        for s = 1:numel(TCs)
            tcs = TCs{s};
            nk  = size(tcs.y,1);

            % --- params: R0 shift and measured inputs for this sample
            params.R0     = params_in.R0 + tcs.initialState;   % (nx×1)
            params.u      = tcs.u';                            % (nu×nk)
            params.tFinal = sys.dt * (nk-1);

            % pad inputs if needed
            du = sys.nrOfInputs - size(params.u,1);
            if du ~= 0
                params.u = cat(1, params.u, zeros(du, size(params.u,2), size(params.u,3)));
            end

            % --- disturbances (required if nrOfDisturbances>0)
            if sys.nrOfDisturbances > 0
                if ~isempty(W_override)
                    params.W = W_override;
                elseif isfield(params_in,'W') && ~isempty(params_in.W)
                    params.W = params_in.W;
                else
                    params.W = zonotope(zeros(sys.nrOfDisturbances,1)); % zero but present
                end
            end

            % --- reach in state space
            R = reach(sys, params, options);      % CORA object
            R = R.timePoint.set;

            % --- map to outputs & accumulate interval width
            for k = 1:nk
                Xk = R{k};
                if isempty(Xk) || ~isa(Xk,'contSet'); continue; end
                Yk = map_to_output(Xk, sys.C);     % robust mapping
                if isa(Yk,'zonotope')
                    acc   = acc + sum(abs(generators(Yk)),'all');
                    count = count + 1;
                end
            end
        end
    end

    sizeI = acc / max(1,count);
end

% ----- helpers -----
function Y = map_to_output(Xset, C)
    try
        Y = linearMap(Xset, C);   % CORA >= certain versions
    catch
        c = center(Xset); G = generators(Xset);
        Y = zonotope(C*c, C*G);
    end
end

function TS1 = split_tc_samples(tc)
    S = size(tc.u,3); TS1 = cell(1,S);
    for s = 1:S
        u_s = tc.u(:,:,s);            % (nk×nu)
        y_s = tc.y(:,:,s);            % (nk×ny)
        TS1{s} = testCase(y_s, u_s, tc.initialState(:,:,s), tc.sampleTime, class(tc));
    end
end

function v = getfielddef(S,f,d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end
