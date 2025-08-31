function sizeI_gray = gray_infer_size_on_VAL(sys, TS, C, params, varargin)
%GRAY_INFER_SIZE_ON_VAL  Aggregated output-interval size over VAL.
% Usage:
%   sizeI_gray = gray_infer_size_on_VAL(sys, TS_val, C, params)
%   sizeI_gray = gray_infer_size_on_VAL(sys, TS_val, C, params, 'overrideW', W) % accepted, ignored
%
% Notes:
%   • Accepts extra name–value args to be API-compatible with callers.
%   • Computes sum_k || interval-width(R_y(k)) ||_1 across all k and testcases.
%   • Uses C.shared.options_reach (with cs stripped) for pure reach.

    % -------- parse optional args (accepted for compatibility) --------
    p = inputParser;
    addParameter(p, 'overrideW', []);       % currently unused for linearSysDT reach
    parse(p, varargin{:});

    % -------- normalize inputs --------
    if ~iscell(TS), TS = {TS}; end

    opts = C.shared.options_reach;
    if isfield(opts, 'cs'), opts = rmfield(opts, 'cs'); end

    sizeI_gray = 0;

    for m = 1:numel(TS)
        tc = TS{m};

        % Build per-case params: deterministic input u and shifted R0
        params_m = params;
        params_m.R0     = params.R0 + tc.initialState;
        Ui = tc.u;            % (n_k × n_u) or (n_k × n_u × s)
        Ui = squeeze(Ui);     % drop singleton 3rd dim if present
        if size(Ui,1) < size(Ui,2)   % (m×n_k) -> transpose to (n_k×m)
            Ui = Ui.';
        end
        params_m.u      = Ui';
        params_m.tFinal = sys.dt * (size(params_m.u,2)-1);

        % Reach in output space
        R = reach(sys, params_m, opts);
        Ysets = R.timePoint.set;     % cell{1..K} of output sets

        % Aggregate interval widths
        for k = 1:numel(Ysets)
            Ik = interval(Ysets{k});
            sizeI_gray = sizeI_gray + sum(abs(Ik.sup(:) - Ik.inf(:)));
        end
    end
end
