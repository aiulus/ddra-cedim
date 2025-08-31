function sizeI_excess = gray_infer_excess_on_VAL(sys, TS, C)
    % Excess size over the data spread (sum over all steps & testcases)
    opts = C.shared.options_reach;
    if isfield(opts,'cs'), opts = rmfield(opts,'cs'); end
    sizeI_excess = 0;
    for m = 1:numel(TS)
        params = struct('R0', C.shared.R0_stub, 'U', C.shared.U_stub); 
        params.R0 = TS{m}.initialState;   % point
        params.u  = TS{m}.u';             % (m × n_k)
        params.tFinal = sys.dt*(size(params.u,2)-1);

        R = reach(sys, params, opts);
        sets = R.timePoint.set;

        Y = TS{m}.y;
        hasY = ~isempty(Y) && ~all(isnan(Y(:)));
        for k = 1:size(sets,2)
            I = interval(sets{k});
            widthR = sum(abs(I.sup - I.inf));

            widthData = 0;
            if hasY
                yk = squeeze(Y(k,:,:)); if isvector(yk), yk = yk(:); end
                ymin = min(yk,[],2); ymax = max(yk,[],2);
                widthData = sum(abs(ymax - ymin));
            end
            sizeI_excess = sizeI_excess + max(0, widthR - widthData);
        end
    end
end
