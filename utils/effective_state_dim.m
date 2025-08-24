function nx = effective_state_dim(sys, options)
    if sys.nrOfDims > 0
        nx = sys.nrOfDims;
        return;
    end
    % ARX/NARX case: outputs * past steps
    p = [];
    try
        p = sys.n_p;                
    catch
        try p = sys.nrOfPastSteps; end
    end
    if isempty(p)
        % fall back to your approximator config if present, else 1
        if isfield(options,'approx') && isfield(options.approx,'p') && ~isempty(options.approx.p)
            p = options.approx.p;
        else
            p = 1;
        end
    end
    nx = sys.nrOfOutputs * p;
end
