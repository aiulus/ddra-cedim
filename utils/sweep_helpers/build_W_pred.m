function Wpred = build_W_pred(gray_sys, ~, W_used)
    % Always return a zonotope when Gray has a disturbance port.
    nd = 0;
    if isprop(gray_sys,'nrOfDisturbances')
        try, nd = gray_sys.nrOfDisturbances; catch, nd = 0; end
    end
    if nd <= 0
        Wpred = [];                    % Gray truly has no disturbance port
        return;
    end

    try
        Wpred = normalizeWForGray(gray_sys, W_used);
        if isempty(Wpred)
            % Defensive: treat empty as zero set.
            Wpred = zonotope(zeros(nd,1));
        end
    catch
        % Any failure → zero disturbance set to keep reach consistent.
        warning('build_W_pred: falling back to zero zonotope (nd=%d).', nd);
        Wpred = zonotope(zeros(nd,1));
    end
end
