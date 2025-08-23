function optTS = ts_options_from_pe(C, pe)
    optTS = struct('p_extr', C.shared.p_extr);
    if isfield(C.shared,'testSuite_mode') && C.shared.testSuite_mode=="ddra_like"
        optTS.inputCurve = "randn";
        optTS.contInput  = false;
        % Do NOT set stateSet here. Caller will set it to R0 if desired.
    end
    if isfield(pe,'mode') && pe.mode=="multisine"
        optTS.inputCurve = "sinWave";   % proxy for PE
    end
end