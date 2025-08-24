function optTS = ts_options_from_pe(C, pe, sys)
% Build CORA testSuite options from our PE spec, mirroring DDRA's shapes
    n_u = sys.nrOfInputs;
    n_k = C.shared.n_k;
    dt  = sys.dt;
    % (n_m not directly used by createTestSuite options, but kept for symmetry)
    % n_m = C.shared.n_m; 

    [optTS, ~] = pe_param_synth(pe, n_u, n_k, dt);

    % Match global knobs
    optTS.p_extr    = getfielddef(C.shared,'p_extr', 0.3);
    optTS.contInput = true;  % CORA default 
end

function v = getfielddef(S, f, d)
    if isstruct(S) && isfield(S,f), v = S.(f); else, v = d; end
end
