function optTS = ts_options_from_pe(C, pe, ~)
%TS_OPTIONS_FROM_PE  Build testSuite options consistent with additional PE spec.

    % -------- defaults --------
    optTS = struct();

    % p_extr default
    if isstruct(C) && isfield(C,'shared') && isfield(C.shared,'p_extr') ...
            && ~isempty(C.shared.p_extr)
        optTS.p_extr = C.shared.p_extr;
    else
        optTS.p_extr = 0.3;
    end

    % inputCurve mapping (CORA expects "rand" / "sin" / ...)
    mode = 'rand';
    if nargin >= 2 && ~isempty(pe) && isstruct(pe) && isfield(pe,'mode') && ~isempty(pe.mode)
        m = lower(char(pe.mode));
        switch m
            case {'rand','randn'}
                mode = 'rand';
            case {'sin','sinwave','sine'}
                mode = 'sin';
            otherwise
                mode = 'rand';
        end
    end
    optTS.inputCurve = string(mode);

    % optional: amplitude passthrough
    if isstruct(pe) && isfield(pe,'strength') && ~isempty(pe.strength)
        optTS.inputStrength = pe.strength;
    end

    % deterministic seeding parity with DDRA
    if isstruct(C) && isfield(C,'shared') && isfield(C.shared,'seed') && ~isempty(C.shared.seed)
        optTS.rngSeed = C.shared.seed;   % your createTestSuite should read this
    end
end
