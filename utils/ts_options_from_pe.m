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

    % inputCurve mapping
    mode = 'rand';  % CORA expects "rand" / "sin" / etc.
    if nargin >= 2 && ~isempty(pe) && isstruct(pe) && isfield(pe,'mode') && ~isempty(pe.mode)
        m = lower(char(pe.mode));
        switch m
            case {'randn','rand'}
                mode = 'rand';
            case {'sin','sinwave','sine'}
                mode = 'sin';
            % add other mappings here if you introduce new PE modes:
            % case {'step','steps'},  mode = 'steps';
            otherwise
                mode = 'rand';
        end
    end
    optTS.inputCurve = string(mode);

    % if isfield(pe,'strength'), optTS.inputStrength = pe.strength; end
end
