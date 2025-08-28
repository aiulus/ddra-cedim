function y = coerce_numeric(v)
% Converts table vars of mixed types (double/string/cellstr/cell/logical)
% to a numeric vector, robustly.

    if istable(v)
        v = table2array(v);
    end

    if isnumeric(v) || islogical(v)
        y = double(v);
        return;
    end

    if isstring(v) || ischar(v)
        y = str2double(string(v));
        return;
    end

    if iscell(v)
        y = cellfun(@(x) ...
            (isnumeric(x) || islogical(x)) .* double(x) + ...
            (ischar(x) || isstring(x))     .* str2double(string(x)) + ...
            (~(isnumeric(x)||islogical(x)||ischar(x)||isstring(x))) .* NaN, ...
            v, 'UniformOutput', true);
        return;
    end

    % Fallback
    try
        y = str2double(string(v));
    catch
        y = NaN(size(v));
    end
end
