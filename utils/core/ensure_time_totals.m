function T = ensure_time_totals(T)
    % create totals if missing, coercing types as needed
    if ~ismember('t_ddra_total', T.Properties.VariableNames)
        L = get_numcol(T,'t_ddra_learn');
        C = get_numcol(T,'t_ddra_check');
        I = get_numcol(T,'t_ddra_infer');
        T.t_ddra_total = L + C + I;
    end
    if ~ismember('t_gray_total', T.Properties.VariableNames)
        Lg = get_numcol(T,'t_gray_learn');
        Vg = get_numcol(T,'t_gray_val');
        Ig = get_numcol(T,'t_gray_infer');
        T.t_gray_total = Lg + Vg + Ig;
    end
end

function v = get_numcol(T, name)
    % returns numeric column or zeros if missing
    if ~ismember(name, T.Properties.VariableNames)
        v = zeros(height(T),1);
        return;
    end
    col = T.(name);
    v = coerce_numeric(col);  % your existing helper; see below fallback
end

% (fallback in case coerce_numeric isn't already in scope)
function v = coerce_numeric(col)
    if isnumeric(col)
        v = col;
    elseif iscell(col)
        v = cellfun(@(x) str2double(string(x)), col);
    elseif isstring(col) || ischar(col) || isa(col,'categorical')
        v = str2double(string(col));
    else
        v = double(col);
    end
end
