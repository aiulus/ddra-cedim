function [T, PS] = load_augmented(results_dir)
    T  = readtable(fullfile(results_dir, 'summary.csv'));
    PS = readtable(fullfile(results_dir, 'summary_perstep.csv'));

    T = coerce_bool_cols(T, {'ddra_ridge','use_noise'});
    PS = coerce_bool_cols(PS, {});  % add names if you introduce any later
end

function T = coerce_bool_cols(T, names)
    for i = 1:numel(names)
        nm = names{i};
        if ~ismember(nm, T.Properties.VariableNames), continue; end
        col = T.(nm);
        if islogical(col), continue; end
        if iscell(col),        col = string(col); end
        if iscategorical(col), col = string(col); end
        if isstring(col) || ischar(col)
            s = lower(string(col));
            T.(nm) = ismember(s, ["true","1","yes","y"]);
        elseif isnumeric(col)
            T.(nm) = col ~= 0;
        else
            T.(nm) = false(height(T),1);
        end
        T.(nm) = logical(T.(nm));
    end
end
