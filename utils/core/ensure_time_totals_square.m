function T = ensure_time_totals_square(T)
    addcol = @(tbl,name,val) iff(~ismember(name, tbl.Properties.VariableNames), ...
        @() addvars(tbl, val, 'NewVariableNames', name), ...
        @() tbl);
    iff = @(cond,a,b) (cond)*a() + (~cond)*b(); %#ok<NASGU> (trick to keep inline)

    if ~ismember('t_black_total', T.Properties.VariableNames)
        T.t_black_total = zeros(height(T),1);
    end
    if ~ismember('t_ddra_total', T.Properties.VariableNames)
        T.t_ddra_total = zeros(height(T),1);
    end
    if all(ismember({'t_black_learn','t_black_val','t_black_infer'}, T.Properties.VariableNames))
        T.t_black_total = T.t_black_learn + T.t_black_val + T.t_black_infer;
    end
    if all(ismember({'t_ddra_build','t_ddra_val','t_ddra_reach_avg'}, T.Properties.VariableNames))
        T.t_ddra_total = T.t_ddra_build + T.t_ddra_val + T.t_ddra_reach_avg;
    end
end