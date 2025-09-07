function fig = boxplot_directional_error(T)
    % Look for E2 columns (median / p90), both methods if available
    V = T.Properties.VariableNames;
    e2g_med = V(contains(V,'E2') & contains(V,'gray') & contains(V,'med'));
    e2d_med = V(contains(V,'E2') & contains(V,'ddra') & contains(V,'med'));
    assert(~isempty(e2g_med) && ~isempty(e2d_med), 'E2 median columns for Gray/DDRA not found.');

    fig = figure('Color','w'); 
    boxchart(categorical(repelem(["Gray","DDRA"], height(T))'), [T.(e2g_med{1}); T.(e2d_med{1})]);
    ylabel('Directional support error (median across k,d)'); title('Shape-aware conservatism (E2)');
    grid on;
end
