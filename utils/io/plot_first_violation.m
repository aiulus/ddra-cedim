function fig = plot_first_violation(T)
    % Try likely column names
    cand = T.Properties.VariableNames(contains(T.Properties.VariableNames, 'first_violation', 'IgnoreCase',true) | ...
                                      contains(T.Properties.VariableNames, 'FV', 'IgnoreCase',true));
    assert(~isempty(cand), 'First-violation metric (E6) not found.');
    col = cand{1};

    fig = figure('Color','w'); 
    histogram(T.(col), 'BinMethod','integers'); grid on;
    xlabel('First violation step'); ylabel('#rows'); title('Distribution of first-violation step (VAL)');
end
