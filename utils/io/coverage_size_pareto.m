function fig = coverage_size_pareto(T)
    % Try to find coverage & size columns
    v = T.Properties.VariableNames;
    c_g = v(contains(v,'cval') & contains(v,'gray'));   % e.g., cval_gray
    c_d = v(contains(v,'cval') & contains(v,'ddra'));   % e.g., cval_ddra
    s_g = v(contains(v,'sizeI') & contains(v,'gray'));  % e.g., sizeI_gray
    s_d = v(contains(v,'sizeI') & contains(v,'ddra'));  % e.g., sizeI_ddra

    assert(~isempty(c_g) && ~isempty(c_d) && ~isempty(s_g) && ~isempty(s_d), 'Missing coverage/size columns.');

    % Fairness highlight (ridge rows)
    ridge = false(height(T),1);
    if ismember('ddra_ridge', T.Properties.VariableNames)
        ridge = logical(T.ddra_ridge);
    end

    fig = figure('Color','w'); hold on; grid on;
    scatter(T.(s_g{1}), T.(c_g{1}), 25, 'filled', 'MarkerFaceAlpha',.7, 'DisplayName','Gray');
    scatter(T.(s_d{1}), T.(c_d{1}), 25, 'filled', 'MarkerFaceAlpha',.7, 'DisplayName','DDRA');

    % Outline ridge-inflated rows if present
    if any(ridge)
        idx = ridge;
        scatter(T.(s_d{1})(idx), T.(c_d{1})(idx), 60, 'o', 'LineWidth', 1.5, 'DisplayName','DDRA (ridge)');
    end

    xlabel('mean interval width (size proxy)'); ylabel('coverage (%)');
    title('Coverage–size tradeoff (VAL)'); legend('Location','southeast'); axis tight;
end
