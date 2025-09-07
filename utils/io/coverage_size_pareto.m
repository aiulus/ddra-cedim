function coverage_size_pareto(T)
% COVERAGE_SIZE_PARETO  Scatter coverage (%) vs size proxy for DDRA and Gray.
% Expects columns: cval_ddra, sizeI_ddra, cval_gray, sizeI_gray.
% Optional: ddra_ridge (logical) to style DDRA points.

    need = {'cval_ddra','sizeI_ddra','cval_gray','sizeI_gray'};
    assert(all(ismember(need, T.Properties.VariableNames)), ...
        'coverage_size_pareto: table must contain %s.', strjoin(need,', '));

    % Build clean masks per method (Nx1 logical)
    idx_d = ~isnan(T.cval_ddra) & ~isnan(T.sizeI_ddra);
    idx_g = ~isnan(T.cval_gray) & ~isnan(T.sizeI_gray);

    % Pull data vectors
    x_d = T.cval_ddra(idx_d);   y_d = T.sizeI_ddra(idx_d);
    x_g = T.cval_gray(idx_g);   y_g = T.sizeI_gray(idx_g);

    rd = T.ddra_ridge(idx_d);
    
    % Robust coercion to logical
    if iscell(rd),           rd = string(rd); end
    if iscategorical(rd),    rd = string(rd); end
    if isstring(rd) || ischar(rd)
        ridge_d = strcmpi(string(rd), "true") | strcmpi(string(rd), "1") | strcmpi(string(rd), "yes");
    elseif isnumeric(rd)
        ridge_d = rd ~= 0;
    elseif islogical(rd)
        ridge_d = rd;
    else
        ridge_d = false(size(x_d));
    end
    ridge_d = logical(ridge_d(:));   % ensure column vector


    figure; hold on; box on; grid on;

    % Style DDRA points; if ddra_ridge exists, split by ridge flag
    if ismember('ddra_ridge', T.Properties.VariableNames) && any(idx_d)
        % Non-ridge
        if any(~ridge_d)
            scatter(x_d(~ridge_d), y_d(~ridge_d), 60, 'o', 'LineWidth', 1.5, 'DisplayName','DDRA');
        end
        % Ridge
        if any(ridge_d)
            scatter(x_d(ridge_d),  y_d(ridge_d),  60, 'o', 'filled', 'LineWidth', 1.5, 'DisplayName','DDRA (ridge)');
        end
    else
        % No ridge info
        if any(idx_d)
            scatter(x_d, y_d, 60, 'o', 'LineWidth', 1.5, 'DisplayName','DDRA');
        end
    end

    % Gray
    if any(idx_g)
        scatter(x_g, y_g, 60, 's', 'LineWidth', 1.5, 'DisplayName','Gray');
    end

    xlabel('coverage on VAL (%)');
    ylabel('mean output-set size (interval width proxy)');
    title('Coverage vs Size (Pareto)');
    xlim([0 100]);
    legend('Location','SouthWest');

    % If you prefer y-log scale when sizes vary a lot:
    % set(gca,'YScale','log');

    % Annotate counts
    text(0.02, 0.98, sprintf('N_DDRA=%d, N_Gray=%d', numel(x_d), numel(x_g)), ...
        'Units','normalized','VerticalAlignment','top');
end
