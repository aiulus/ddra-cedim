function audit_augmented(T, PS)
    fprintf('Rows in summary: %d\n', height(T));
    uk = unique(PS.row);
    fprintf('Rows in per-step: %d (unique row keys)\n', numel(uk));

    % Check per-step coverage consistency
    has_cov = any(contains(PS.Properties.VariableNames,'cov_'));
    fprintf('Per-step coverage columns present? %d\n', has_cov);

    % Ensure per-step k starts at 1
    if any(PS.k < 1)
        warning('Some per-step k < 1. Check VAL generation; k is expected to be 1..n_k_val.');
    end

    % Spot-check NaNs in E* columns
    Ecols = T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'E'));
    if ~isempty(Ecols)
        nn = varfun(@(x) sum(isnan(x)), T(:,Ecols));
        disp('NaNs per E* column:'); disp(nn);
    end
end
