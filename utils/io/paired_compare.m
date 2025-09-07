function paired_compare(T, metric_stub)
% PAIRED_COMPARE  Paired comparison scatter + diff stats for a metric.
% Robust to multiple column-name variants and missing columns.

    if nargin < 2
        error('paired_compare: need T and metric_stub (E1/E2/E3).');
    end
    metric_stubU = upper(string(metric_stub));
    valid = ["E1","E2","E3"];
    if ~ismember(metric_stubU, valid)
        error('paired_compare: metric_stub must be one of %s.', strjoin(valid,", "));
    end

    [g, d, label_text] = resolve_metric_columns(T, char(metric_stubU));

    assert(~isempty(g) && ~isempty(d), ...
        'Could not find matching Gray/DDRA columns for %s', metric_stub);

    % Vectorize and mask finite pairs
    g = double(g(:)); d = double(d(:));
    m = isfinite(g) & isfinite(d);
    g = g(m); d = d(m);
    if isempty(g)
        warning('No finite paired entries for %s; skipping.', metric_stubU);
        return;
    end

    % Scatter: DDRA (x) vs Gray (y)
    figure('Name',['Paired compare: ' char(metric_stubU)]);
    scatter(d, g, 60, 'filled'); hold on; grid on;
    lims = [min([g; d; 0]), max([g; d; 1])];
    plot(lims, lims, 'k--', 'LineWidth', 1);
    axis([lims lims]);
    xlabel(sprintf('DDRA — %s', label_text), 'Interpreter','tex');
    ylabel(sprintf('Gray — %s', label_text), 'Interpreter','tex');
    title(sprintf('Paired comparison: %s (n=%d)', metric_stubU, numel(g)));

    % Paired deltas
    delta = g - d;
    mu  = mean(delta,'omitnan');
    sd  = std(delta,'omitnan');
    med = median(delta,'omitnan');
    p90 = prctile(delta,90);

    fprintf('[%s] Gray - DDRA: mean=%.3g, std=%.3g, median=%.3g, p90=%.3g\n', ...
        metric_stubU, mu, sd, med, p90);

    % Histogram of deltas (guard small n)
    figure('Name',['Delta (Gray - DDRA): ' char(metric_stubU)]);
    if numel(delta) > 1
        histogram(delta, 'BinMethod','fd');
    else
        histogram(delta);
    end
    grid on;
    xlabel(sprintf('Gray - DDRA (%s)', label_text), 'Interpreter','tex');
    ylabel('count');
    title(sprintf('Paired deltas: %s', metric_stubU));
end

function [g, d, label_text] = resolve_metric_columns(T, stub)
    V = T.Properties.VariableNames;

    switch stub
        case 'E1'
            % Coverage AUC (%); accept both legacy and new names
            g = pickvar(T, {'cov_auc_gray','auc_cov_gray','cval_gray'});
            d = pickvar(T, {'cov_auc_ddra','auc_cov_ddra','cval_ddra'});
            label_text = 'coverage AUC (%)';

        case 'E2'
            % Directional support median (shape-aware conservatism)
            g = pickvar(T, {'supp_med_gray','dir_eps_med_gray','dir_eps_med'});
            d = pickvar(T, {'supp_med_ddra','dir_eps_med_ddra','dir_eps_med'});
            label_text = 'support(\hat{Y}) / support(Y_{true})';

        case 'E3'
            % Hausdorff (outer); prefer mean, fallback to p90
            g = pickvar(T, {'hdist_mean_gray','hdist_p90_gray'});
            d = pickvar(T, {'hdist_mean_ddra','hdist_p90_ddra'});
            label_text = 'outer Hausdorff (Ŷ vs Y_{true})';
    end

    if isempty(g) || isempty(d)
        fprintf('Available columns:\n  %s\n', strjoin(V, ', '));
    end
end

function x = pickvar(T, candidates)
    x = [];
    for i = 1:numel(candidates)
        if ismember(candidates{i}, T.Properties.VariableNames)
            x = T.(candidates{i});
            return;
        end
    end
end
