function plot_perstep_coverage(PS)
% Plot per-step coverage for Gray and DDRA with 95% CI bands.

    assert(all(ismember({'row','k','cov_gray','cov_ddra'}, PS.Properties.VariableNames)), ...
        'PS must contain row,k,cov_gray,cov_ddra.');

    rows = unique(PS.row);
    nk   = max(PS.k);

    % Pivot into matrices (#rows x nk)
    Cg = nan(numel(rows), nk);
    Cd = nan(numel(rows), nk);
    for i = 1:numel(rows)
        idx = PS.row == rows(i);
        kk  = PS.k(idx);
        Cg(i, kk) = PS.cov_gray(idx);
        Cd(i, kk) = PS.cov_ddra(idx);
    end

    % Means (omit NaNs)
    mu_g = mean(Cg, 1, 'omitnan');
    mu_d = mean(Cd, 1, 'omitnan');

    % Per-step sample counts actually used
    n_g  = sum(~isnan(Cg), 1);
    n_d  = sum(~isnan(Cd), 1);

    % Std dev with 'omitnan'; use population normalization so 1 sample -> 0
    sd_g = std(Cg, 1, 1, 'omitnan');  % 2nd arg=1 => divide by N, not N-1
    sd_d = std(Cd, 1, 1, 'omitnan');

    % Standard error and 95% CI; guard n<1 and n<2
    se_g = sd_g ./ sqrt(max(n_g, 1));
    se_d = sd_d ./ sqrt(max(n_d, 1));
    ci_g = 1.96 * se_g;
    ci_d = 1.96 * se_d;

    kvec = 1:nk;

    figure; hold on; box on;
    % Shaded CI for Gray
    fill([kvec, fliplr(kvec)], [mu_g-ci_g, fliplr(mu_g+ci_g)], [0.9 0.9 1], ...
         'EdgeColor','none','FaceAlpha',0.4);
    % Shaded CI for DDRA
    fill([kvec, fliplr(kvec)], [mu_d-ci_d, fliplr(mu_d+ci_d)], [0.9 1 0.9], ...
         'EdgeColor','none','FaceAlpha',0.4);
    % Means
    plot(kvec, mu_g, '-', 'LineWidth', 2);
    plot(kvec, mu_d, '-', 'LineWidth', 2);

    ylim([0 100]);
    xlabel('time step k');
    ylabel('coverage (%)');
    legend({'Gray 95% CI','DDRA 95% CI','Gray mean','DDRA mean'}, 'Location','SouthWest');
    title('Per-step coverage with 95% CI (omit NaNs)');

    % Optional: annotate effective sample sizes
    txt = sprintf('Gray n_min=%d, n_max=%d | DDRA n_min=%d, n_max=%d', ...
                  min(n_g), max(n_g), min(n_d), max(n_d));
    annotation('textbox',[.15 .80 .7 .1],'String',txt,'LineStyle','none');
end
