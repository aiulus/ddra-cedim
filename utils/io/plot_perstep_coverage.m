function fig = plot_perstep_coverage(PS)
    % Try to detect coverage columns
    vars = PS.Properties.VariableNames;
    cgray = vars(contains(vars,'cov_gray') | contains(vars,'coverage_gray'));
    cddra = vars(contains(vars,'cov_ddra') | contains(vars,'coverage_ddra'));

    assert(~isempty(cgray) && ~isempty(cddra), 'Coverage columns not found in per-step CSV.');

    % reshape into [nRows x nK]
    [klist,~,kid] = unique(PS.k,'stable');
    [rlist,~,rid] = unique(PS.row,'stable');
    nk = numel(klist); nr = numel(rlist);

    % helper to pivot
    pivot = @(colname) accumarray([rid,kid], PS.(colname), [nr,nk], @mean, NaN);

    Cg = pivot(cgray{1});     % Gray per-row, per-k
    Cd = pivot(cddra{1});     % DDRA per-row, per-k

    mu_g = nanmean(Cg,1);  se_g = nanstd(Cg,[],1)/sqrt(nr);  ci_g = 1.96*se_g;
    mu_d = nanmean(Cd,1);  se_d = nanstd(Cd,[],1)/sqrt(nr);  ci_d = 1.96*se_d;

    fig = figure('Color','w'); hold on; grid on;
    k = klist(:)';
    h1 = shadedErrorBar(k, mu_g, ci_g); h1.mainLine.DisplayName = 'Gray';
    h2 = shadedErrorBar(k, mu_d, ci_d); h2.mainLine.DisplayName = 'DDRA';
    xlabel('step k'); ylabel('coverage (%)'); title('Per-step coverage (VAL)'); legend('Location','southwest');
end
