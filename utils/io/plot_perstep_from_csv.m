function plot_perstep_from_csv(csv_perstep, rows)
    T = readtable(csv_perstep);
    rows = rows(:)';
    figure('Color','w'); tiledlayout(numel(rows),1,'TileSpacing','compact');
    for i = 1:numel(rows)
        R = T(T.row==rows(i),:);
        nexttile; hold on; grid on; box on;
        plot(R.k, R.wid_gray, '-', 'DisplayName','Gray width');
        plot(R.k, R.wid_ddra, '-', 'DisplayName','DDRA width');
        if any(isfinite(R.cov_gray)), yyaxis right; plot(R.k, R.cov_gray,'--','DisplayName','Gray cov'); end
        if any(isfinite(R.cov_ddra)), yyaxis right; plot(R.k, R.cov_ddra,':','DisplayName','DDRA cov'); end
        xlabel('step k'); ylabel('width'); title(sprintf('Row %d', rows(i)));
        legend('Location','best'); set(gca,'FontSize',11);
    end
end
