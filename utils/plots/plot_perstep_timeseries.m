function plot_perstep_timeseries(cfg)
    [plots_dir, results_dir] = init_io(cfg);
    csv_perstep = fullfile(results_dir,'summary_perstep.csv');
    if ~exist(csv_perstep,'file'), warning('No per-step CSV found'); return; end
    T = readtable(csv_perstep);

    % pick the first finished row or use cfg.io.plot_rows(1) if present
    if isfield(cfg,'io') && isfield(cfg.io,'plot_rows') && ~isempty(cfg.io.plot_rows)
        rid = cfg.io.plot_rows(1);
    else
        rid = T.row(1);
    end
    Tk = T(T.row==rid, :);
    k  = Tk.k;

    % --- Widths vs k (DDRA, Gray, +True if artifact present)
    fW = figure('Name','Widths vs k','Color','w'); hold on; grid on;
    plot(k, Tk.wid_ddra, '-o', 'LineWidth', 1.6, 'DisplayName', 'DDRA (mean interval width)');
    plot(k, Tk.wid_gray, '-s', 'LineWidth', 1.6, 'DisplayName', 'RCSI/Gray (mean interval width)');

    % try to overlay true mean-width from artifact
    art = fullfile(results_dir, 'artifacts', sprintf('row_%04d.mat', rid));
    if exist(art, 'file')
        S = load(art, 'metrics');
        if isfield(S,'metrics') && isfield(S.metrics,'mw_true_k') && ~isempty(S.metrics.mw_true_k)
            plot(k, S.metrics.mw_true_k(:), '-^', 'LineWidth', 1.6, 'DisplayName', 'True (mean width)');
        end
    end
    xlabel('k (time step)'); ylabel('Mean output interval width');
    title(sprintf('Per-step widths — row %d', rid)); legend('Location','best');
    save_plot(fW, plots_dir, sprintf('perstep_widths_row_%04d', rid), ...
        'Formats', {'png','pdf'}, 'Resolution', 200);

    % --- Containment vs k (if enhanced metrics recorded)
    has_cov = any(~isnan(Tk.cov_ddra)) || any(~isnan(Tk.cov_gray));
    if has_cov
        fC = figure('Name','Containment vs k','Color','w'); hold on; grid on;
        if any(~isnan(Tk.cov_ddra)), plot(k, Tk.cov_ddra, '-o', 'LineWidth', 1.6, 'DisplayName', 'DDRA cover%'); end
        if any(~isnan(Tk.cov_gray)), plot(k, Tk.cov_gray, '-s', 'LineWidth', 1.6, 'DisplayName', 'RCSI/Gray cover%'); end
        xlabel('k (time step)'); ylabel('Containment (%)');
        title(sprintf('Per-step containment — row %d', rid)); legend('Location','best');
        save_plot(fC, plots_dir, sprintf('perstep_containment_row_%04d', rid), ...
            'Formats', {'png','pdf'}, 'Resolution', 200);
    end
end
