function plot_size_profiles_vs_k(kvec, prof_gray_mean, prof_ddra_mean, varargin)
    p = inputParser;
    addParameter(p,'title','Normalized interval size vs time');
    addParameter(p,'label_gray','RCSI');
    addParameter(p,'label_ddra','DDRA');
    parse(p,varargin{:});
    f = figure('Color','w'); hold on; grid on;
    plot(kvec, prof_gray_mean, '-s', 'Color',[0.85 0.33 0.10], 'LineWidth',1.6, 'DisplayName', p.Results.label_gray);
    plot(kvec, prof_ddra_mean, '-o', 'Color',[0.23 0.49 0.77], 'LineWidth',1.6, 'DisplayName', p.Results.label_ddra);
    xlabel('time step k'); ylabel('normalized set size');
    title(p.Results.title); legend('Location','best');
end
