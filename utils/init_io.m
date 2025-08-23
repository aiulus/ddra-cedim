function [plots_dir, results_dir] = init_io(cfg)
    if ~isfield(cfg,'io') || ~isfield(cfg.io,'save_tag')
        cfg.io.save_tag = datestr(now,'yyyymmdd_HHMMSS');
    end
    plots_dir   = fullfile('plots',   cfg.io.save_tag + "_sweeps");
    results_dir = fullfile('results', cfg.io.save_tag + "_sweeps");
    if ~exist(plots_dir, 'dir'),   mkdir(plots_dir);   end
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
end
