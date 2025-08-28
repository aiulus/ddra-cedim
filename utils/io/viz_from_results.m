function viz_from_results(data_dir, which_row, b, dims, tikz)
    % viz_from_results  Make CSV-based panels + one reach-set figure from saved artifacts.
    %
    % Usage:
    %   viz_from_results('experiments/results/data/kMSD_sample_size_sweep_grayLS_sweeps', 3, 1, [1 2], true)
    
    if nargin < 5, tikz = false; end
    
    csv_path = fullfile(data_dir, 'summary.csv');
    plots_dir = strrep(data_dir, fullfile('results','data'), fullfile('results','plots'));
    if ~exist(plots_dir,'dir'), mkdir(plots_dir); end
    
    % 1) CSV panels
    plot_sweeps_from_csv(csv_path, ...
        'SaveDir', fullfile(plots_dir, 'from_csv'), ...
        'TikZ', tikz, ...
        'Which', {'runtime','fc'});
    
    % 2) Reach-set figure from artifact (row index = which_row)
    artdir = fullfile(data_dir, 'artifacts');
    artfile = fullfile(artdir, sprintf('row_%04d.mat', which_row));
    if ~exist(artfile,'file')
        error('Artifact not found: %s. Enable cfg.io.save_artifacts = true and re-run the sweep.', artfile);
    end
    
    C = colorscheme('tum');
    plot_reachsets_2d_from_artifact(artfile, b, dims, ...
        'Colors', C, ...
        'SaveDir', fullfile(plots_dir, 'reach'), ...
        'TikZ', tikz, ...
        'Name', sprintf('row%04d_block%d', which_row, b));
end
