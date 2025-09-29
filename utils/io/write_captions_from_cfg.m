function caps = write_captions_from_cfg(cfg)
%WRITE_CAPTIONS_FROM_CFG  Build captions from results/<save_tag>_sweeps/summary.csv
% Usage (end of sample_size.m / noise_model.m / scalability.m):
%   caps = write_captions_from_cfg(cfg);
%   % -> TeX file saved alongside figures; also returns strings in caps.*

    [plots_dir, results_dir] = init_io(cfg);
    sum_csv  = fullfile(results_dir, 'summary.csv');

    tex_out  = fullfile(plots_dir, sprintf('%s_captions.tex', cfg.io.save_tag));
    caps = make_sweep_captions(sum_csv, cfg, 'writeTex', true, 'texPath', tex_out);

    fprintf('[captions] Wrote %s\n', tex_out);
end
