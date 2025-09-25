function sweepio_save_artifact(IO, row_index, artifact)
% Persist one .mat file for plotting/debug.
    artdir = fullfile(IO.results_dir, 'artifacts');
    if ~exist(artdir,'dir'), mkdir(artdir); end
    save(fullfile(artdir, sprintf('row_%04d.mat', row_index)), ...
         '-struct','artifact','-v7.3');
end
