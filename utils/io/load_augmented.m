function [T, PS] = load_augmented(results_dir)
    T  = readtable(fullfile(results_dir,'summary.csv'));
    PS = readtable(fullfile(results_dir,'summary_perstep.csv'));

    % Drop rows that were intentionally skipped
    if ismember('skipped', T.Properties.VariableNames)
        T = T(~T.skipped, :);
    end

    % Ensure expected keys exist
    assert(all(ismember({'row'}, PS.Properties.VariableNames)), 'per-step CSV missing ''row''.');
    assert(all(ismember({'row'},  T.Properties.VariableNames)), 'summary CSV missing ''row''.');

    % Quick duplicate check
    assert(height(T) == numel(unique(T.row)), 'Duplicate rows in summary.csv');

    % Quick glance at the new E* fields
    disp('Row-level E* columns:');
    disp(T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'E')));
end
