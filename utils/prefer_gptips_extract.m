function cleanup = prefer_gptips_extract()
% Put GPTIPS' extract.m ahead of MATLAB's and clear cached resolution.

    % 1) Find desired GPTIPS extract
    gp_all = which('extract','-all');
    gp_idx = find(contains(gp_all, ...
        fullfile('global','thirdparty','gptips2','extract.m'), ...
        'IgnoreCase',true), 1);

    if isempty(gp_idx)
        error('prefer_gptips_extract:NotFound', ...
            'Could not find GPTIPS extract.m on the MATLAB path.');
    end
    gp_dir = fileparts(gp_all{gp_idx});

    % 2) Clear cached resolution & push GPTIPS path to the top
    clear extract          % clear function cache for 'extract'
    oldpath = path;
    cleanup = onCleanup(@() path(oldpath));  % auto-restore on function exit
    addpath(gp_dir, '-begin');
    rehash toolboxcache    % refresh path cache

    % 3) Sanity check
    all_now = which('extract','-all');
    assert(contains(all_now{1}, 'gptips2'), ...
        'Failed to prioritize GPTIPS extract.m; first is %s', all_now{1});
end
