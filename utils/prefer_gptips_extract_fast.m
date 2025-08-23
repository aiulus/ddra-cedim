function prefer_gptips_extract_fast()
% Ensure GPTIPS' extract.m is the one MATLAB resolves, but do it only once.

    persistent done okFirst
    if ~isempty(done) && done && okFirst
        return; % already prioritized and verified
    end

    % Find GPTIPS extract.m
    allx = which('extract','-all');
    gp_idx = find(contains(allx, fullfile('global','thirdparty','gptips2','extract.m'), ...
                           'IgnoreCase', true), 1);
    if isempty(gp_idx)
        error('prefer_gptips_extract_fast:NotFound', ...
              'Could not find GPTIPS extract.m on the MATLAB path.');
    end
    gp_dir = fileparts(allx{gp_idx});

    % If GPTIPS is already first, skip heavy work
    if contains(allx{1}, 'gptips2')
        okFirst = true; done = true; return;
    end

    % Otherwise put it in front and clear cache ONCE
    clear extract          % clear cached resolution for the symbol
    addpath(gp_dir, '-begin');
    rehash toolboxcache

    % sanity check
    allx2 = which('extract','-all');
    okFirst = contains(allx2{1}, 'gptips2');
    assert(okFirst, 'Failed to prioritize GPTIPS extract.m'); 
    done = true;
end
