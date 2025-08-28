function plot_artifacts_dir(artdir, rows, dims, savedir)
if nargin<2 || isempty(rows)
    S = dir(fullfile(artdir,'row_*.mat')); rows = 1:numel(S);
else
    S = dir(fullfile(artdir,'row_*.mat'));
end
if nargin<3 || isempty(dims), dims = [1 2]; end
if nargin<4 || isempty(savedir), savedir = fullfile(artdir,'..','..','plots','reach_offline'); end
if ~exist(savedir,'dir'), mkdir(savedir); end

for r = rows(:)'
    f = fullfile(artdir, sprintf('row_%04d.mat', r));
    if ~isfile(f), warning('Artifact missing: %s', f); continue; end
    art = load(f);
    % tolerate older runs without M_AB
    MAB = []; if isfield(art,'M_AB'), MAB = art.M_AB; end

    plot_reach_all_onepanel(art.sys_ddra, art.sys_gray, art.VAL, ...
        'Dims', dims, 'MAB', MAB, 'W', art.W_eff, ...
        'Save', fullfile(savedir, sprintf('row_%04d_y%d%d', r, dims(1), dims(2))));
end
end
