% RCSI Report Builder
% --------------------
% One-shot script that assembles a PDF report from the sweep outputs
% produced by your RCSI/DDRA/Gray experiments.
%
% What it does:
%   1) Loads results_dir/summary.csv and summary_perstep.csv.
%   2) Creates a few key figures (coverage, width trends, ratios, case overlays).
%   3) Emits a standalone LaTeX report (report.tex) that embeds tables & figures.
%   4) Optionally compiles to PDF if pdflatex is on PATH.
%
% Usage:
%   - Set the CONFIG section below and run this file in MATLAB.
%   - The script is defensive: if a resource is missing, it skips that section.
%
% Notes:
%   - Requires only MATLAB + your CORA path (for optional overlay plots that
%     re-use reach sets from artifacts). If artifacts are present (saved via
%     sweepio_save_artifact), the script will render overlays; otherwise they
%     are skipped gracefully.
%   - You can re-run this anytime; figures and .tex are overwritten.

%% ============ CONFIG ============
results_dir   = fullfile(pwd, 'ddra-cedim/experiments/results');     % where summary.csv lives
plots_dir     = fullfile(results_dir, 'plots');
artifacts_dir = fullfile(results_dir, 'artifacts');
report_name   = 'rcsi_report';                % basename (no extension)

% Figures to produce (toggle individually)
MAKE.figure_coverage   = true;   % scatter / bars of containment & ratios
MAKE.figure_perstep    = true;   % per-step width trends (uses summary_perstep.csv)
MAKE.figure_overlays   = true;   % overlays from artifacts (if present)
MAKE.figure_system_map = true;   % count of runs per system (sanity viz)

% How many overlay panels to include (if artifacts exist)
MAX_OVERLAY = 3;                 % pick up to N artifacts for case figures

% Which columns to show in the main summary table
TABLE_COLS = {'dyn','type','n_k','n_m','n_s', ...
              'cval_gray','cval_ddra','ctrain_gray', ...
              'ratio_gray_true','sizeI_gray','sizeI_ddra', ...
              'Tvalidate_g','T_ddra'};  %#ok<*NASGU>

%% ============ PREP ============
if ~exist(results_dir,'dir')
    error('results_dir not found: %s', results_dir);
end
if ~exist(plots_dir,'dir'), mkdir(plots_dir); end

sum_csv  = fullfile(results_dir, 'summary.csv');
per_csv  = fullfile(results_dir, 'summary_perstep.csv');
tex_path = fullfile(results_dir, [report_name '.tex']);

T  = table(); TP = table();
if exist(sum_csv,'file')
    T = readtable(sum_csv);
else
    warning('Missing %s — the report will have stubs only.', sum_csv);
end
if exist(per_csv,'file')
    TP = readtable(per_csv);
else
    if MAKE.figure_perstep
        warning('Missing %s — per-step figure will be skipped.', per_csv);
        MAKE.figure_perstep = false;
    end
end

% Column resolver (case-insensitive, many aliases)
pick = @(tbl, opts) i_pickcol(tbl, opts);

% Canonical metrics (resolve best-effort from whatever headers are present)
col.dyn      = pick(T, {'dyn','system','sys','dynamics'});
col.type     = pick(T, {'type','uncertainty','set_type'});
col.n_k      = pick(T, {'n_k','nk','steps'});
col.n_m      = pick(T, {'n_m','nm','blocks','cases'});
col.n_s      = pick(T, {'n_s','ns','samples'});
col.cval_g   = pick(T, {'cval_gray','contain_val_gray','cov_gray_val','containment_gray_val','containment_gray'});
col.ctrain_g = pick(T, {'ctrain_gray','contain_train_gray','cov_gray_train'});
col.cval_d   = pick(T, {'cval_ddra','contain_val_ddra','cov_ddra_val','containment_ddra_val','containment_ddra'});
col.size_g   = pick(T, {'sizeI_gray','width_gray','size_gray','mean_width_gray'});
col.size_d   = pick(T, {'sizeI_ddra','width_ddra','size_ddra','mean_width_ddra'});
col.ratio    = pick(T, {'ratio_gray_true','ratio_gray_ddra','ratio'});
col.Tg       = pick(T, {'Tvalidate_g','time_validate_gray','t_val_gray'});
col.Td       = pick(T, {'T_ddra','time_ddra','t_ddra'});

%% ============ FIGURES ============
fig_paths = struct();

% --- Figure 1: Coverage vs Ratio (scatter) ---
if MAKE.figure_coverage && ~isempty(T)
    try
        x = i_col_or_nan(T, col.ratio);
        y = i_col_or_nan(T, col.cval_g);
        z = i_col_or_nan(T, col.cval_d);
        dynlab = i_str_or_blank(T, col.dyn);

        f1 = figure('Color','w','Name','Coverage_vs_Ratio'); hold on; grid on; box on;
        s1 = scatter(x, y, 28, 'filled', 'MarkerFaceAlpha', 0.7); %#ok<NASGU>
        s2 = scatter(x, z, 28, 'filled', 'MarkerFaceAlpha', 0.7);
        xlabel('Width ratio (gray vs true or ddra)');
        ylabel('Containment on VAL [%]');
        legend({'Gray','DDRA'}, 'Location','best');
        title('Containment vs. width ratio');

        % annotate a few extremes by dyn name
        [~,idx] = maxk(y, min(6, numel(y)));
        for k = idx(:)'
            text(x(k), y(k), ['  ' dynlab{k}], 'FontSize', 8, 'Interpreter','none');
        end
        fig_paths.coverage = fullfile(plots_dir, 'fig_coverage.png');
        exportgraphics(f1, fig_paths.coverage, 'Resolution', 200);
        close(f1);
    catch ME
        fprintf('Coverage figure failed: %s', ME.message);
    end
end

% --- Figure 2: Per-step width trends ---
if MAKE.figure_perstep && ~isempty(TP)
    try
        % Choose first row index that has non-NaN series
        row_ids = unique(TP.row);
        pick_rows = row_ids(1:min(4, numel(row_ids)));

        f2 = figure('Color','w','Name','PerStep_Widths');
        tiledlayout(numel(pick_rows), 1, 'TileSpacing','compact','Padding','compact');
        for i = 1:numel(pick_rows)
            rr = pick_rows(i);
            S = TP(TP.row==rr, :);
            nexttile; hold on; grid on; box on;
            plot(S.k, S.wid_gray, '-', 'LineWidth', 1.2);
            plot(S.k, S.wid_ddra, '-', 'LineWidth', 1.2);
            ylabel(sprintf('row %d', rr));
            if i==1, title('Per-step width (gray vs ddra)'); end
            if i==numel(pick_rows), xlabel('time step k'); end
            legend({'gray','ddra'}, 'Location','best');
        end
        fig_paths.perstep = fullfile(plots_dir, 'fig_perstep.png');
        exportgraphics(f2, fig_paths.perstep, 'Resolution', 200);
        close(f2);
    catch ME
        fprintf('Per-step figure failed: %s', ME.message);
    end
end

% --- Figure 3: Overlays from artifacts (best-effort) ---
if MAKE.figure_overlays && exist(artifacts_dir,'dir')
    mats = dir(fullfile(artifacts_dir, 'row_*.mat'));
    mats = mats(~[mats.isdir]);
    if ~isempty(mats)
        take = min(MAX_OVERLAY, numel(mats));
        for i = 1:take
            matpath = fullfile(mats(i).folder, mats(i).name);
            try
                % Use your helper to draw a compact 3-panel overlay
                overlay_reach_sets(matpath, 'Dims', [1 2], 'Reduce', red_order(i_guess_cfg()), ...
                    'SaveBase', fullfile(plots_dir, ['overlay_' i_strip_ext(mats(i).name)]), 'Show', false);
            catch ME
                warning('Overlay failed for %s: %s', mats(i).name, ME.message);
            end
        end
    else
        warning('No artifacts found in %s; skipping overlays.', artifacts_dir);
    end
end

% --- Figure 4: System frequency map ---
if MAKE.figure_system_map && ~isempty(T) && ~isempty(col.dyn)
    try
        dyns = i_str_or_blank(T, col.dyn);
        u = unique(dyns);
        count = zeros(numel(u),1);
        for i=1:numel(u)
            count(i) = sum(strcmp(dyns, u{i}));
        end
        [count,ord] = sort(count,'descend'); u = u(ord);
        f4 = figure('Color','w','Name','Systems');
        bar(count);
        set(gca,'XTick',1:numel(u), 'XTickLabel', u, 'XTickLabelRotation', 30);
        ylabel('# runs'); title('Runs per system'); grid on; box on;
        fig_paths.systems = fullfile(plots_dir, 'fig_systems.png');
        exportgraphics(f4, fig_paths.systems, 'Resolution', 200);
        close(f4);
    catch ME
        fprintf('System map figure failed: %s', ME.message);
    end
end

%% ============ TABLE SNIPPETS (LATEX) ============
% Build a compact top-N table (by highest VAL containment for gray)
latex_table_path = fullfile(plots_dir, 'table_top.tex');
if ~isempty(T)
    try
        n = height(T);
        score = i_col_or_nan(T, col.cval_g);
        [~,ord] = sort(score, 'descend', 'MissingPlacement','last');
        keep = ord(1:min(12,n));
        Sub = T(keep, :);

        % Columns to export (only those that exist)
        cols_exist = i_cols_exist(Sub, TABLE_COLS);
        Sub2 = Sub(:, cols_exist);
        i_write_tabular_tex(Sub2, latex_table_path, 'Top runs by Gray containment (VAL)');
    catch ME
        fprintf('LaTeX table build failed: %s', ME.message);
    end
end

%% ============ WRITE LATEX REPORT ============
tex = i_latex_preamble(report_name);

% Title block
tex = [tex, sprintf('\\title{%s: RCSI Reachset Conformance Report}\\author{Auto-generated}\\date{\\today}\\begin{document}\\maketitle\n', i_tex_esc(report_name))];

% Abstract
tex = [tex, ['\\begin{abstract}\n' ...
             'This report summarizes reachset-conformance experiments comparing a gray-box RCSI approach against a direct discrete reachable-approximation (DDRA) baseline. ' ...
             'We report identification/validation containment, output-set width proxies, and per-step trends across systems and uncertainty regimes. ' ...
             'All results are generated from CSV logs and artifacts produced by the sweep harness.\\end{abstract}\n\n']];

% Overview
tex = [tex, '\\section{Setup \\& Metrics}\n'];
tex = [tex, '\\begin{itemize}\n'];
tex = [tex, '\\item \\textbf{Containment (\\%)} on train/val via CORA\\texttt{\\char`\_}validateReach with robustness margin inherited from the experiment config.\n'];
tex = [tex, '\\item \\textbf{Size proxy} is the per-step output-interval width aggregated (sum/mean) across outputs, then averaged across time and blocks.\n'];
tex = [tex, '\\item \\textbf{Ratio} compares Gray vs. True (or DDRA when True absent), lower is tighter.\n'];
tex = [tex, '\\end{itemize}\n'];

% Quick counts
if ~isempty(T)
    tex = [tex, sprintf('There are %d runs in \\texttt{%s}.\\\n', height(T), i_tex_esc(results_dir))];
end

% Insert main table (if created)
if exist(latex_table_path,'file')
    tex = [tex, '\\section{Top Results (by VAL containment)}\n'];
    tex = [tex, sprintf('\\input{%s}\n', i_tex_relpath(tex_path, latex_table_path))];
end

% Insert figures (only if present)
tex = i_add_fig(tex, fig_paths, 'coverage',  'Containment vs. width ratio');
tex = i_add_fig(tex, fig_paths, 'perstep',   'Per-step width trends');
tex = i_add_fig(tex, fig_paths, 'systems',   'Runs per system');

% Overlay figures (if any were saved by overlay_reach_sets)
over_pngs = dir(fullfile(plots_dir, 'overlay_*_overlays.png'));
if ~isempty(over_pngs)
    tex = [tex, '\\section{Case Overlays (True vs DDRA vs Gray)}\n'];
    for i = 1:numel(over_pngs)
        p = fullfile(over_pngs(i).folder, over_pngs(i).name);
        tex = [tex, i_fig_block(i_tex_relpath(tex_path, p), over_pngs(i).name, 0.9)]; %#ok<AGROW>
    end
end

% Appendix: raw schema
if ~isempty(T)
    tex = [tex, '\\appendix\\section{CSV Schema}\n'];
    tex = [tex, '\\begin{verbatim}\n'];
    tex = [tex, strjoin(T.Properties.VariableNames, ', '), '\n'];
    tex = [tex, '\\end{verbatim}\n'];
end

tex = [tex, '\\end{document}\n'];

% Write .tex
fid = fopen(tex_path,'w'); fprintf(fid, '%s', tex); fclose(fid);
fprintf('LaTeX written to: %s\n', tex_path);

% Try to compile
[ok, log] = i_try_pdflatex(tex_path);
if ok
    fprintf('PDF generated: %s\n', strrep(tex_path, '.tex', '.pdf'));
else
    fprintf('pdflatex not found or failed. You can compile manually.\n');
    if ~isempty(log), fprintf('pdflatex output:\n%s\n', log); end
end

%% ================= helpers =================
function name = i_strip_ext(name)
    [~,name] = fileparts(name);
end

function rp = i_tex_relpath(texfile, other)
% Return a relative path from texfile's folder to 'other'
    [ptex,~] = fileparts(texfile);
    rp = char(string(relpath(other, ptex)));
end

function esc = i_tex_esc(s)
    esc = string(s);
    esc = strrep(esc, '_', '\\_');
    esc = char(esc);
end

function tex = i_latex_preamble(jobname)
    tex = ['% Auto-generated by rcsi\_report\_builder.m\n' ...
           '\\documentclass[11pt,a4paper]{article}\n' ...
           '\\usepackage[margin=1in]{geometry}\n' ...
           '\\usepackage{graphicx}\n' ...
           '\\usepackage{booktabs}\n' ...
           '\\usepackage{siunitx}\n' ...
           '\\usepackage{hyperref}\n' ...
           '\\hypersetup{colorlinks=true, allcolors=blue}\n' ...
           '\\setlength{\\parskip}{4pt}\n' ...
           '\\setlength{\\parindent}{0pt}\n'];
end

function tex = i_add_fig(tex, fig_paths, key, caption)
    if isfield(fig_paths, key) && exist(fig_paths.(key),'file')
        tex = [tex, sprintf('\\section{%s}\n', caption)];
        tex = [tex, i_fig_block(fig_paths.(key), caption, 0.9)];
    end
end

function blk = i_fig_block(path_rel_or_abs, caption, scale)
    if nargin<3, scale=0.95; end
    blk = sprintf(['\\begin{figure}[t]\n' ...
                   '  \\centering\\includegraphics[width=%0.2f\\linewidth]{%s}\n' ...
                   '  \\caption{%s}\n' ...
                   '\\end{figure}\n'], scale, i_tex_esc(path_rel_or_abs), i_tex_esc(caption));
end

function c = i_pickcol(T, names)
% Find the first variable in names (case-insensitive). Return '' if none.
    c = '';
    if isempty(T), return; end
    v = string(T.Properties.VariableNames);
    for i = 1:numel(names)
        j = find(strcmpi(v, names{i}), 1, 'first');
        if ~isempty(j)
            c = T.Properties.VariableNames{j}; return;
        end
    end
end

function v = i_col_or_nan(T, cname)
    if isempty(T) || isempty(cname) || ~ismember(cname, T.Properties.VariableNames)
        v = nan(height(T),1); return;
    end
    v = T.(cname);
    if iscell(v)
        try v = cellfun(@double, v); catch, v = nan(height(T),1); end
    end
end

function v = i_str_or_blank(T, cname)
    if isempty(T) || isempty(cname) || ~ismember(cname, T.Properties.VariableNames)
        v = repmat({''}, height(T),1); return;
    end
    raw = T.(cname);
    if isstring(raw), v = cellstr(raw); elseif iscellstr(raw), v = raw; else, v = cellstr(string(raw)); end
end

function cols_exist = i_cols_exist(T, cols)
    cols_exist = cols(ismember(cols, T.Properties.VariableNames));
end

function i_write_tabular_tex(Tsub, outpath, caption)
    if nargin<3, caption=''; end
    f = fopen(outpath,'w');
    vars = Tsub.Properties.VariableNames;
    n = numel(vars);
    fprintf(f, '%% Auto-generated table\n');
    fprintf(f, '\\begin{table}[t]\n');
    fprintf(f, '\\centering\n');
    fprintf(f, '\\small\n');
    fprintf(f, '\\begin{tabular}{%s}\n', repmat('l',1,n));
    fprintf(f, '\\toprule\n');
    % header
    for j=1:n
        lab = regexprep(vars{j}, '_', '\\_');
        if j<n, fprintf(f, '%s & ', lab); else, fprintf(f, '%s \\ \\ \n', lab); end
    end
    fprintf(f, '\\midrule\n');
    % rows (stringify)
    for i=1:height(Tsub)
        for j=1:n
            val = Tsub{i,j};
            if iscell(val), val = val{1}; end
            if isstring(val), val = char(val); end
            if isnumeric(val)
                if isfinite(val)
                    s = sprintf('%g', val);
                else
                    s = '';
                end
            else
                s = char(string(val));
            end
            s = regexprep(s, '_', '\\_');
            if j<n, fprintf(f, '%s & ', s); else, fprintf(f, '%s \\ \\ \n', s); end
        end
    end
    fprintf(f, '\\bottomrule\n');
    if ~isempty(caption)
        fprintf(f, '\\caption{%s}\n', regexprep(caption,'_','\\_'));
    end
    fprintf(f, '\\end{tabular}\n');
    fprintf(f, '\\end{table}\n');
    fclose(f);
    fprintf('Table written: %s\n', outpath);
end

function [ok, out] = i_try_pdflatex(texfile)
    [p,~,b] = fileparts(texfile);
    cwd = pwd; cd(p);
    try
        [s,o] = system(sprintf('pdflatex -interaction=nonstopmode -halt-on-error %s.tex', b));
        ok = (s==0);
        out = o;
    catch ME
        ok = false; out = ME.message;
    end
    cd(cwd);
end

function cfg = i_guess_cfg()
% Best-effort to retrieve a config for red_order(). If not available,
% return a struct with nested defaults.
    cfg = struct();
    cfg.shared = struct('options_reach', struct('zonotopeOrder', 60));
    cfg.lowmem = struct('zonotopeOrder_cap', 60);
end
