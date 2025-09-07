function export_summary_tables(T, outdir)
    if ~exist(outdir,'dir'), mkdir(outdir); end

    % Example: keep a curated list of columns for the paper
    keep = T.Properties.VariableNames( ...
        contains(T.Properties.VariableNames, {'D','alpha_w','n_m','n_s','n_k','use_noise'}) | ...
        contains(T.Properties.VariableNames, {'cval_gray','cval_ddra','sizeI_gray','sizeI_ddra'}) | ...
        startsWith(T.Properties.VariableNames,'E') );

    Tp = T(:, keep);

    % Write CSV
    writetable(Tp, fullfile(outdir,'paper_summary.csv'));

    % Minimal LaTeX (quick & dirty)
    fid = fopen(fullfile(outdir,'paper_summary.tex'),'w');
    fprintf(fid,'\\begin{tabular}{%s}\\toprule\n', repmat('c',1,numel(Tp.Properties.VariableNames)));
    fprintf(fid,'%s\\\\\\midrule\n', strjoin(strrep(Tp.Properties.VariableNames,'_','\_'),' & '));
    for i=1:height(Tp)
        row = Tp(i,:);
        cells = strings(1,width(row));
        for j=1:width(row)
            v = row{1,j};
            if ismissing(v), cells(j) = '---';
            elseif isnumeric(v), cells(j) = sprintf('%.3g', v);
            else, cells(j) = string(v);
            end
        end
        fprintf(fid,'%s\\\\\n', strjoin(cells,' & '));
    end
    fprintf(fid,'\\bottomrule\\end{tabular}\n'); fclose(fid);
end
