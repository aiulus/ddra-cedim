function save_overview_text(outdir, tag, row)
    fid=fopen(fullfile(outdir, tag + "_summary.txt"),'w');
    fprintf(fid,'dyn=%s D=%d nm=%d ns=%d nk=%d pe=%s', row.dyn, row.D, ...
        row.n_m, row.n_s, row.n_k, row.pe_mode);
    fprintf(fid,'DDRA contain_val=%.2f%% sizeI=%.4g | rank=%g cond=%g', ...
        row.cval_ddra, row.sizeI_ddra, row.rankZ_ddra, row.condZ_ddra);
    fprintf(fid,'GRAY  ctrain=%.2f%% cval=%.2f%% sizeI=%.4g', ...
        row.ctrain_gray, row.cval_gray, row.sizeI_gray);
    fclose(fid);
end