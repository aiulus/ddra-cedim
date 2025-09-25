function stream_perstep(csv_perstep, row_index, wid_ddra_k, wid_gray_k, ratio_k, cov_ddra_k, cov_gray_k)
    ensure_perstep_header(csv_perstep);
    nkv = max([numel(wid_gray_k), numel(wid_ddra_k), numel(ratio_k), numel(cov_ddra_k), numel(cov_gray_k)]);
    if isempty(wid_ddra_k), wid_ddra_k = nan(nkv,1); end
    if isempty(wid_gray_k), wid_gray_k = nan(nkv,1); end
    if isempty(ratio_k),    ratio_k    = nan(nkv,1); end
    if isempty(cov_ddra_k), cov_ddra_k = nan(nkv,1); end
    if isempty(cov_gray_k), cov_gray_k = nan(nkv,1); end

    fid_ps = fopen(csv_perstep,'a');
    for kk = 1:nkv
        fprintf(fid_ps, '%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g\n', ...
            row_index, kk, wid_ddra_k(min(kk,end)), wid_gray_k(min(kk,end)), ...
            ratio_k(min(kk,end)), cov_ddra_k(min(kk,end)), cov_gray_k(min(kk,end)));
    end
    fclose(fid_ps);
end