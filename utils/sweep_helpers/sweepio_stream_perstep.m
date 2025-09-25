function sweepio_stream_perstep(IO, row_index, wid_ddra_k, wid_gray_k, ratio_k, cov_ddra_k, cov_gray_k)
% Append per-step metrics with header on first use.

    i_ensure_perstep_header(IO.csv_perstep);

    nkv = max([numel(wid_gray_k), numel(wid_ddra_k), numel(ratio_k), ...
               numel(cov_ddra_k), numel(cov_gray_k), 1]);

    wid_ddra_k = i_padv(wid_ddra_k, nkv);
    wid_gray_k = i_padv(wid_gray_k, nkv);
    ratio_k    = i_padv(ratio_k,    nkv);
    cov_ddra_k = i_padv(cov_ddra_k, nkv);
    cov_gray_k = i_padv(cov_gray_k, nkv);

    fid = fopen(IO.csv_perstep,'a');
    for kk = 1:nkv
        fprintf(fid, '%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g\n', ...
            row_index, kk, wid_ddra_k(kk), wid_gray_k(kk), ratio_k(kk), ...
            cov_ddra_k(kk), cov_gray_k(kk));
    end
    fclose(fid);
end

function i_ensure_perstep_header(csv_perstep)
    if ~exist(csv_perstep,'file') || dir(csv_perstep).bytes==0
        fid = fopen(csv_perstep,'w');
        fprintf(fid,'row,k,wid_ddra,wid_gray,ratio_gray_true,cov_ddra,cov_gray\n');
        fclose(fid);
    end
end

function v = i_padv(v, n)
    if isempty(v), v = nan(n,1); return; end
    v = v(:);
    if numel(v) < n, v(end+1:n,1) = NaN; end
end
