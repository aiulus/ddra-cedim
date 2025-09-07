function sanity_check_row(results_dir, row_id)
    art = load(fullfile(results_dir,'artifacts', sprintf('row_%04d.mat', row_id)));

    % Recompute Gray reach widths for VAL and compare to per-step CSV
    sysG = art.sys_gray; sysT = art.sys_ddra; VAL = art.VAL;
    params = struct('R0', art.sys_gray.C\ (art.sys_gray.C*art.VAL.R0)); 

    % Gray widths using your helper (from project): gray_infer_size_on_VAL
    [~, wid_gray_k] = gray_infer_size_on_VAL(sysG, ...
        testSuite_fromVAL(sysT, VAL), struct('shared',struct('n_k_val', numel(VAL.u{1}))), ...
        struct('R0', art.VAL.R0, 'W', art.W_eff));   % adapt to your signature

    PS = readtable(fullfile(results_dir,'summary_perstep.csv'));
    wps = PS.wid_gray(PS.row==row_id);
    fprintf('Row %d: max |recomputed - CSV| = %.3g\n', row_id, max(abs(wid_gray_k(:) - wps(:))));
end
