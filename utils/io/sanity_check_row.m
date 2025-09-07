function sanity_check_row(results_dir, row_id)
    art = load(fullfile(results_dir,'artifacts', sprintf('row_%04d.mat', row_id)));

    % Recompute Gray reach widths for VAL and compare to per-step CSV
    sysG = art.sys_gray; sysT = art.sys_ddra; VAL = art.VAL;
    b  = 1;                                      % or whichever VAL block you want to inspect
    nk = size(art.VAL.y{b},1);
    
    % Gray initial set: zero-center VAL.R0 and shift by the block's x0
    R0_gray = zonotope([zeros(size(center(art.VAL.R0),1),1), generators(art.VAL.R0)]) + art.VAL.x0{b};
    paramsG = struct('R0', R0_gray, 'u', art.VAL.u{b}', 'tFinal', art.sys_gray.dt*(nk-1));
    
    % (optional, but common) match input dims
    if isprop(art.sys_gray,'nrOfInputs')
        du = art.sys_gray.nrOfInputs - size(paramsG.u,1);
        if du>0, paramsG.u = [paramsG.u; zeros(du,size(paramsG.u,2))]; end
        if du<0, paramsG.u = paramsG.u(1:art.sys_gray.nrOfInputs,:); end
    end
    
    % If Gray has disturbance channels, pass W in its disturbance space
    W_pred = []; 
    if isprop(art.sys_gray,'nrOfDisturbances') && art.sys_gray.nrOfDisturbances>0
        try
            if ~isempty(art.sys_gray.E) && isequal(size(art.sys_gray.E),size(art.sys_gray.B)) ...
                    && norm(art.sys_gray.E - art.sys_gray.B,'fro')<1e-12
                W_pred = art.VAL.U;  % E==B: pass U so E*W = B*U
            else
                W_pred = normalizeWForGray(art.sys_gray, art.W_eff);
            end
            paramsG.W = coerceWToSys(art.sys_gray, W_pred);
        catch
            % leave empty if anything goes wrong
        end
    end


    % Gray widths using your helper (from project): gray_infer_size_on_VAL
    [~, wid_gray_k] = gray_infer_size_on_VAL(sysG, ...
        testSuite_fromVAL(sysT, VAL), struct('shared',struct('n_k_val', numel(VAL.u{1}))), ...
        struct('R0', art.VAL.R0, 'W', art.W_eff));   % adapt to your signature

    PS = readtable(fullfile(results_dir,'summary_perstep.csv'));
    wps = PS.wid_gray(PS.row==row_id);
    fprintf('Row %d: max |recomputed - CSV| = %.3g\n', row_id, max(abs(wid_gray_k(:) - wps(:))));
end
