function sizeI = gray_infer_size_interval(cfg_gray, C)
    sys_id = cfg_gray.sys;
    R0_id  = cfg_gray.params.R0;
    U_id   = cfg_gray.params.U;

    sizeI = 0; nK = C.shared.n_k_val;
    params = struct('R0',R0_id,'U',U_id);
    opts   = C.shared.options_reach;

    params.tFinal = sys_id.dt * (nK-1);
    Rall = reach(sys_id, params, opts);
    for k=1:nK
        Zk = Rall.timePoint.set{k};
        if ~isa(Zk,'contSet'), Zk = zonotope(Zk); end
        Iv = interval(Zk);
        sizeI = sizeI + sum(abs(Iv.sup(:) - Iv.inf(:)));
    end
end
