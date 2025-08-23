function sizeI = gray_infer_size(sys_id, R0, U, C)
    sizeI = 0; nK = C.shared.n_k_val;
    for k=1:nK
        Rk = reach(sys_id, struct('R0',R0,'U',U,'tFinal',sys_id.dt*(k-1)), C.shared.options_reach);
        S = Rk.timePoint.set{k}; if ~isa(S,'contSet'), S = zonotope(S); end
        sizeI = sizeI + sum(abs(S.generators),'all');
    end
end