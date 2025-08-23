function [Xsets, sizeI] = ddra_infer(sys, R0, U, W, M_AB, C)
    Xsets = cell(C.shared.n_k,1); Xk = R0; sizeI=0;
    for k=1:C.shared.n_k
        Xk = reduce(Xk,'girard',400);
        Xk = M_AB * cartProd(Xk,U) + W;
        Xsets{k} = Xk; sizeI = sizeI + sum(abs(generators(Xk)),'all');
    end
end