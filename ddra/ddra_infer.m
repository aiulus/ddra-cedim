function [Xsets, sizeI] = ddra_infer(sys, R0, U, W, M_AB, C)
    K = red_order(C);
    Xsets = cell(C.shared.n_k,1); Xk = R0; sizeI=0;
    for k=1:C.shared.n_k
        Xk = reduce(Xk,'girard',K);                             
        Xnext = M_AB * cartProd(Xk,U) + W;
        Xnext = reduce(Xnext,'girard',K);                       
        Xsets{k} = Xnext; 
        sizeI = sizeI + sum(abs(generators(Xnext)),'all');
        Xk = Xnext;
    end
end
