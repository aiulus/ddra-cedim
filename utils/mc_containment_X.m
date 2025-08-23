function pct = mc_containment_X(sys, R0, U, W, Xsets, C)
    % Monte Carlo state containment on validation-like draws
    n_m=C.shared.n_m; n_s=C.shared.n_s; n_k=C.shared.n_k; contain=0; total=0;
    for b=1:n_m*n_s
        x = randPoint(R0);
        for k=1:n_k
            u = randPoint(U); w = randPoint(W);
            x = sys.A*x + sys.B*u + w; total=total+1;
            if contains(Xsets{k}, x, 'approx', 1e-6), contain=contain+1; end
        end
    end
    pct = 100*contain/total;
end