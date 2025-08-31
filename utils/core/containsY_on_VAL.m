function contain_pct = containsY_on_VAL(Ysets, VAL, tol)
    if nargin<3 || isempty(tol), tol = 1e-6; end
    num_in = 0; num_all = 0;
    B = numel(VAL.y);
    for b = 1:B
        Yb = VAL.y{b};  % (n_k×ny)
        for k = 1:min(size(Yb,1), numel(Ysets))
            yk = Yb(k,:).';
            if contains_interval(yk, Ysets{k}, tol), num_in = num_in + 1; end
            num_all = num_all + 1;
        end
    end
    contain_pct = (num_all>0)*100*num_in/max(1,num_all);
    if num_all==0, contain_pct = NaN; end
end
