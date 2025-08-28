function [nin, nout] = contain_points_in_sets(Y, Xsets)
    [~, nk, ns] = size(Y);
    nin = 0; nout = 0;
    for i = 1:ns
        for k = 2:nk
            S = Xsets{k};
            if ~isa(S,'contSet'); S = zonotope(S); end
            y = Y(:,k,i);
            if contains(S, y, 'approx', 1e-6)
                nin = nin + 1;
            else
                nout = nout + 1;
            end
        end
    end
end
