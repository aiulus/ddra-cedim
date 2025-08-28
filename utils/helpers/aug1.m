function A = aug1(S)
    c = center(S); G = generators(S);  % S in R^{n+m}
    A = zonotope([1; c], blkdiag(zeros(1,size(G,2)), G));
end