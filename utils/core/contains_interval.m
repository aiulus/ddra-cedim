function tf = contains_interval(y, Yset, tol)
    if nargin<3 || isempty(tol), tol = 1e-6; end
    I = interval(Yset);
    try, lo = infimum(I); hi = supremum(I); catch, lo = I.inf; hi = I.sup; end
    tf = all(y <= hi + tol) && all(y >= lo - tol);
end
