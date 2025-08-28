function tf = contains_interval(y, Z, tol)
% True iff point y lies in the interval hull of set Z (zonotope/matZonotope/contSet).
% Conservative, fast, and identical across methods.
    if nargin < 3 || isempty(tol), tol = 1e-6; end
    I = interval(Z);                       % CORA interval() works on most contSets
    tf = all(y(:) >= I.inf(:) - tol & y(:) <= I.sup(:) + tol);
end