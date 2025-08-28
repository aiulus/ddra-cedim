function tf = contains_interval(y, Z, tol)
% True iff point y lies in the interval hull of set Z.
% Works for CORA contSets (zonotope/matZonotope/etc.) and coerces others.

    if nargin < 3 || isempty(tol), tol = 1e-6; end
    y = y(:);

    % Coerce to a CORA set if needed
    if ~isa(Z,'contSet')
        Z = zonotope(Z);  % handles numeric [c G] forms and other inputs
    end

    % Empty set => never contained
    if representsa(Z,'emptySet')
        tf = false;
        return
    end

    Iv = interval(Z);                % CORA returns struct with .inf, .sup
    tf = all(y >= Iv.inf(:) - tol & y <= Iv.sup(:) + tol);
end
