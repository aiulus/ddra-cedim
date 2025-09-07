function Wd = coerceWToSys(sys_, W_in)
    % Conservative preimage: build W_d in disturbance space s.t. E*W_d ≳ W_in
    nw = sys_.nrOfDisturbances;
    if nw == 0 || isempty(W_in), Wd = []; return; end
    if ~isa(W_in,'zonotope'), Wd = zonotope(zeros(nw,1)); return; end
    E = getfielddef(sys_,'E',[]);
    if isempty(E), Wd = zonotope(zeros(nw,1)); return; end
    c  = center(W_in); Gx = generators(W_in);
    rc = lsqminnorm(E, c);
    R  = zeros(nw, size(Gx,2)); res = 0;
    for j = 1:size(Gx,2)
        rj = lsqminnorm(E, Gx(:,j));
        R(:,j) = rj;
        res = max(res, norm(E*rj - Gx(:,j), inf));
    end
    eps_inf = res + 1e-12;                 % tiny padding to enforce inclusion
    if eps_inf > 0
        Wd = zonotope([rc, [R, eps_inf*eye(nw)]]);
    else
        Wd = zonotope([rc, R]);
    end
end
