function Wd = coerceWToSys(sys_, W_in)
    % Conservative preimage: build W_d in disturbance space s.t. E*W_d \supseteq W_in
    nw = sys_.nrOfDisturbances;
    if nw == 0 || isempty(W_in), Wd = []; return; end
    if ~isa(W_in,'zonotope'), Wd = zonotope(zeros(nw,1)); return; end

    d_in = size(center(W_in),1);
    if d_in == nw
        % Already in disturbance space
        Wd = W_in; 
        return;
    end

    E = getfielddef(sys_,'E',[]);
    if isempty(E)
        Wd = zonotope(zeros(nw,1)); 
        return;
    end
    if d_in ~= size(E,1)
        warning('coerceWToSys: dim(W_in)=%d incompatible with rows(E)=%d. Using zero set.', d_in, size(E,1));
        Wd = zonotope(zeros(nw,1)); 
        return;
    end

    % LSQ + residual padding to enforce inclusion
    c  = center(W_in); Gx = generators(W_in);
    rc = lsqminnorm(E, c);
    R  = zeros(nw, size(Gx,2)); res = 0;
    for j = 1:size(Gx,2)
        rj = lsqminnorm(E, Gx(:,j));
        R(:,j) = rj;
        res = max(res, norm(E*rj - Gx(:,j), inf));
    end
    eps_inf = res + 1e-12;
    if eps_inf > 0
        Wd = zonotope([rc, [R, eps_inf*eye(nw)]]);
    else
        Wd = zonotope([rc, R]);
    end
end
