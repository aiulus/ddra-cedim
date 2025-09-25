function tf = contains_points_approx(S, P)
% P: (p x N). Fast approximate membership for (con/poly)zonotope S.
% Try exact generator feasibility when possible; else fallback to directional test.
    Z = toZono(S); c = center(Z); G = generators(Z);
    V = P - c;                         % (p x N)
    if ~isempty(G)
        % Solve least squares for beta, check residual and infinity norm
        B = G \ V;                     % minimal-norm if square/overdetermined
        res = vecnorm(V - G*B, 2, 1);
        tf  = (res < 1e-8) & (max(abs(B),[],1) <= 1+1e-6);
    else
        tf  = all(abs(V) <= 1e-12, 1);
    end
    tf = tf(:).';
end
