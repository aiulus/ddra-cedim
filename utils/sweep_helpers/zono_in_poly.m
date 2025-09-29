function ok = zono_in_poly(Z, H, h, tau)
% Exact Z ⊆ {y | H y <= h - tau*||H||_1}, Z is a (possibly translated) zonotope
    if nargin < 4, tau = 0; end
    c = center(Z); G = generators(zonotope(Z));   % accept pure linearMap results too
    if isempty(G), G = zeros(numel(c),0); end
    row1 = sum(abs(H),2);                   % ||H_i||_1
    rhs  = h(:) - tau*row1;
    lhs  = H*c + sum(abs(H*G),2);
    ok   = all(lhs <= rhs + 1e-9);          % small numerical slack
end
