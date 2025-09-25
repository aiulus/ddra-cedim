function s = support_of_zono(Z, Ddirs)
% Ddirs: (p x Nd); Z is (con/poly)zonotope – we coerce to zonotope
    Zz = toZono(Z);
    c  = center(Zz); G = generators(Zz);
    s  = support_zono_vec(Ddirs, c, G);    % (Nd x 1), you already have this
end
