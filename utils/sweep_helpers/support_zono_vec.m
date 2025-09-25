function s = support_zono_vec(Ddirs, c, G)
% Ddirs: (p x Nd), c: (p x 1), G: (p x g)
    if isempty(G), s = Ddirs' * c; return; end
    s = Ddirs' * c + sum(abs(Ddirs' * G), 2);
end