function dm = signed_margin_point(y, S, Ddirs)
    % y: (p x 1), S set; approximate margin via support function
    % δ(y,S) ≈ max_u [u' y - h_S(u)] with sign
    su = support_of_zono(S, Ddirs);          % (Nd x 1)
    gap= (Ddirs.' * y) - su;                 % (Nd x 1)
    m  = max(gap);
    % if m<=0 point is inside; margin is negative distance-to-boundary proxy
    dm = m;   % you may record both dm and max(-gap) for inside depth
end
