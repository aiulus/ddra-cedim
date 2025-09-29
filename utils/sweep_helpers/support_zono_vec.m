function s = support_zono_vec(Ddirs, c, G)
% Ddirs: either (ny × Nd) or (Nd × ny)
% c: (ny × 1), G: (ny × g)
    ny = numel(c);
    [r,cD] = size(Ddirs);
    if r == ny
        D = Ddirs;          % already (ny × Nd)
    elseif cD == ny
        D = Ddirs.';        % transpose to (ny × Nd)
    else
        error('support_zono_vec: Ddirs has size %dx%d; expected one dimension to equal ny=%d.', r, cD, ny);
    end
    if isempty(G), s = D.' * c; return; end
    s = D.' * c + sum(abs(D.' * G), 2);
end
