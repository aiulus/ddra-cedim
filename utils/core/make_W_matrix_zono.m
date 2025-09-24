function MZ = make_W_matrix_zono(W, T)
    n = size(center(W), 1);
    C = repmat(center(W), 1, T);          % <-- tile actual center

    % robustly fetch W generators as an (n x h) matrix
    Gmat = [];
    try, Gmat = generators(W); catch, end
    if isempty(Gmat)
        if isprop(W, 'G') && ~isempty(W.G), Gmat = W.G;
        elseif isprop(W, 'Z') && ~isempty(W.Z), Gmat = W.Z(:,2:end);
        else, Gmat = zeros(n,0);
        end
    end
    h = size(Gmat,2);

    % time-sliced 3D generators (independent disturbances per step)
    G3 = zeros(n, T, h*T);
    for t = 1:T
        cols = (t-1)*h + (1:h);
        G3(:, t, cols) = Gmat;
    end

    MZ = matZonotope(C, G3);
end
