function MZ = make_W_matrix_zono(W, T)
% Build a matrix zonotope representing stacked disturbances across T samples.
% New CORA API: matZonotope(C, G3) with G3 of size (n x T x h).
% Falls back to old {cell of generators} signature if needed.

    % center: n x T zeros (disturbance mean is 0)
    n = size(center(W), 1);
    C = zeros(n, T);

    % fetch generator matrix of W robustly
    Gmat = [];
    try
        Gmat = generators(W);                 
    catch
        % fallbacks for older CORA objects
        if isprop(W, 'G') && ~isempty(W.G)
            Gmat = W.G;                        % n x h
        elseif isprop(W, 'Z') && ~isempty(W.Z)
            Z = W.Z;                           % [c | G1 ... Gh]
            if size(Z,2) >= 2
                Gmat = Z(:, 2:end);
            else
                Gmat = [];
            end
        else
            Gmat = [];
        end
    end

    h = size(Gmat, 2);                         % number of W generators

    % new-API path: pack as 3D array (n x T x h)
    G3 = zeros(n, T, h);
    for j = 1:h
        gj = Gmat(:, j);                      
        G3(:, :, j) = repmat(gj, 1, T);
    end

    % try new CORA constructor
    try
        MZ = matZonotope(C, G3);
        return
    catch
        % fallback for older CORA: cell of (n x T) generator matrices
        GW = cell(1, h*T);
        idx = 1;
        for j = 1:h
            gj = Gmat(:, j);
            for t = 1:T
                Gt = zeros(n, T);
                Gt(:, t) = gj;
                GW{idx} = Gt;
                idx = idx + 1;
            end
        end
        MZ = matZonotope(C, GW);
    end
end
