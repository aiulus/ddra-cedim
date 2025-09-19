function MZ = aux_buildMatrixZonotope(W, dim_x, totalsamples)
    G = extractGenerators(W);           % dim_x x nG
    nG = size(G,2);
    C  = zeros(dim_x, totalsamples);

    if nG == 0
        try MZ = matZonotope(C, zeros(dim_x, totalsamples, 0)); return; end
        MZ = matZonotope(C, {}); return;  % fallback for older CORA
    end

    H  = nG * totalsamples;             % total generators
    G3 = zeros(dim_x, totalsamples, H); % new signature: (n x m x h)

    k = 1;
    for i = 1:nG
        vec  = G(:,i);
        base = [vec, zeros(dim_x, totalsamples-1)];  % vec in col 1
        for s = 0:totalsamples-1
            G3(:,:,k) = base(:, mod((0:totalsamples-1)-s, totalsamples)+1);
            k = k + 1;
        end
    end

    try
        MZ = matZonotope(C, G3);
    catch
        % Fallback to deprecated cell-of-generators
        GW = cell(1,H); k = 1;
        for i = 1:nG
            vec  = G(:,i);
            base = [vec, zeros(dim_x, totalsamples-1)];
            for s = 0:totalsamples-1
                GW{k} = base(:, mod((0:totalsamples-1)-s, totalsamples)+1);
                k = k + 1;
            end
        end
        MZ = matZonotope(C, GW);
    end
end

function G = extractGenerators(W)
    try
        G = generators(W);     % CORA ≥ 2023
    catch
        try G = W.generators;           % older field
        catch, ZW = W.Z; G = ZW(:,2:end);% very old
        end
    end
end
