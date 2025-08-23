function MZ = build_Mw_matrix_zonotope(W, dim_x, T)
    CMw = zeros(dim_x, T);

    Z   = [W.center, W.generators];
    nG  = size(W.generators,2);
    GW  = cell(1, nG*T);

    idx = 1;
    for i=1:nG
        g = Z(:, i+1);
        G0 = [g, zeros(dim_x, T-1)]; GW{idx} = G0;
        for j=1:(T-1)
            GW{idx+j} = [GW{idx+j-1}(:,2:end), GW{idx+j-1}(:,1)];
        end
        idx = idx + T;
    end
    MZ = matZonotope(CMw, GW);
end
