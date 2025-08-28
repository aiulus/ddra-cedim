function MZ = build_Mw_matrix_zonotope(W, dim_x, T)
    % center
    C = zeros(dim_x, T);

    % 3-D generator array (n x m x h)
    nG = size(W.G,2);
    if nG == 0
        G3 = zeros(dim_x, T, 0);
    else
        G3 = zeros(dim_x, T, nG*T);
        h = 0;
        for i = 1:nG
            g = W.G(:,i);
            col = [g, zeros(dim_x, T-1)];
            for j = 0:(T-1)
                h = h+1;
                G3(:, :, h) = circshift(col, [0, j]);
            end
        end
    end

    MZ = matZonotope(C, G3);  
end
