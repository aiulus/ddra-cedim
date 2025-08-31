function ok = check_PE_order(TS, L)
    % Extract U as m x T across all training blocks
    if isfield(TS, 'U_blocks')
        U = cat(2, TS.U_blocks{:});
    elseif isfield(TS, 'U')
        U = TS.U;
    else
        error('TS missing U/U_blocks');
    end
    [m,T] = size(U);
    if T < 2*L
        ok = false; return;
    end
    % Block-Hankel of depth L
    H = [];
    for i = 1:L
        H = [H; U(:, i : (T-L+i))]; 
    end
    ok = (rank(H) == m*L);
end
