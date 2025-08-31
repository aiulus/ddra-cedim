function ok = check_PE_order(TS, L)
%CHECK_PE_ORDER  Verify mL-order PE on the training inputs (block-Hankel rank test).
% Accepts either TS.U_blocks{b} (m x n_k) or TS.U (m x T).
    if isfield(TS, 'U_blocks') && ~isempty(TS.U_blocks)
        U = cat(2, TS.U_blocks{:});    % m x T (concatenate all training blocks)
    elseif isfield(TS, 'U') && ~isempty(TS.U)
        U = TS.U;                      % m x T
    else
        error('TS missing U/U_blocks for PE check.');
    end
    [m,T] = size(U);
    if T < 2*L
        ok = false; return;
    end
    % Build block-Hankel of depth L
    H = [];
    for i = 1:L
        H = [H; U(:, i:(T-L+i))]; 
    end
    ok = (rank(H) == m*L);
end
