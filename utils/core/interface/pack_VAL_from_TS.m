function VAL = pack_VAL_from_TS(TS)
%PACK_VAL_FROM_TS Canonical VAL: unroll (m x s) into a flat cell list.
% Stable order: for b=1..n_m, for s=1..n_s.
%
% Invariants:
% - VAL.x0, VAL.u, VAL.y are cell arrays of length n_m*n_s
% - Each y{i} is (n_k x q)
% - Each u{i} is (n_k x p)
% - Each x0{i} is (nx x 1)

    n_m = numel(TS);
    n_k = size(TS{1}.y,1);
    q   = size(TS{1}.y,2);
    n_s = size(TS{1}.y,3);
    p   = size(TS{1}.u,2);
    nx  = size(TS{1}.initialState,1);

    VAL = struct();
    VAL.x0 = cell(1, n_m*n_s);
    VAL.u  = cell(1, n_m*n_s);
    VAL.y  = cell(1, n_m*n_s);

    idx = 0;
    for b = 1:n_m
        % replicate x0 across s (x0 has no 3rd dim)
        x0b = TS{b}.initialState;
        ub  = TS{b}.u;
        for s = 1:n_s
            idx = idx + 1;
            VAL.x0{idx} = x0b;                      % (nx x 1)
            VAL.u{idx}  = ub;                       % (n_k x p)
            VAL.y{idx}  = TS{b}.y(:, :, s);         % (n_k x q)
        end
    end

    % (Optional) attach meta
    VAL.meta = struct('n_m', n_m, 'n_s', n_s, 'n_k', n_k, 'q', q, 'p', p, 'nx', nx);
end
