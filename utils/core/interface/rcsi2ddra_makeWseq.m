function W_seq = rcsi2ddra_makeWseq(A, B, X0, Usets)
% Returns {W1,...,WT} with W1 = A*X0 \oplus B*U1,  Wk = B*Uk (k>=2).
% Usets can be:
%   - cell: Usets{k} is a zonotope in R^m
%   - numeric m×T (columns are time steps)
%   - numeric T×m (rows are time steps)

    m = size(B,2);

    % Auto-wrap numeric input sequences as point zonotopes
    if isnumeric(Usets)
        [r,c] = size(Usets);
        if r == m
            T = c;
            Ucell = cell(1,T);
            for k = 1:T, Ucell{k} = zonotope(Usets(:,k)); end
            Usets = Ucell;
        elseif c == m
            T = r;
            Ucell = cell(1,T);
            for k = 1:T, Ucell{k} = zonotope(Usets(k,:).'); end
            Usets = Ucell;
        else
            error('rcsi2ddra_makeWseq: numeric Usets has size %dx%d, which is incompatible with B (m=%d).', r, c, m);
        end
    end

    % Build W sequence
    T = numel(Usets);
    W_seq = cell(1,T);
    W_seq{1} = A*X0 + B*Usets{1};
    for k = 2:T
        W_seq{k} = B*Usets{k};
    end
end
