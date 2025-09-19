% randPoint() - Wrapper for CORA versions that don't support it, currently
% not being used
function M = sample_U_points(U, n_k, extreme)
    if nargin < 3, extreme = false; end
    try
        if extreme, M = randPoint(U, n_k, 'extreme');
        else,       M = randPoint(U, n_k);
        end
    catch
        % fallback: draw column-by-column
        c = center(U); d = numel(c);
        M = zeros(d, n_k);
        for tt=1:n_k
            if extreme, M(:,tt) = randPoint(U,1,'extreme');
            else,       M(:,tt) = randPoint(U);
            end
        end
    end
end
