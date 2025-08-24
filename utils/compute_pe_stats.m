function PE = compute_pe_stats(Uminus, m, L)
    % Uminus: (m x N) stacked inputs across one big concatenation
    % m: #inputs; L: desired PE order
    N = size(Uminus,2);
    PE = struct('L',L,'N',N,'rankZ',NaN,'condZ',NaN, ...
                'hankel_rank',NaN,'hankel_full',false,'rank_frac',NaN);
    
    % Block-Hankel over the whole sequence
    if N >= L
        H = block_hankel(Uminus, L);
        r = rank(H);   PE.hankel_rank = r;
        PE.hankel_full = (r >= min(size(H))) && (r >= m*L);
        PE.rank_frac   = r / min(size(H));
    else
        PE.hankel_rank = 0; PE.rank_frac = 0; PE.hankel_full = false;
    end
    
    % Z = [X_minus; U_minus] rank already tracked elsewhere; we only add cond
    % (leave 'rankZ' to your existing Zinfo if you prefer)
    try
        PE.condZ = cond(double(Uminus*Uminus')); % proxy
    catch
        PE.condZ = NaN;
    end
end

function H = block_hankel(U, L)
    % U: (m x N) -> H: (mL x (N-L+1))
    [m, N] = size(U);
    cols = N - L + 1;
    H = zeros(m*L, max(cols,0));
    for i=1:L
        H((i-1)*m+1:i*m, :) = U(:, i:i+cols-1);
    end
end
