 function Wg = normalizeWForGray(sys, W_in)
    % Return a zonotope in disturbance space (nd = sys.nrOfDisturbances).
    % Backwards compatible: same signature, same call sites.
    % Improvements:
    %  - Never returns [] when nd>0 (avoids silent zeroing by caller).
    %  - Preserves magnitude when nd == 1; otherwise uses pinv(E).

    nd = sys.nrOfDisturbances;

    % If no disturbance channels, nothing to do.
    if nd == 0
        Wg = []; 
        return;
    end

    % If W_in is empty, return 0-disturbance of correct dim (explicit).
    if isempty(W_in)
        Wg = zonotope(zeros(nd,1));
        return;
    end

    % If already in disturbance space, pass through.
    d = size(center(W_in),1);
    if d == nd
        Wg = W_in; 
        return;
    end

    % Need the E matrix to map state-noise -> disturbance space
    if ~isprop(sys,'E') && ~isfield(sys,'E')
        warning('normalizeWForGray: missing E; using zero-disturbance.');
        Wg = zonotope(zeros(nd,1));
        return;
    end
    E = sys.E;              % size nx × nd
    nx = size(E,1);

    if d ~= nx || size(E,2) ~= nd
        warning('normalizeWForGray: dim mismatch (W_in:%d, E:%dx%d, nd:%d); using zero-disturbance.', d, size(E,1), size(E,2), nd);
        Wg = zonotope(zeros(nd,1));
        return;
    end

    % Map generators
    Gx = generators(W_in);          % nx × g (g can be 0)
    if isempty(Gx)
        Wg = zonotope(zeros(nd,1)); % explicit zero, correct dim
        return;
    end

    if nd == 1
        % Magnitude-preserving 1-D pullback:
        % choose scalars r_j so that ||E*r_j||_1 matches ||Gx(:,j)||_1.
        denom = max(sum(abs(E),1), eps);      % scalar since nd=1
        r = sum(abs(Gx),1) / denom;           % 1 × g
        Wg = zonotope(0, r);                  % 1-D zonotope with g gens
    else
        % General case: least-squares pullback
        R = pinv(E) * Gx;                     % nd × g

        % (Optional) gentle rescale to keep ||E*R|| close to ||Gx||:
        s_num = sum(abs(Gx),'all');
        s_den = max(sum(abs(E*R),'all'), eps);
        s = min( max(s_num/s_den, 0), 10 );   % cap to avoid blow-ups
        R = R * s;

        Wg = zonotope(zeros(nd,1), R);
    end
end
