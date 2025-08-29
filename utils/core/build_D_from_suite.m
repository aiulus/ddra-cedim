function D = build_D_from_suite(suite, p)
% Returns D.Xminus (nx×N), D.Uminus (m*(p+1)×N), D.Xplus (nx×N).

    if nargin < 2, p = 0; end
    Xm = []; Um = []; Xp = [];

    for b = 1:numel(suite)
        [Y, U] = case_to_YU(suite{b});     % (ny×T), (m×T)
        if isempty(Y) || isempty(U), continue; end
        [ny, T] = size(Y); m = size(U,1);
        if T < 2, continue; end

        for k = 1:(T-1)
            uk_stack = [];
            for lag = 0:p
                idx = max(1, k-lag);        % pad with first input
                uk_stack = [uk_stack; U(:,idx)];
            end
            Xm = [Xm, Y(:,k)];
            Um = [Um, uk_stack];
            Xp = [Xp, Y(:,k+1)];
        end
    end

    D = struct('Xminus', Xm, 'Uminus', Um, 'Xplus', Xp, ...
               'meta', struct('N', size(Xm,2), 'p', p));
end

