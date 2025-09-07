function TS1 = pack_for_gray(TS)
% Pack {n_m} testCase objects (each: y n_k×q×n_s, u n_k×p, x0 n×1)
% into a single testCase (y n_k×q×(n_m n_s), u n_k×p×(n_m n_s), x0 n×1×(n_m n_s))
% which is what CORA gray expects for linearSysDT.

    n_m = numel(TS);
    n_k = size(TS{1}.y,1);
    q   = size(TS{1}.y,2);
    n_s = size(TS{1}.y,3);
    p   = size(TS{1}.u,2);
    nx  = size(TS{1}.initialState,1);

    y  = zeros(n_k, q, 0);
    u  = zeros(n_k, p, 0);
    x0 = zeros(nx, 1, 0);

    for m = 1:n_m
        % stack measurements across samples
        y  = cat(3, y,  TS{m}.y);                            % (n_k×q×n_s)
        % repeat the same nominal input and x0 for each s
        u  = cat(3, u,  repmat(TS{m}.u,  1, 1, n_s));         % (n_k×p×n_s)
        x0 = cat(3, x0, repmat(TS{m}.initialState, 1, 1, n_s)); % (nx×1×n_s)
    end

    % name is just a label; keep CORA's convention
    TS1 = { testCase(y, u, x0, TS{1}.sampleTime, 'linearSysDT') };
end
