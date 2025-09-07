function TS1 = pack_for_gray(TS)
% Pack {n_m} testCase objects into one (CORA gray expects this for linearSysDT).

    n_m = numel(TS);
    n_k = size(TS{1}.y,1);
    q   = size(TS{1}.y,2);
    n_s = size(TS{1}.y,3);
    p   = size(TS{1}.u,2);
    nx  = size(TS{1}.initialState,1);
    Ts  = TS{1}.sampleTime;

    % sanity checks
    for m = 1:n_m
        assert(size(TS{m}.y,1) == n_k && size(TS{m}.y,2) == q, 'pack_for_gray: y dims mismatch at m=%d', m);
        assert(size(TS{m}.u,1) == n_k && size(TS{m}.u,2) == p, 'pack_for_gray: u dims mismatch at m=%d', m);
        assert(size(TS{m}.initialState,1) == nx, 'pack_for_gray: x0 dim mismatch at m=%d', m);
        assert(abs(TS{m}.sampleTime - Ts) < 1e-12, 'pack_for_gray: sampleTime mismatch at m=%d', m);
    end

    y  = zeros(n_k, q, 0);
    u  = zeros(n_k, p, 0);
    x0 = zeros(nx, 1, 0);

    for m = 1:n_m
        y  = cat(3, y,  TS{m}.y);
        u  = cat(3, u,  repmat(TS{m}.u,  1, 1, n_s));
        x0 = cat(3, x0, repmat(TS{m}.initialState, 1, 1, n_s));
    end

    TS1 = { testCase(y, u, x0, Ts, 'linearSysDT') };
end
