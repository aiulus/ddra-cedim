function TS = gray_testSuite_from_dataset(sys, R0, DATASET)
% Returns a cell array of "testCase-like" structs consumable by validateReach.
% Each TS{i} has:
%   .u             (n_k x n_u)          -- matches validateReach's expectation (transpose happens inside)
%   .y             (n_k x n_y x 1)      -- measurement samples for containment counting
%   .initialState  (dim_x x 1)          -- added to R0 inside validateReach

    assert(isstruct(DATASET) && isfield(DATASET,'n_k'), 'DATASET must be struct with n_k (got %s)', typestr(DATASET));


    n_k = DATASET.n_k;  M = DATASET.n_blocks;
    n_y = DATASET.n_y;

    TS = cell(M,1);
    for i = 1:M
        ui = permute(DATASET.U_blocks(:,:,i), [2 1]);      % (n_k x n_u)
        yi = permute(DATASET.Y_blocks(:,:,i), [2 1]);      % (n_k x n_y)
        xi0 = DATASET.x0_list(:,i);                        % (dim_x x 1)

        tc = struct();
        tc.u = ui;
        tc.y = reshape(yi, [n_k, n_y, 1]);                 % (n_k x n_y x 1)
        tc.initialState = xi0;                             % (dim_x x 1)
        tc.sampleTime   = sys.dt;  

        % minimal fields used by validateReach; R0 comes via configs{i}.params.R0
        TS{i} = tc;
    end
end
