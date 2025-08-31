function TS1 = split_tc_samples(tc)
% Return cell array with S single-sample testCases from multi-sample tc
    S = size(tc.u,3);
    TS1 = cell(1,S);
    for s = 1:S
        u_s = tc.u(:,:,s);                      % (n_k × n_u)
        y_s = tc.y(:,:,s);                      % (n_k × n_y)
        TS1{s} = testCase(y_s, u_s, tc.initialState, tc.sampleTime, class(tc));
    end
end
