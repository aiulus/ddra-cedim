function VAL = build_validation_bundle(TS_val, sys_cora)
% Collect the exact validation samples Gray will score and the signals needed
% to evaluate DDRA sets in OUTPUT space.
%
% VAL fields:
%   .y      : [n_k_val × dim_y × S]   stacked measured outputs across testcases
%   .u      : [n_k_val × dim_u × S]   stacked nominal inputs used in those testcases
%   .C, .D  : output matrices from sys_cora (linearSysDT)
%   .n_k    : number of time steps (n_k_val)
%
    S_total = 0;
    n_k = size(TS_val{1}.y,1);
    dim_y = size(TS_val{1}.y,2);
    dim_u = size(TS_val{1}.u,1);

    % count S across all testcases
    for m = 1:numel(TS_val)
        S_total = S_total + size(TS_val{m}.y, 3);
    end

    Y = zeros(n_k, dim_y, S_total);
    U = zeros(n_k, dim_u, S_total);

    idx = 0;
    for m = 1:numel(TS_val)
        Sm = size(TS_val{m}.y,3);
        Y(:, :, idx+(1:Sm)) = TS_val{m}.y;      % (k × y × s)
        % note: TS_val{m}.u is (dim_u × n_k) for one traj, replicate per sample s
        Um = permute(repmat(TS_val{m}.u, 1, 1, Sm), [2 1 3]); % (n_k × dim_u × Sm)
        U(:, :, idx+(1:Sm)) = Um;
        idx = idx + Sm;
    end

    VAL = struct();
    VAL.y  = Y;
    VAL.u  = U;
    VAL.n_k = n_k;

    % output map from linear system (RCSI reach works in output space already)
    if isa(sys_cora,'linearSysDT')
        VAL.C = sys_cora.C;
        VAL.D = sys_cora.D;
    else
        VAL.C = []; VAL.D = [];
    end
end
