function TS = testSuite_fromDDRA(sys, R0, DATASET, n_k, n_m, n_s)
% testSuite_fromDDRA
% Maps DDRA DATASET blocks to a CORA-compatible test suite.
%
% Contract this function enforces for CORA/gray:
%   TS{m}.y : (n_k × n_y × n_s)  measured outputs for test case m
%   TS{m}.u : (n_k × n_u)        single *nominal* input for test case m
%   TS{m}.initialState : (n_x × 1) nominal initial state (center(R0) or 0)
%
% All *per-sample* inputs remain available in DATASET.U_blocks and are
% consumed when building VAL for DDRA; CORA/gray uses only the nominal u.

ny = DATASET.n_y; nu = DATASET.n_u; dt = sys.dt;
TS = cell(n_m,1);

use_center = (exist('center','file') && isa(R0,'zonotope'));
if use_center
    x0_nom = center(R0);
else
    x0_nom =  zeros(DATASET.dim_x,1);
end

for m = 1:n_m
    % y: gather all samples for this nominal test case
    y_m = zeros(n_k, ny, n_s);
    for s = 1:n_s
        idx = (m-1)*n_s + s;                    
        y_m(:,:,s) = DATASET.Y_blocks(:,:,idx).';   % transpose to (n_k×n_y)
    end

    % u: pick the *nominal* trajectory (2D)
    if isfield(DATASET,'U_nom_all') && numel(DATASET.U_nom_all) >= m && ~isempty(DATASET.U_nom_all{m})
        u_m = DATASET.U_nom_all{m}.';           % (n_k×n_u)
    else
        % fallback: use the first sample's input as nominal
        idx0 = (m-1)*n_s + 1;
        u_m  = DATASET.U_blocks(:,:,idx0).';     % (n_k×n_u)
    end

    TS{m} = testCase(y_m, u_m, x0_nom, dt, class(sys));
end
end
