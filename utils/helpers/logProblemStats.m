function stats = logProblemStats(P, tag)
% P: struct for linprog/quadprog (A, b, Aeq, beq, lb, ub), or prob2struct output
% tag: char label for the experiment
    if isfield(P,'A');    A = P.A; else; A = []; end
    if isfield(P,'Aeq'); Aeq = P.Aeq; else; Aeq = []; end
    nvar = numel(P.f);
    mI  = size(A,1);  mE = size(Aeq,1);
    M   = [A; Aeq];   % full constraint matrix
    stats = struct('tag',tag,'nvar',nvar,'mI',mI,'mE',mE,...
        'nnz',nnz(M),'density',nnz(M)/numel(M),'time',nan);
    fprintf('[%s] nvar=%d, mI=%d, mE=%d, nnz=%d (dens=%.4f)\n', ...
        tag, nvar, mI, mE, stats.nnz, stats.density);
end
