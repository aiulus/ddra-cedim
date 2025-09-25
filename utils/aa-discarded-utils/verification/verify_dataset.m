function verify_dataset(DATASET, sys)
    nx = DATASET.dim_x; nu = DATASET.n_u; nk = DATASET.n_k; M = DATASET.n_blocks;
    A = sys.A; B = sys.B;
    bad = 0;
    
    for b = 1:M
        X = DATASET.X_blocks(:,:,b);        % nx × (nk+1)
        U = DATASET.U_blocks(:,:,b);        % nu × nk
    
        % 1) first two columns should NOT be identical after the fix
        if norm(X(:,2) - X(:,1)) < 1e-9
            fprintf('WARN: duplicate x0 at block %d (col1==col2)\n', b); bad = bad+1; 
        end
    
        % 2) dynamics consistency: A*X(:,t)+B*U(:,t) = X(:,t+1)  for t=1..nk
        for t = 1:nk
            x_next = A*X(:,t) + B*U(:,t);
            if norm(x_next - X(:,t+1), Inf) > 1e-6
                fprintf('FAIL dyn at block %d, t=%d: ||Ax+Bu - x_next||=%.3g\n', b, t, norm(x_next - X(:,t+1)));
                bad = bad + 1; break
            end
        end
    end
    
    if bad==0
        disp('OK: X_blocks/U_blocks dynamics and indexing consistent.');
    else
        fprintf('Found %d issues.\n', bad);
    end
end
