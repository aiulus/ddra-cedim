function [Lvec,gamma] = estimate_Lipschitz_from_data(fun, U, R0, steps, n_init, dim_x)
% Vector-wise L_i estimate using max secant slope across randomized sample pairs.
% fun(x,u) returns x_next; R0,U are zonotopes. Keep steps small (1-2), n_init modest.

tot = steps * n_init;
x0 = zeros(dim_x, n_init);
for j=1:n_init, x0(:,j) = randPoint(R0); end
u  = zeros(size(center(U),1), tot);
for i=1:tot, u(:,i) = randPoint(U); end

% Simulate free (noise-free) trajectories
X0 = zeros(dim_x, tot); X1 = zeros(dim_x, tot);
idx = 1;
for j=1:n_init
    x = x0(:,j);
    for s=1:steps
        z = x;                       % X_- per sample
        up = u(:,idx);
        x  = fun(x,up);              % X_+ per sample
        X0(:,idx) = z;
        X1(:,idx) = x;
        idx = idx + 1;
    end
end

% Pairwise slope (sample random subset for speed)
M = size(X0,2);
pick = min(5000, M*(M-1)/2);
Lvec  = zeros(dim_x,1);
gamma = zeros(dim_x,1);

for i=1:dim_x
    bestL = 0; bestG = 0;
    for t=1:pick
        a = randi(M); b = randi(M); if a==b, continue; end
        dz = [X0(i,a); u(:,a)] - [X0(i,b); u(:,b)];
        df = X1(i,a) - X1(i,b);
        nz = norm(dz,2);
        if nz > 0
            candL = abs(df) / nz;
            if candL > bestL, bestL = candL; end
            if nz > bestG, bestG = nz; end
        end
    end
    Lvec(i)  = bestL;
    gamma(i) = bestG;
end
end
