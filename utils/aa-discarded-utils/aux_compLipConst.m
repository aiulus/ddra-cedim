function [gamma, L] = aux_compLipConst(fun, U, X0, steps, initpoints, dim_x, normType, opts)
% AUX_COMPLIPCONST  Estimate coord-wise Lipschitz constants for f and the sampled diameter.
% Usage:
%   [gamma,L] = aux_compLipConst(fun,U,X0,steps,initpoints,dim_x,normType)
%   [gamma,L] = aux_compLipConst(..., normType, opts) with opts.max_pairs (approx) etc.
%
% Inputs:
%   fun        : @(x,u) -> x_next   (noise-free dynamics)
%   U          : zonotope for inputs
%   X0         : zonotope for initial states
%   steps      : rollout length per trajectory
%   initpoints : number of independent trajectories
%   dim_x      : state dimension
%   normType   : 1, 2, or inf (default 2)
%   opts       : struct (optional) with fields:
%                   .max_pairs  = integer cap on random pairs per d (speeds up O(N^2))
%                   .rng_seed   = integer for reproducibility
%
% Outputs:
%   gamma(d)   : sampled diameter in z-space for coord d
%   L(d)       : Lipschitz estimate for coord d

    if nargin < 7 || isempty(normType), normType = 2; end
    if nargin < 8, opts = struct(); end

    % -------- sample inputs --------
    dim_u = numel(center(U));
    totalsamples = initpoints * steps;

    u = zeros(dim_u, totalsamples);
    for i = 1:totalsamples
        u(:,i) = randPoint(U);
    end

    % -------- roll out noise-free trajectories --------
    x_free = zeros(dim_x*initpoints, steps+1);
    idx = 1;
    for j = 1:dim_x:initpoints*dim_x
        x_free(j:j+dim_x-1,1) = randPoint(X0);
        for t = 1:steps
            x_free(j:j+dim_x-1,t+1) = fun(x_free(j:j+dim_x-1,t), u(:,idx));
            idx = idx + 1;
        end
    end

    % -------- vectorize into X0/X1 --------
    X0s = zeros(dim_x, totalsamples);  % x(k)
    X1s = zeros(dim_x, totalsamples);  % x(k+1)
    p0 = 1; p1 = 1;
    for j = 1:dim_x:initpoints*dim_x
        for t = 2:steps+1, X1s(:,p1) = x_free(j:j+dim_x-1,t);   p1 = p1 + 1; end
        for t = 1:steps,   X0s(:,p0) = x_free(j:j+dim_x-1,t);   p0 = p0 + 1; end
    end

    % -------- compute L and gamma per coordinate --------
    L = zeros(dim_x,1);
    gamma = zeros(dim_x,1);

    % Optional approximate mode
    use_approx = isfield(opts,'max_pairs') && ~isempty(opts.max_pairs) && opts.max_pairs > 0;
    if isfield(opts,'rng_seed') && ~isempty(opts.rng_seed)
        rng(opts.rng_seed,'twister');
    end

    if ~use_approx
        % Full O(N^2) scan
        for d = 1:dim_x
            Ld = 0; gd = 0;
            for i = 1:totalsamples
                z1 = [X0s(d,i); u(:,i)]; f1 = X1s(d,i);
                for j = 1:totalsamples
                    if i == j, continue; end
                    z2 = [X0s(d,j); u(:,j)]; f2 = X1s(d,j);
                    dz = norm(z1 - z2, normType);
                    if dz <= 0, continue; end
                    newL = norm(f1 - f2, normType) / dz;
                    if newL > Ld, Ld = newL; end
                    if dz > gd,   gd = dz;   end
                end
            end
            L(d) = Ld; gamma(d) = gd;
        end
        return;
    end

    % Approximate: sample random pairs (without replacement, across all pairs)
    N = totalsamples;
    % draw indices with replacement on [1..N] if pairs > N*(N-1)
    maxp = min(opts.max_pairs, N*(N-1));
    Ii = randi(N, maxp, 1);
    Jj = randi(N, maxp, 1);
    same = (Ii == Jj);
    if any(same)
        % ensure i != j
        Jj(same) = 1 + mod(Jj(same), N);
    end

    for d = 1:dim_x
        z1 = [X0s(d,Ii); u(:,Ii)];
        z2 = [X0s(d,Jj); u(:,Jj)];
        f1 =  X1s(d,Ii);
        f2 =  X1s(d,Jj);

        dz = vecnorm(z1 - z2, normType, 1);
        mask = (dz > 0);
        if any(mask)
            L(d) = max(abs(f1(mask) - f2(mask)) ./ dz(mask));
            gamma(d) = max(dz(mask));
        else
            L(d) = 0; gamma(d) = 0;
        end
    end
end
