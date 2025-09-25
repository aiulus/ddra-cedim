function D = sample_unit_dirs(p, Nd, seed)
    if nargin<3 || isempty(seed), seed = 12345; end
    rng(seed,'twister');
    X = randn(p, Nd);
    n = sqrt(sum(X.^2,1)) + eps;
    D = X ./ n;
end