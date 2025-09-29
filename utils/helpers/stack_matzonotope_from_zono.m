function MZ = stack_matzonotope_from_zono(Z, N)
% Build a matrix zonotope of size (nx × N) by sliding each generator col
% across time. Z is a zonotope(center,generators) in R^nx.
    c = center(Z); G = generators(Z); nx = size(c,1);
    if isempty(G), MZ = matZonotope(zeros(nx,N), {}); return; end
    Gcells = {};
    idx = 1;
    for ig = 1:size(G,2)
        g = G(:,ig);
        % seed with g at col 1
        G0 = [g, zeros(nx,N-1)];
        Gcells{idx} = G0; idx = idx+1;
        for j = 2:N
            Gcells{idx} = [Gcells{idx-1}(:,2:end), Gcells{idx-1}(:,1)];
            idx = idx+1;
        end
    end
    MZ = matZonotope(zeros(nx,N), Gcells);
end
