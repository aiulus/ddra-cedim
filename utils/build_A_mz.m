function A_mz = build_A_mz(A, mult)
% mult.type: 'diag' or 'full'
% mult.alpha_A: scalar generator scale
    if ~isfield(mult,'alpha_A'), mult.alpha_A = 0.0; end
    if mult.alpha_A == 0
        A_mz = matZonotope(A, {}); return;  % degenerate (no mult noise)
    end
    n = size(A,1); m = size(A,2);
    G = {};
    switch string(getfielddef(mult,'type','diag'))
        case "diag"
            % n generators on the diagonal
            for i=1:min(n,m)
                Gi = zeros(n,m); Gi(i,i) = mult.alpha_A;
                G{end+1} = Gi;
            end
        otherwise % "full"
            % one generator per entry (n*m) — can be large
            for i=1:n
                for j=1:m
                    Gi = zeros(n,m); Gi(i,j) = mult.alpha_A;
                    G{end+1} = Gi;
                end
            end
    end
    A_mz = matZonotope(A, G);
end
