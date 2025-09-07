function Mw = rcsi2ddra_buildMw(A, B, X0, Usets)
% Build matrix zonotope Mw (size n x T) from RCSI sets:
%   W0 = A*X0 \oplus B*U0,  Wk = B*Uk (k>=1),
% and place each Wk as column k of Mw.
%
% Inputs:
%   A,B   : system matrices (n×n and n×m)
%   X0    : zonotope (n-dim)
%   Usets : 1×T cell, each Usets{k} is a zonotope in R^m
%
% Output:
%   Mw    : matZonotope of size n×T

    assert(iscell(Usets) && ~isempty(Usets), 'Usets must be a nonempty cell array');
    T = numel(Usets);

    % per-step noise column-sets
    W = cell(1,T);
    W{1} = A*X0 + B*Usets{1};     % W0
    for k = 2:T
        W{k} = { B*Usets{k} };    % Wk
    end

    % Assemble matZonotope Mw whose column k equals zonotope W{k}
    n  = length(center(W{1}));
    Cw = zeros(n, T);
    Gs = {};                      % collect generator "slices" (each n×T)

    for k = 1:T
        ck = center(W{k});
        Gk = generators(W{k});            % n×gk (gk generators, possibly 0)
        Cw(:,k) = ck;

        gk = size(Gk,2);
        for i = 1:gk
            Gi = zeros(n, T);
            Gi(:,k) = Gk(:,i);           % put generator i only in column k
            Gs{end+1} = Gi;               
        end
    end

    if isempty(Gs)
        % Degenerate (deterministic) case: no generators
        Gw = zeros(n, T, 0);
    else
        Gw = cat(3, Gs{:});              % n×T×Ng
    end

    Mw = matZonotope(Cw, Gw);
end
