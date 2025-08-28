function D = build_D_from_suite(S, p)
% Build data for DDRA from CORA test suites (NARX-friendly).
% Treat outputs y as "state". Include p input lags: Uminus = [u_k; u_{k-1}; ...; u_{k-p}].
%
% Returns:
%   D.Xminus : [ny x N]     with y_k
%   D.Uminus : [nu*(p+1) x N] with [u_k; u_{k-1}; ...; u_{k-p}]
%   D.Xplus  : [ny x N]     with y_{k+1}
%   D.meta   : struct with sizes and p

    if nargin < 2, p = 1; end

    Xm = []; Um = []; Xp = [];

    for m = 1:numel(S)
        Y = S{m}.y;                 % [nk x ny x ns]
        U_nom = S{m}.u;             % [nu x nk]  (nominal, same for all s)
        [nk, ny, ns] = size(Y);
        nu = size(U_nom,1);

        % build stacked inputs for all k (1..nk)
        Ustack = zeros(nu*(p+1), nk);
        for k = 1:nk
            ublock = zeros(nu*(p+1),1);
            for lag = 0:p
                kk = k - lag;
                if kk < 1, kk = 1; end  % pad with first input (mirrors common NARX handling)
                ublock(lag*nu + (1:nu)) = U_nom(:, kk);
            end
            Ustack(:,k) = ublock;
        end

        % collect samples for k = 1..nk-1  (y_{k+1} vs y_k, [u_k; ...; u_{k-p}])
        for s = 1:ns
            y_s = squeeze(Y(:,:,s))';  % [ny x nk]
            Xm  = [Xm,  y_s(:,1:nk-1)];
            Um  = [Um,  Ustack(:,1:nk-1)];
            Xp  = [Xp,  y_s(:,2:nk)];
        end
    end

    D = struct('Xminus', Xm, 'Uminus', Um, 'Xplus', Xp);
    D.meta = struct('ny', size(Xm,1), 'nu', nu, 'p', p);
end
