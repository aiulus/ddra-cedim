function [Lvec, delta] = estimate_L_and_delta_from_data(D)
% Estimate per-dimension Lipschitz constants for the residual and a cover radius.
% D has fields Xminus (nx×N), Uminus (m'×N), Xplus (nx×N).

    nx = size(D.Xminus,1);
    N  = size(D.Xminus,2);
    if N < 2
        Lvec  = zeros(nx,1);
        delta = 0;
        return;
    end

    % local LS fit to get residuals r = Xplus - Mhat*[1;Xminus;Uminus]
    Phi  = [ones(1,N); D.Xminus; D.Uminus];            % (1+nx+m') × N
    Mhat = (D.Xplus) * pinv(Phi);                      % nx × (1+nx+m')
    R    = D.Xplus - Mhat*Phi;                         % residuals, nx × N
    Z    = [D.Xminus; D.Uminus];                       % regressor, (nx+m') × N

    % nearest-neighbor differences (O(N^2) simplified with one NN per point)
    Lvec  = zeros(nx,1);
    dmins = zeros(1,N);
    for j = 1:N
        dz  = Z - Z(:,j);
        d2  = sum(dz.^2,1); d2(j) = inf;
        [dmin, k] = min(d2);
        d  = sqrt(dmin);
        if d == 0, continue; end
        dmins(j) = d;

        dr = abs(R(:,j) - R(:,k));         % nx×1
        Lvec = max(Lvec, dr / d);          % elementwise max slope
    end

    % cover radius proxy: a robust typical neighbor distance
    good = dmins > 0 & isfinite(dmins);
    if any(good)
        delta = median(dmins(good));
    else
        delta = 0;
    end
end
