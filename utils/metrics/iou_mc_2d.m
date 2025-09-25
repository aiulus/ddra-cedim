function [iou, hitsA, hitsB, hitsI] = iou_mc_2d(SA, SB, Ndims, Nmc, seed)
    if nargin<4, Nmc=3e4; end
    if nargin<5, seed=1234; end
    rng(seed,'twister');

    % bounding box from interval hull of the union
    IA = interval(SA); IB = interval(SB);
    try, lo = min(infimum(IA), infimum(IB)); hi = max(supremum(IA), supremum(IB));
    catch, lo = min(IA.inf, IB.inf); hi = max(IA.sup, IB.sup);
    end
    lo = lo(Ndims); hi = hi(Ndims);

    X = lo(1) + (hi(1)-lo(1))*rand(Nmc,1);
    Y = lo(2) + (hi(2)-lo(2))*rand(Nmc,1);
    P = [X Y].';

    inA = contains_points_approx(SA, P);
    inB = contains_points_approx(SB, P);

    hitsA = mean(inA); hitsB = mean(inB);
    hitsI = mean(inA & inB);
    denom = hitsA + hitsB - hitsI + eps;
    iou   = hitsI ./ denom;
end
