function s = size_interval_sum(Z)
% Standardized "conservatism" proxy = sum of interval widths of Z.
% Works for zonotope and matZonotope images (via interval()).
    I = interval(Z);
    s = sum(abs(I.sup(:) - I.inf(:)));
end