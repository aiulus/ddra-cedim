function dH = hausdorff_dir(SA, SB, Ddirs)
% symmetric Hausdorff via support gaps over sampled directions
    hA  = support_of_zono(SA, Ddirs);
    hB  = support_of_zono(SB, Ddirs);
    dH  = max(abs(hA - hB));   % L_infty over directions
end
