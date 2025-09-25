function mw = mean_width_dir(Z, Ddirs)
% mean width w(S) ≈ mean_u [h_S(u)+h_S(-u)]
    su  = support_of_zono(Z, Ddirs);
    sneg= support_of_zono(Z, -Ddirs);
    mw  = mean(su + sneg, 'omitnan');
end
