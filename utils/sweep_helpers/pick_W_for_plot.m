function Wplot = pick_W_for_plot(W_used, use_noise)
    if use_noise, Wplot = W_used;
    else, Wplot = zonotope(zeros(size(center(W_used),1),1));
    end
end