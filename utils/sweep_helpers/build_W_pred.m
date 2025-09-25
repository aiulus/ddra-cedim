function Wpred = build_W_pred(gray_sys, ~, W_used)
    if isprop(gray_sys,'nrOfDisturbances') && gray_sys.nrOfDisturbances > 0
        Wpred = normalizeWForGray(gray_sys, W_used);
    else
        Wpred = [];
    end
end