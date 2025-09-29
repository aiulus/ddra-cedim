function Wpred = build_W_pred(sys_gray, ~, W_used)
% Same pass-through under E = I.
    if isprop(sys_gray,'nrOfDisturbances') && sys_gray.nrOfDisturbances > 0
        must( size(center(W_used),1) == size(sys_gray.A,1), ...
            'W dim must equal state dim when E=I');
        Wpred = W_used;
    else
        Wpred = [];
    end
end
