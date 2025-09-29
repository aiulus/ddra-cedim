function Wfg = build_W_for_gray(sys_gray, ~, W_used)
% With E = I, simply pass W through. If Gray has no disturbance channel, return [].
    if isprop(sys_gray,'nrOfDisturbances') && sys_gray.nrOfDisturbances > 0
        must( size(center(W_used),1) == size(sys_gray.A,1), ...
            'W dim must equal state dim when E=I');
        Wfg = W_used;
    else
        Wfg = [];
    end
end
