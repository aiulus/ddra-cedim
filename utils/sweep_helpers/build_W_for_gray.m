function Wfg = build_W_for_gray(sys_cora, ~, W_used)
    % Map state-noise W_used -> disturbance space so E*Wfg \supseteq W_used.
    if isprop(sys_cora,'nrOfDisturbances') && sys_cora.nrOfDisturbances > 0
        Wfg = normalizeWForGray(sys_cora, W_used);
    else
        Wfg = [];
    end
end