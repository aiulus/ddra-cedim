function lbl = rcsi_label_from_cfg(cfg)
% Returns a compact label like "graySeq", "grayLS", or "graySeq+grayLS".
    try
        m = string(cfg.gray.methodsGray);
        if numel(m) == 0
            lbl = 'graySeq';  % default fallback
        else
            lbl = strjoin(cellstr(m), '+');
        end
    catch
        lbl = 'graySeq';
    end
end
