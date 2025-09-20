function use_noise = resolve_use_noise(S)
    if isfield(S,'use_noise') && ~isempty(S.use_noise)
        use_noise = logical(S.use_noise); return;
    end
    g = ~isfield(S,'noise_for_gray') || logical(S.noise_for_gray);
    d = ~isfield(S,'noise_for_ddra') || logical(S.noise_for_ddra);
    use_noise = g && d;
end
