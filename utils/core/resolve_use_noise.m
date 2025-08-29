function use_noise = resolve_use_noise(S)
    if isfield(S,'noise_for_gray') || isfield(S,'noise_for_ddra')
        g = ~isfield(S,'noise_for_gray') || logical(S.noise_for_gray);
        d = ~isfield(S,'noise_for_ddra') || logical(S.noise_for_ddra);
        use_noise = g && d;
    else
        use_noise = true;
    end
end
