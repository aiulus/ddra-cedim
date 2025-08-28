function use_noise = resolve_use_noise(shared)
% Prefer shared.use_noise if present; otherwise fall back to the old flags.
% If both old flags exist and disagree, throw to avoid asymmetrical comparisons.
    if isfield(shared,'use_noise')
        use_noise = logical(shared.use_noise);
        return;
    end
    has_g = isfield(shared,'noise_for_gray');
    has_d = isfield(shared,'noise_for_ddra');
    if ~has_g && ~has_d
        use_noise = true;   % default: noise ON 
        return;
    end
    if has_g && has_d && (logical(shared.noise_for_gray) ~= logical(shared.noise_for_ddra))
        error('Config conflict: noise_for_gray ~= noise_for_ddra. Set shared.use_noise instead.');
    end
    if has_g, use_noise = logical(shared.noise_for_gray); else, use_noise = logical(shared.noise_for_ddra); end
end
