function uk = get_uk(Ublk, k, m)
    if m==0, uk = zeros(0,1); return; end
    Ublk = squeeze(Ublk);
    if size(Ublk,1) == m
        uk = Ublk(:, min(k, size(Ublk,2)));
    elseif size(Ublk,2) == m
        uk = Ublk(min(k, size(Ublk,1)), :)';
    else
        error('get_uk: U has incompatible size for m=%d.', m);
    end
end