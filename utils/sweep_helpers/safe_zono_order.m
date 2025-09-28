function o = safe_zono_order(Z)
    try
        G = generators(Z); d = size(center(Z),1);
        o = max(1, size(G,2) / max(1,d));
    catch
        o = NaN;
    end
end
