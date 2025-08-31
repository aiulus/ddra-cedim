function Y = map_to_output(Xset, C)
    try
        Y = linearMap(Xset, C);
    catch
        c = center(Xset); G = generators(Xset);
        Y = zonotope(C*c, C*G);
    end
end
