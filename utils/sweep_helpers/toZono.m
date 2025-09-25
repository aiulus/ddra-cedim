function Z = toZono(S)
% Convert various CORA set types to a (conservative) zonotope.
    if isa(S,'zonotope')
        Z = S;
    elseif isa(S,'polyZonotope') || isa(S,'conZonotope')
        Z = zonotope(S);   % CORA outer-approx
    else
        try
            Z = zonotope(S);
        catch
            error('toZono: unsupported set type %s', class(S));
        end
    end
end