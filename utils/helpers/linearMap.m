function Zout = linearMap(arg1, arg2)
% Accept linearMap(M,Z) or linearMap(Z,M); unwrap cells; delegate to CORA's mtimes.
    if iscell(arg1), arg1 = arg1{1}; end
    if iscell(arg2), arg2 = arg2{1}; end

    isSet = @(x) isa(x,'zonotope') || isa(x,'interval') || isa(x,'conZonotope') ...
               || isa(x,'polyZonotope') || isa(x,'zonoBundle') || isa(x,'ellipsoid');

    if isSet(arg1) && ~isSet(arg2), Z = arg1; M = arg2; else, M = arg1; Z = arg2; end
    Zout = M * Z;   % CORA overloading
end


