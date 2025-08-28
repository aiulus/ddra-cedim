function tf = hasGenerators(Z)
    tf = isa(Z,'zonotope') && ~isempty(generators(Z)) && size(generators(Z),2) > 0;
end
