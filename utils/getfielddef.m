function val = getfielddef(S, name, default)
    if isstruct(S) && isfield(S,name)
        val = S.(name); 
    else
        val = default;
    end
end
