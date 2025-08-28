function v = getfielddef(S, name, def)
    if nargin < 3, def = []; end
    if iscell(name)
        v = S;
        for i=1:numel(name)
            if isstruct(v) && isfield(v, name{i})
                v = v.(name{i});
            else
                v = def; return;
            end
        end
    else
        if isstruct(S) && isfield(S, name)
            v = S.(name);
        else
            v = def;
        end
    end
end
