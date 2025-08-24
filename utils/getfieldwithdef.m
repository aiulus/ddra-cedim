function v = getfieldwithdef(S, key, def)
% Return S.(key) if it exists; otherwise return def.
% - key: string/char for a single field, or cellstr for nested access {'a','b','c'}.
% - Works for structs; for non-structs returns def.

    % Nested path support
    if iscell(key)
        cur = S;
        for i = 1:numel(key)
            k = key{i};
            if isstruct(cur) && isfield(cur, k)
                cur = cur.(k);
            else
                v = def; return;
            end
        end
        v = cur; return;
    end

    % Single field
    if isstruct(S) && isfield(S, key)
        v = S.(key);
    else
        v = def;
    end
end
