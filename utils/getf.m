function v = getf(S, key, def)
% GETF  Fetch struct field (optionally nested) with default.
%   v = getf(S, 'name', def)
%   v = getf(S, {'a','b','c'}, def)   % nested path
    if nargin < 3, def = []; end
    v = def;
    if isempty(S), return; end
    if ischar(key) || isstring(key)
        if isstruct(S) && isfield(S, key), v = S.(key); end
        return;
    end
    % nested path
    cur = S;
    for i = 1:numel(key)
        k = key{i};
        if isstruct(cur) && isfield(cur, k)
            cur = cur.(k);
        else
            return;
        end
    end
    v = cur;
end
