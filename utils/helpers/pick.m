function v = pick(x)
    % normalize scalar -> vector / passthrough cell
    if iscell(x), v = x; return; end
    if isscalar(x), v = x; else, v = x(:)'; end
end