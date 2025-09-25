function idxGray = pick_gray_config(configs, C)
    want = "graySeq";
    try
        if isfield(C,'gray') && isfield(C.gray,'methodsGray') && ~isempty(C.gray.methodsGray)
            want = string(C.gray.methodsGray(1));
        end
    catch
    end
    idxGray = find(cellfun(@(c) isfield(c,'name') && want==string(c.name), configs), 1, 'first');
    if isempty(idxGray), idxGray = min(2, numel(configs)); end
end