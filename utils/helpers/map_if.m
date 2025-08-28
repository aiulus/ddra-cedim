function map_if(options, apx, fname)
    if isfield(apx, fname)
        options.approx.(fname) = apx.(fname);
    end
end