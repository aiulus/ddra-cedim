function sweep_feed(md, v)
    if isnumeric(v) || islogical(v)
        % shape + class + raw bytes
        updateStr(md, class(v));
        updateIntVec(md, int64(size(v)));
        if ~isempty(v)
            b = typecast(v(:), 'uint8');  % column-major
            md.update(b);
        end
    elseif ischar(v)
        updateStr(md, 'char'); md.update(uint8(v(:)));
    elseif isstring(v)
        updateStr(md, 'string');
        for s = v(:).'
            sweep_feed(md, char(s));
        end
    elseif iscell(v)
        updateStr(md, 'cell'); updateIntVec(md, int64(size(v)));
        for i = 1:numel(v), sweep_feed(md, v{i}); end
    elseif isstruct(v)
        updateStr(md, 'struct'); updateIntVec(md, int64(size(v)));
        f = sort(fieldnames(v));
        for k = 1:numel(f)
            fn = f{k}; md.update(uint8(fn));
            for i = 1:numel(v)
                sweep_feed(md, v(i).(fn));
            end
        end
    elseif isa(v,'zonotope')
        % If hashing sets, reduce to center/generators
        try
            sweep_feed(md, center(v)); sweep_feed(md, generators(v));
        catch
            updateStr(md, class(v));
        end
    else
        % Fallback to class+num2str (last resort)
        updateStr(md, ['<', class(v), '>']);
        try
            sweep_feed(md, num2str(v));
        catch
        end
    end
end