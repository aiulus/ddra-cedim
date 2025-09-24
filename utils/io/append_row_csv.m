function append_row_csv(csv_path, hdr, row)
    out = cell(1, numel(hdr));
    
    for j = 1:numel(hdr)
        if ~isfield(row, hdr{j}), v = []; else, v = row.(hdr{j}); end
        % v = row.(hdr{j});
        if isstring(v), v = char(v); end
        out{j} = local_to_str(v);
    end

    fid = fopen(csv_path,'a');
    assert(fid>0, 'append_row_csv: could not open %s for append', csv_path);
    fprintf(fid, '%s\n', strjoin(out, ','));
    fclose(fid);

    function s = local_to_str(v)
        if isempty(v)
            s = '';
        elseif isnumeric(v) && isscalar(v)
            s = num2str(v, '%.10g');
        elseif islogical(v) && isscalar(v)
            s = char(string(v));                
        elseif ischar(v)
            s = v;
        elseif isstring(v)
            s = char(v);
        elseif isnumeric(v) || islogical(v)
            % collapse small arrays into a single cell as space-separated
            s = strtrim(sprintf('%.10g ', v(:)));
        else
            % last-resort stringify
            s = char(string(v));
        end
        % sanitize commas/newlines so CSV stays 1-cell wide
        s = strrep(s, sprintf('\n'), ' ');
        s = strrep(s, ',', ';');
    end
end
