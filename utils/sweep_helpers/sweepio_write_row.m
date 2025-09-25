function IO = sweepio_write_row(IO, row)
% Write one summary row (append or in-memory). Initializes header on first call.

    if isempty(IO.hdr)
        IO.hdr = fieldnames(orderfields(row))';
        if IO.append
            fid = fopen(IO.csv_path,'w');  fprintf(fid, '%s\n', strjoin(IO.hdr, ','));  fclose(fid);
        else
            if isempty(IO.Ntot), error('sweepio_write_row: Ntot must be provided for in-memory mode.'); end
            IO.cells = cell(IO.Ntot, numel(IO.hdr));
        end
    end

    vals = i_row_to_cells(IO.hdr, row);

    if IO.append
        line = i_join_cell_as_csv(vals);
        fid  = fopen(IO.csv_path,'a');  fprintf(fid, '%s\n', line);  fclose(fid);
    else
        IO.cells(IO.rowi+1, :) = vals;
    end

    IO.rowi = IO.rowi + 1;
end

function vals = i_row_to_cells(hdr, row)
    vals = cell(1, numel(hdr));
    for j = 1:numel(hdr)
        key = hdr{j};
        if isfield(row, key)
            v = row.(key);
            if isstring(v), v = char(v); end
            vals{j} = v;
        else
            vals{j} = NaN;
        end
    end
end

function line = i_join_cell_as_csv(vals)
    parts = cell(1,numel(vals));
    for i = 1:numel(vals)
        v = vals{i};
        if isnumeric(v) && isscalar(v)
            parts{i} = sprintf('%.12g', v);
        elseif islogical(v) && isscalar(v)
            parts{i} = char(string(v));
        elseif ischar(v)
            parts{i} = i_csv_quote(v);
        elseif isstring(v)
            parts{i} = i_csv_quote(char(v));
        else
            try
                parts{i} = i_csv_quote(jsonencode(v));
            catch
                parts{i} = i_csv_quote('<unserializable>');
            end
        end
    end
    line = strjoin(parts, ',');
end

function s = i_csv_quote(c)
    if isempty(c), s = ''; return; end
    c = strrep(c, '"', '""');           % escape quotes
    s = ['"', c, '"'];
end
