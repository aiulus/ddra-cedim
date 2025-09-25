function write_row_init_if_needed(row, csv_path, LM)
    % Initializes header and in-memory cells (if needed) using caller scope.
    % Relies on outer-scope variables: hdr, cells, Ntot, rowi
    persistent initialized
    if isempty(initialized) || rowi == 1
        hdr = fieldnames(orderfields(row))';
        if LM.append_csv
            fid = fopen(csv_path, 'w'); fprintf(fid, '%s\n', strjoin(hdr, ',')); fclose(fid);
        else
            cells = cell(Ntot, numel(hdr));
        end
        initialized = true;
    end
end