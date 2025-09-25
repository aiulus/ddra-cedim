function SUMMARY = sweepio_finalize(IO)
% Close out and return summary table.
    if IO.append
        fprintf('Sweeps done. Rows: %d\nCSV -> %s\n', IO.rowi, IO.csv_path);
        SUMMARY = readtable(IO.csv_path);
    else
        writecell([IO.hdr; IO.cells(1:IO.rowi,:)], IO.csv_path);
        fprintf('Sweeps done. Rows: %d\nCSV -> %s\n', IO.rowi, IO.csv_path);
        SUMMARY = cell2table(IO.cells(1:IO.rowi,:), 'VariableNames', IO.hdr);
    end
end
