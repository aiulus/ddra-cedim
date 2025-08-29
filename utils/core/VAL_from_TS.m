function VAL = VAL_from_TS(TS_val, DATASET_val)
    % Build synchronized validation record for EXACT (x0,u,y[,w]) points.

    B = numel(TS_val);
    VAL = struct('x0',{cell(1,B)}, 'u',{cell(1,B)}, 'y',{cell(1,B)}, 'w',{cell(1,B)});

    have_DS = (nargin >= 2) && ~isempty(DATASET_val);

    for b = 1:B
        % x0
        VAL.x0{b} = TS_val{b}.initialState;

        Ui = squeeze(TS_val{b}.u);        % (n_k × n_u)

        if size(Ui,2) > size(Ui,1)
            Ui = Ui.';                    % ensure (n_k × n_u)
        end
        VAL.u{b} = Ui;

        % y : prefer DATASET_val if present; otherwise use the TS field
        yi = [];
        if have_DS
            if isfield(DATASET_val,'Y_blocks')           % new layout: (ny × n_k × M)
                yi = DATASET_val.Y_blocks(:,:,b).';      % (n_k × ny)
            elseif isfield(DATASET_val,'y_blocks')       % rare legacy naming
                yi = DATASET_val.y_blocks(:,:,b).';
            end
        end
        if isempty(yi) && isfield(TS_val{b},'y') && ~isempty(TS_val{b}.y)
            ti = TS_val{b}.y;                            % (n_k × ny × 1 or s)
            yi = squeeze(ti);                            % (n_k × ny) if s==1
            if isvector(yi), yi = yi(:); end
        end
        VAL.y{b} = yi;  

        if have_DS && isfield(DATASET_val,'W_blocks')
            VAL.w{b} = DATASET_val.W_blocks(:,:,b);      % (nx × n_k)
        end
    end
end
