function rowi = emit_skip_row(rowi, reason, varargin)
    % Common skip pathway to avoid copy-paste
    % varargin lets allows passing name/value that end up on row fields
    % Relies on outer-scope vars already in run_sweeps: C, D, alpha_w, pe,
    % Zinfo, Tlearn, Tcheck, use_noise, csv_path, LM, rowi
    if exist('Tlearn','var'),        Tlearn_val = Tlearn;        else, Tlearn_val = NaN; end
    if exist('Tcheck','var'),        Tcheck_val = Tcheck;        else, Tcheck_val = NaN; end
    if exist('Tinfer','var'),        Tinfer_val = Tinfer;        else, Tinfer_val = NaN; end
    if exist('Tlearn_g','var'),      Tlearn_g_val = Tlearn_g;    else, Tlearn_g_val = NaN; end
    if exist('Tvalidate_g','var'),   Tvalidate_g_val = Tvalidate_g; else, Tvalidate_g_val = NaN; end
    if exist('Tinfer_g','var'),      Tinfer_g_val = Tinfer_g;    else, Tinfer_g_val = NaN; end

    rowi = rowi + 1;
    row = pack_row(C, D, alpha_w, pe, ...
        NaN, NaN, NaN, NaN, NaN, ...
        Zinfo.rankZ, Zinfo.condZ, ...
        Tlearn_val, Tcheck_val, Tinfer_val, Tlearn_g_val, Tvalidate_g_val, Tinfer_g_val);
    if ~isfield(row,'skipped'),     row.skipped = true; end
    if ~isfield(row,'skip_reason'), row.skip_reason = string(reason); end
    row.use_noise   = use_noise;

    % Ensure schema stability for columns that appear in non-skip rows
    row.cov_auc_gray = NaN;
    row.cov_auc_ddra = NaN;
    row.fv_gray_med  = NaN;
    row.fv_ddra_med  = NaN;
    
    row.dir_eps_med  = NaN;
    row.dir_eps_p90  = NaN;
    row.hout_med     = NaN;
    row.hout_p90     = NaN;
    
    % Ridge defaults (match non-skip rows)
    row.ddra_ridge        = false;
    row.ddra_lambda       = 0;
    row.ddra_kappa        = NaN;
    row.ddra_ridge_policy = "none";


    % attach any extra fields (e.g., ridge info)
    for ps = 1:2:numel(varargin)
        row.(varargin{ps}) = varargin{ps+1};
    end

    sweep_write_row(row, csv_path, LM);
end