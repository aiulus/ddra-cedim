function row = build_skip_row(C, D, alpha_w, pe, Zinfo, timers, use_noise, reason, extra_pairs)
% Construct a schema-stable skip row.
    if nargin < 9, extra_pairs = {}; end
    if ~isstruct(Zinfo), Zinfo = struct('rankZ',NaN,'condZ',NaN); end
    nf = @(s,f,def) (isfield(s,f) && ~isempty(s.(f))) * s.(f) + (~(isfield(s,f) && ~isempty(s.(f))))*def;
    Tlearn      = nf(timers,'Tlearn',NaN);
    Tcheck      = nf(timers,'Tcheck',NaN);
    Tinfer      = nf(timers,'Tinfer',NaN);
    Tlearn_g    = nf(timers,'Tlearn_g',NaN);
    Tvalidate_g = nf(timers,'Tvalidate_g',NaN);
    Tinfer_g    = nf(timers,'Tinfer_g',NaN);

    row = pack_row(C, D, alpha_w, pe, ...
            NaN, NaN, NaN, NaN, NaN, ...
            Zinfo.rankZ, Zinfo.condZ, ...
            Tlearn, Tcheck, Tinfer, Tlearn_g, Tvalidate_g, Tinfer_g);

    row.skipped            = true;
    row.skip_reason        = string(reason);
    row.use_noise          = logical(use_noise);
    row.cov_auc_gray       = NaN;  row.cov_auc_ddra = NaN;
    row.fv_gray_med        = NaN;  row.fv_ddra_med  = NaN;
    row.dir_eps_med        = NaN;  row.dir_eps_p90  = NaN;
    row.hout_med           = NaN;  row.hout_p90     = NaN;
    row.ddra_ridge         = false;
    row.ddra_lambda        = 0;    row.ddra_kappa   = NaN;
    row.ddra_ridge_policy  = "none";

    for k = 1:2:numel(extra_pairs)
        row.(extra_pairs{k}) = extra_pairs{k+1};
    end
end
