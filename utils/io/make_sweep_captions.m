function caps = make_sweep_captions(summary_csv, cfg, varargin)
%MAKE_SWEEP_CAPTIONS  Create LaTeX-ready figure captions from a sweep summary CSV.
% Returns a struct of strings and (optionally) writes a TeX file with macros.
%
% Required:
%   summary_csv : path to experiments/results/data/<tag>_sweeps/summary.csv
%   cfg         : the cfg used for the run (for budgets / noise / reduction text)
%
% Name-Value:
%   'writeTex'  : true/false (default true)
%   'texPath'   : path to write (default: same folder as CSV, <save_tag>_captions.tex)
%
% Output fields (strings):
%   caps.fidcons      : caption for Fidelity/Conservatism vs <axis>
%   caps.runtime      : caption for Runtime panels vs <axis>
%   caps.shapeaware   : caption for Shape-aware metrics vs <axis> (if available)
%   caps.perstep      : caption for Per-step panel (if available)

    p = inputParser;
    addParameter(p,'writeTex', true);
    addParameter(p,'texPath', "");
    parse(p, varargin{:});
    writeTex = p.Results.writeTex;

    T = readtable(summary_csv);

    % --- detect sweep axis (the one that varies) ---
    cand = {'n_m','n_s','n_k','alpha_w','D'};
    axisName = '';
    for i = 1:numel(cand)
        if any(strcmpi(T.Properties.VariableNames, cand{i}))
            v = unique(double(T.(cand{i})));
            if numel(v) > 1, axisName = cand{i}; break; end
        end
    end
    if isempty(axisName)
        % fallbacks: take first numeric col with >1 unique
        ncols = T.Properties.VariableNames;
        for i = 1:numel(ncols)
            col = ncols{i};
            if isnumeric(T.(col)) && numel(unique(T.(col)))>1
                axisName = col; break;
            end
        end
    end
    axisVals = [];
    if ~isempty(axisName), axisVals = unique(roundn(double(T.(axisName)), -12)); end

    % --- grab fixed budgets & settings from cfg (robust to missing fields) ---
    sh  = getf(cfg,'shared', struct());
    io  = getf(cfg,'io', struct());
    lm  = getf(cfg,'lowmem', struct());
    gr  = getf(cfg,'gray', struct());
    dd  = getf(cfg,'ddra', struct());

    save_tag  = getf(getf(cfg,'io',struct()), 'save_tag', 'run');
    n_m  = picknum(sh,'n_m');
    n_s  = picknum(sh,'n_s');
    n_k  = picknum(sh,'n_k');
    n_mv = picknum(sh,'n_m_val');
    n_sv = picknum(sh,'n_s_val');
    n_kv = picknum(sh,'n_k_val');
    D    = tryget_scalar(T,'D', picknum(getf(cfg,'sweep_grid',struct()),'D_list'));
    useN = resolve_use_noise(getf(cfg,'shared',struct()));

    % reduction text
    cap_stream = picknum(lm,'zonotopeOrder_cap', 50);
    cap_off    = picknum(getf(sh,'options_reach',struct()), 'zonotopeOrder', 100);

    % noise text (process noise only per current codebase)
    etaW  = picknum(dd,'eta_w', NaN);
    alpha = picknum(dd,'alpha_w', NaN);
    noise_str = 'no process noise (W=0)';
    if useN
        if isnan(etaW) || isnan(alpha)
            noise_str = 'process noise enabled (W\neq 0)';
        else
            noise_str = sprintf('process noise $W=\\langle 0,\\alpha_W G_W\\rangle$ with $\\eta_W=%g$, $\\alpha_W=%g$', etaW, alpha);
        end
    end

    % gray variant & ridge text
    gray_var = '(graySeq)';
    if isfield(gr,'methodsGray') && ~isempty(gr.methodsGray)
        gray_var = sprintf('(%s)', strjoin(string(gr.methodsGray), ','));
    end
    ridge_str = '';
    if isfield(dd,'allow_ridge') && dd.allow_ridge
        lam  = picknum(dd,'lambda', NaN);
        gam  = picknum(dd,'ridge_gamma', NaN);
        pol  = getf(dd,'ridge_policy','MAB');
        ridge_str = sprintf('DDRA ridge guard on ($\\lambda=%g$, $\\gamma=%g$, policy=%s); ', lam, gam, string(pol));
    end

    % axis labelling
    axis_label = pretty_axis_name(axisName);
    axis_vals_str = format_values(axisVals);

    % fixed budgets string (protect missing)
    budg = join_nonempty({kv('n_m',n_m), kv('n_s',n_s), kv('n_k',n_k), kv('D',D), ...
                          kv('n_m^{val}',n_mv), kv('n_s^{val}',n_sv), kv('n_k^{val}',n_kv)}, ', ');

    % reduction string
    red_str = sprintf('Girard reduction (stream cap=%g; offline order=%g)', cap_stream, cap_off);

    % seeds/PE 
    seed_str = 'rng(1), deterministic PE inputs';

    % support whether shape-aware columns exist
    have_shape = any(ismember(lower(T.Properties.VariableNames), ...
        {'dir_eps_med','dir_eps_p90','hout_p90','haus_sym_med','mw_gray_mean','mw_ddra_mean'}));

    % --- Build captions ---
    caps = struct();

    caps.fidcons = sprintf([ ...
        '\\textbf{Fidelity and conservatism vs %s.} ', ...
        'Grid: $%s=\\bigsqbra{%s}$. Fixed: %s. Methods: DDRA and RCSI/Gray %s. ', ...
        '%sNoise: %s. Reduction: %s. %s. ', ...
        'Fidelity uses interval-hull coverage (shape-aware/support metrics reported alongside). ', ...
        'Conservatism is aggregated output-interval width (and mean width where noted).'], ...
        axis_label, axis_label, axis_vals_str, budg, gray_var, ridge_str, noise_str, red_str, seed_str);

    caps.runtime = sprintf([ ...
        '\\textbf{Runtime vs %s.} ', ...
        'Total, learning, validation, and inference times vs $%s$. ', ...
        'Grid: $%s=\\bigsqbra{%s}$. Fixed: %s. Methods: DDRA and RCSI/Gray %s. ', ...
        '%sNoise: %s. Reduction: %s. %s.'], ...
        axis_label, axis_label, axis_label, axis_vals_str, budg, gray_var, ridge_str, noise_str, red_str, seed_str);

    if have_shape
        caps.shapeaware = sprintf([ ...
            '\\textbf{Shape-aware summaries vs %s.} ', ...
            'Directional support ratios (median/p90), one-sided support gap $H_{\\mathrm{out}}$ (p90), ', ...
            'and symmetric Hausdorff (median) vs $%s$. ', ...
            'Grid: $%s=\\bigsqbra{%s}$. Fixed: %s. Reduction: %s. %s.'], ...
            axis_label, axis_label, axis_label, axis_vals_str, budg, red_str, seed_str);
    else
        caps.shapeaware = '';
    end

    % per-step caption is generic
    caps.perstep = sprintf([ ...
        '\\textbf{Per-step diagnostics at a representative sweep point.} ', ...
        'Interval width (DDRA/Gray), coverage (%%), and Gray/True size ratio vs step $k$. ', ...
        'Same settings as main plot; see text for selected row and directions.']);

    % --- optionally write a TeX file with macros ---
    if writeTex
        if strlength(p.Results.texPath) > 0
            tex_out = p.Results.texPath;
        else
            tex_out = fullfile(fileparts(summary_csv), sprintf('%s_captions.tex', save_tag));
        end
        fid = fopen(tex_out,'w');
        fprintf(fid,'%% Auto-generated from %s\n', summary_csv);
        fprintf(fid,'\\newcommand{\\capFIDCONS}{%s}\n\n', escape_tex(caps.fidcons));
        fprintf(fid,'\\newcommand{\\capRUNTIME}{%s}\n\n',  escape_tex(caps.runtime));
        if ~isempty(caps.shapeaware)
            fprintf(fid,'\\newcommand{\\capSHAPE}{%s}\n\n', escape_tex(caps.shapeaware));
        end
        fprintf(fid,'\\newcommand{\\capPERSTEP}{%s}\n', escape_tex(caps.perstep));
        fclose(fid);
    end
end

% ----------------- helpers -----------------
function s = getf(S, f, d)
    if isstruct(S) && isfield(S,f), s = S.(f); else, s = d; end
end
function v = picknum(S, f, d)
    if nargin<3, d = NaN; end
    v = d;
    if ~isstruct(S) || ~isfield(S,f) || isempty(S.(f)), return; end
    x = S.(f);
    if isnumeric(x) && isscalar(x), v = x; return; end
    if isnumeric(x) && ~isscalar(x), v = x(1); return; end
    if iscell(x) && ~isempty(x) && isnumeric(x{1}), v = x{1}(1); end
end
function D = tryget_scalar(T, col, fallback)
    D = fallback;
    if istable(T) && any(strcmpi(T.Properties.VariableNames, col))
        v = unique(double(T.(col))); if ~isempty(v), D = v(1); end
    end
end
function s = format_values(v)
    v = v(:).';
    if isempty(v), s = ''; return; end
    if numel(v)>=3
        d = diff(v);
        if all(abs(d - d(1)) < 1e-12)
            s = sprintf('%g:%g:%g', v(1), d(1), v(end)); return;
        end
    end
    s = strjoin(compose('%g', v), ', ');
end
function lab = pretty_axis_name(ax)
    if isempty(ax), lab = 'axis'; return; end
    switch lower(ax)
        case 'n_m',    lab = 'n_m';
        case 'n_s',    lab = 'n_s';
        case 'n_k',    lab = 'n_k';
        case 'alpha_w',lab = '\alpha_W';
        case 'd',      lab = 'D';
        otherwise,     lab = ax;
    end
end
function s = kv(name, val)
    if isnan(val) || isempty(val), s = ''; else, s = sprintf('%s=%g', name, val); end
end
function s = join_nonempty(cellstrs, sep)
    cellstrs = cellstrs(~cellfun(@isempty, cellstrs));
    if isempty(cellstrs), s = '—'; else, s = strjoin(cellstrs, sep); end
end
function s = escape_tex(in)
    % Minimal escaping (leave math intact); assume input already uses $...$ where needed.
    s = strrep(in, '%', '\%');
    s = strrep(s, '&', '\&');
    s = strrep(s, '_', '\_');
end
function use_noise = resolve_use_noise(S)
    if isfield(S,'use_noise') && ~isempty(S.use_noise)
        use_noise = logical(S.use_noise); return;
    end
    g = ~isfield(S,'noise_for_gray') || logical(S.noise_for_gray);
    d = ~isfield(S,'noise_for_ddra') || logical(S.noise_for_ddra);
    use_noise = g && d;
end
function out = roundn(x, n)
    % round to 10^n; n = -12 useful for cleaning table floats
    out = round(x .* 10^(-n)) ./ 10^(-n);
end
