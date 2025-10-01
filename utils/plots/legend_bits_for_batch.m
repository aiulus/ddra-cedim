function [ddra_name, gray_name, hdr] = legend_bits_for_batch(cfg, sweep_grid, vary_key, rcsi_lbl)
% Build informative legend labels and a header summarizing fixed settings.
% vary_key ∈ {'n_k','n_m','alpha_w'}

    % --- DDRA variant (std/meas), robust to missing fields
    variant = "std";
    if isstruct(cfg) && isfield(cfg,'ddra') && isfield(cfg.ddra,'variant') && ~isempty(cfg.ddra.variant)
        variant = string(cfg.ddra.variant);
    end
    ddra_name = sprintf('DDRA-%s', char(variant));

    % --- Gray label from rcsi_lbl
    if nargin < 4 || isempty(rcsi_lbl), rcsi_lbl = "gray"; end
    gray_name = sprintf('RCSI-%s', char(rcsi_lbl));

    % --- Extract fixed values from sweep_grid (use NaN if vector/missing)
    D  = scalar_or_nan(getfield_safe(sweep_grid,'D_list'));
    nm = scalar_or_nan(getfield_safe(sweep_grid,'n_m_list'));
    ns = scalar_or_nan(getfield_safe(sweep_grid,'n_s_list'));
    nk = scalar_or_nan(getfield_safe(sweep_grid,'n_k_list'));
    aw = scalar_or_nan(getfield_safe(sweep_grid,'alpha_w_list'));

    % --- PE order L (if provided)
    L = NaN;
    pel = getfield_safe(sweep_grid,'pe_list');
    try
        if iscell(pel) && ~isempty(pel) && isstruct(pel{1}) && isfield(pel{1},'order') && ~isempty(pel{1}.order)
            L = pel{1}.order;
        elseif isstruct(pel) && ~isempty(pel) && isfield(pel(1),'order') && ~isempty(pel(1).order)
            L = pel(1).order;
        end
    catch
    end

    % --- Build tail string, omitting the swept axis
    kv = struct('D',D,'n_m',nm,'n_s',ns,'n_k',nk,'alpha_W',aw,'L',L);
    if isstring(vary_key) || ischar(vary_key), vk = char(vary_key); else, vk = ''; end
    switch vk
        case 'n_k',     rm = {'n_k'};
        case 'n_m',     rm = {'n_m'};
        case 'alpha_w', rm = {'alpha_W'};
        otherwise,      rm = {};
    end
    keys = {'D','n_m','n_s','n_k','alpha_W','L'};
    keys = keys(~ismember(keys, rm));

    parts = strings(0,1);
    for i = 1:numel(keys)
        val = kv.(keys{i});
        if ~isnan(val)
            if strcmp(keys{i}, 'alpha_W')
                parts(end+1) = sprintf('%s=%.3g', keys{i}, val);
            else
                parts(end+1) = sprintf('%s=%g', keys{i}, val);
            end
        end
    end
    tail = strjoin(parts, ', ');
    if strlength(tail)==0, tail = "fixed settings"; end

    ddra_name = sprintf('%s (%s)', ddra_name, tail);
    gray_name = sprintf('%s (%s)', gray_name, tail);

    % --- Header (system name + tail)
    sysname = "system";
    if isfield(cfg,'shared') && isfield(cfg.shared,'dyn') && ~isempty(cfg.shared.dyn)
        sysname = string(cfg.shared.dyn);
    end
    hdr = sprintf('%s — %s', char(sysname), char(tail));
end

% ------ helpers (local) ------
function v = getfield_safe(S, f)
    if isstruct(S) && isfield(S,f)
        v = S.(f);
    else
        v = NaN;
    end
end

function out = scalar_or_nan(v)
    if isnumeric(v) && isscalar(v)
        out = v;
    else
        out = NaN;
    end
end
