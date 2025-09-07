function augment_metrics(results_dir)
% AUGMENT_METRICS Compute additional metrics rom artifacts stored as .mat files
% and append to summary.csv's.
% Usage:
%   augment_metrics(fullfile(cfg.io.results_root, cfg.io.save_tag))

if nargin<1 || ~isfolder(results_dir)
    error('Provide results_dir containing summary.csv and artifacts/');
end
artdir = fullfile(results_dir, 'artifacts');
csv_sum = fullfile(results_dir, 'summary.csv');
csv_ps  = fullfile(results_dir, 'summary_perstep.csv');

% --- Load/prepare CSVs
S = readtable(csv_sum);
if ~ismember('row', S.Properties.VariableNames)
    S.row = (1:height(S))';    % stable row index = artifact number
end
if exist(csv_ps,'file'), PS = readtable(csv_ps); else, PS = table(); end
ps_has = @(v) ismember(v, PS.Properties.VariableNames);

% --- Enumerate artifacts
files = dir(fullfile(artdir, 'row_*.mat'));
if isempty(files)
    warning('No artifacts found at %s', artdir);
    return;
end
fprintf('Found %d artifact files.\n', numel(files));

% --- Back up old CSVs (once)
backup_once(csv_sum); backup_once(csv_ps);

% --- Parameters for E3/E4 directional sampling
L_dirs = 64;     % number of unit directions in output space
rng(42);         % deterministic sampling for reproducibility

for f = files'
    row_id = sscanf(f.name,'row_%d.mat');
    art = load(fullfile(f.folder,f.name));
    fprintf('Row %04d: computing E1–E7...\n', row_id);

    % ----- Common quantities
    sysG = art.sys_gray;     % linearSysDT (identified Gray)
    sysT = art.sys_ddra;     % linearSysDT (normalized "true")
    VAL  = art.VAL;          % struct with x0{b}, u{b}, y{b}, R0, U
    W_eff = art.W_eff;       % DDRA effective W actually used
    M_AB  = art.M_AB;        % matrix-zonotope model (DDRA)
    B = numel(VAL.x0);
    ny = size(VAL.y{1},2);
    nk = size(VAL.y{1},1);   % n_k_val (assumed consistent across blocks)
    Kred =  min(100,  getfield_def(art, {'meta','options_reach','zonotopeOrder'}, 100));

    % --- Build Gray disturbance override (same logic as run_sweeps)
    W_pred = choose_W_pred(sysG, sysT, VAL.U, W_eff);

    % =========================
    % E1: Per-step coverage, AUC
    % =========================
    [cov_gray_k, cov_ddra_k, fv_gray, fv_ddra] = ...
        perstep_coverage(sysG, sysT, VAL, M_AB, W_eff, W_pred, Kred);

    auc_cov_gray = mean(cov_gray_k);
    auc_cov_ddra = mean(cov_ddra_k);
    fv_med_gray  = median(fv_gray(~isinf(fv_gray)));
    fv_med_ddra  = median(fv_ddra(~isinf(fv_ddra)));

    % Append per-step coverages to perstep CSV image (merge-safe)
    PS = upsert_perstep(PS, row_id, (1:nk)', ...
        {'cov_gray','cov_ddra'}, [cov_gray_k(:), cov_ddra_k(:)]);

    % =========================
    % E2: Axis-wise width ratios (vs True)
    % =========================
    [waxis_gray_k_i, waxis_ddra_k_i] = axis_widths_vs_true(sysG, sysT, VAL, M_AB, W_eff, W_pred, Kred);
    % Summaries (median, 90th percentile across steps) for each output axis
    waxis_gray_med90 = summarize_axes(waxis_gray_k_i);
    waxis_ddra_med90 = summarize_axes(waxis_ddra_k_i);

    % Persist per-step axis widths (optional, can be large). We add only medians p90 to summary.csv
    % To persist per-step into perstep.csv, uncomment:
    % for j=1:ny
    %   PS = upsert_perstep(PS, row_id, (1:nk)', ...
    %      {sprintf('waxis_gray_%d',j), sprintf('waxis_ddra_%d',j)}, ...
    %      [waxis_gray_k_i(:,j), waxis_ddra_k_i(:,j)]);
    % end

    % =========================
    % E3/E4: Directional support ratios & Hausdorff (outer, sampled)
    % =========================
    dirs = sample_unit_dirs(ny, L_dirs);
    [supp_med_gray, supp_p90_gray, hdist_mean_gray, hdist_p90_gray] = ...
        support_and_hausdorff(sysG, sysT, VAL, W_pred, dirs, Kred);
    [supp_med_ddra, supp_p90_ddra, hdist_mean_ddra, hdist_p90_ddra] = ...
        support_and_hausdorff_DDRA(sysT, VAL, M_AB, W_eff, dirs, Kred);

    % =========================
    % E6: Nominal trajectory error (center dynamics)
    % =========================
    [enom_mean, enom_p90] = nominal_error(sysT, sysG, VAL);

    % =========================
    % E7: Ridge-induced widening at output (absolute & multiplier if base saved)
    % =========================
    [ridge_abs, ridge_mult] = ridge_widening(sysT, W_eff, art);

    % ----- Write into summary.csv (row-aligned)
    rid = S.row == row_id;
    if ~any(rid)
        warning('Row %d not found in summary.csv; skipping write.', row_id);
    else
        S = addcol(S, 'auc_cov_gray', auc_cov_gray, rid);
        S = addcol(S, 'auc_cov_ddra', auc_cov_ddra, rid);
        S = addcol(S, 'fv_med_gray',  fv_med_gray,  rid);
        S = addcol(S, 'fv_med_ddra',  fv_med_ddra,  rid);

        % axis summaries per output j
        for j=1:ny
            S = addcol(S, sprintf('waxis_gray_med_%d',j), waxis_gray_med90(j,1), rid);
            S = addcol(S, sprintf('waxis_gray_p90_%d',j), waxis_gray_med90(j,2), rid);
            S = addcol(S, sprintf('waxis_ddra_med_%d',j), waxis_ddra_med90(j,1), rid);
            S = addcol(S, sprintf('waxis_ddra_p90_%d',j), waxis_ddra_med90(j,2), rid);
        end

        % support/hausdorff summaries
        S = addcol(S, 'supp_med_gray',   supp_med_gray,   rid);
        S = addcol(S, 'supp_p90_gray',   supp_p90_gray,   rid);
        S = addcol(S, 'hdist_mean_gray', hdist_mean_gray, rid);
        S = addcol(S, 'hdist_p90_gray',  hdist_p90_gray,  rid);

        S = addcol(S, 'supp_med_ddra',   supp_med_ddra,   rid);
        S = addcol(S, 'supp_p90_ddra',   supp_p90_ddra,   rid);
        S = addcol(S, 'hdist_mean_ddra', hdist_mean_ddra, rid);
        S = addcol(S, 'hdist_p90_ddra',  hdist_p90_ddra,  rid);

        % nominal error and ridge widening
        S = addcol(S, 'enom_mean', enom_mean, rid);
        S = addcol(S, 'enom_p90',  enom_p90,  rid);
        S = addcol(S, 'ridge_output_abs',  ridge_abs,  rid);
        S = addcol(S, 'ridge_output_mult', ridge_mult, rid);
    end
end

% --- Save updated CSVs
writetable(S, csv_sum);
if ~isempty(PS)
    writetable(PS, csv_ps);
end
fprintf('Augmentation done.\n');
end

% --------- Helpers ----------
function backup_once(p)
if exist(p,'file') && ~exist([p '.bak'],'file')
    copyfile(p, [p '.bak']);
    fprintf('Backed up %s -> %s\n', p, [p '.bak']);
end
end

function S = addcol(S, name, val, rowsel)
if ~ismember(name, S.Properties.VariableNames)
    S.(name) = nan(height(S),1);
end
S.(name)(rowsel) = val;
end

function v = getfield_def(S, path, defaultVal)
% getfield with nested path = cellstr
v = defaultVal;
try
    t = S;
    for i=1:numel(path)
        if isstruct(t) && isfield(t, path{i})
            t = t.(path{i});
        else
            return;
        end
    end
    v = t;
catch
end
end

function W_pred = choose_W_pred(sysG, sysT, Uset, W_eff)
% Mirror disturbance policy used in run_sweeps
W_pred = [];
try
    if isprop(sysG,'nrOfDisturbances') && sysG.nrOfDisturbances>0
        if ~isempty(sysG.E) && isequal(size(sysG.E), size(sysG.B)) && norm(sysG.E-sysG.B,'fro') < 1e-12
            W_pred = Uset;   % E == B
        else
            W_pred = normalizeWForGray(sysG, W_eff);
        end
    end
catch
    % Leave empty if anything fails
end
end

function [cov_gray_k, cov_ddra_k, fv_gray, fv_ddra] = perstep_coverage(sysG, sysT, VAL, M_AB, W_eff, W_pred, Kred)
% Per-step containment (fraction of components contained), and FV per block.
B = numel(VAL.x0); ny = size(VAL.y{1},2); nk = size(VAL.y{1},1);

% Gray reaches per block
cov_gray_k_num = zeros(nk,1); cov_gray_k_den = B*ny*ones(nk,1);
FVg = inf(B,1);

for b=1:B
    params = struct();
    params.R0 = zero_center(VAL.R0) + VAL.x0{b};
    params.u  = VAL.u{b}';
    params.tFinal = sysG.dt*(nk-1);
    if isprop(sysG,'nrOfInputs')
        du = sysG.nrOfInputs - size(params.u,1);
        if du>0, params.u = [params.u; zeros(du,size(params.u,2))]; end
        if du<0, params.u = params.u(1:sysG.nrOfInputs,:); end
    end
    if isprop(sysG,'nrOfDisturbances') && sysG.nrOfDisturbances>0 && ~isempty(W_pred)
        params.W = coerceWToSys(sysG, W_pred);
    end
    % opt = struct(); Rg = reach(sysG, params, opt); Rg = Rg.timePoint.set;
    Rg = reach_with_options(sysG, params, Kred);

    hit_b = false;
    for k=1:nk
        Yg = safe_linearMap(Rg{k}, sysG.C);
        [lo,hi] = bounds_interval(Yg);
        y = VAL.y{b}(k,:)';
        inside = (y<=hi+1e-6) & (y>=lo-1e-6);
        cov_gray_k_num(k) = cov_gray_k_num(k) + sum(inside);
        if ~all(inside) && ~hit_b
            FVg(b) = k;
            hit_b = true;
        end
    end
end
cov_gray_k = 100*(cov_gray_k_num./cov_gray_k_den);
fv_gray = FVg;

% DDRA reaches per block (matrix-zonotope recursion)
cov_ddra_k_num = zeros(nk,1); cov_ddra_k_den = B*ny*ones(nk,1);
FVd = inf(B,1);

for b=1:B
    X = VAL.R0 + VAL.x0{b};
    for k=1:nk
        u_k = VAL.u{b}(k,:)';
        Yd = sysT.C * X + sysT.D * u_k;
        [lo,hi] = bounds_interval(Yd);
        y = VAL.y{b}(k,:)';
        inside = (y<=hi+1e-6) & (y>=lo-1e-6);
        cov_ddra_k_num(k) = cov_ddra_k_num(k) + sum(inside);
        if ~all(inside) && isinf(FVd(b)), FVd(b)=k; end

        % propagate to k+1
        U_pt = zonotope(u_k);
        if iscell(W_eff), Wc = W_eff{1}; else, Wc = W_eff; end
        X = M_AB * cartProd(X, U_pt) + Wc;
        X = reduce(X,'girard',Kred);
    end
end
cov_ddra_k = 100*(cov_ddra_k_num./cov_ddra_k_den);
fv_ddra = FVd;
end

function [waxis_gray_k_i, waxis_ddra_k_i] = axis_widths_vs_true(sysG, sysT, VAL, M_AB, W_eff, W_pred, Kred)
% Compute per-step, per-axis widths normalized by True widths (ratios).
B = numel(VAL.x0); ny = size(VAL.y{1},2); nk = size(VAL.y{1},1);

% Precompute True reaches per block (used for both Gray and DDRA)
Ytrue = cell(B,1);
for b=1:B
    paramsT = struct('R0', VAL.R0, 'U', VAL.U, 'u', VAL.u{b}', 'tFinal', sysT.dt*(nk-1));
    % Rt = reach(sysT, paramsT, struct()); Rt = Rt.timePoint.set;
    Rt = reach_with_options(sysT, paramsT, Kred);

    Ytrue{b} = cell(nk,1);
    for k=1:nk, Ytrue{b}{k} = safe_linearMap(Rt{k}, sysT.C); end
end

% Gray widths
waxis_gray_k_i = zeros(nk,ny);
for b=1:B
    paramsG = struct('R0', zero_center(VAL.R0)+VAL.x0{b}, 'u', VAL.u{b}', 'tFinal', sysG.dt*(nk-1));
    if isprop(sysG,'nrOfInputs')
        du = sysG.nrOfInputs - size(paramsG.u,1);
        if du>0, paramsG.u = [paramsG.u; zeros(du,size(paramsG.u,2))]; end
        if du<0, paramsG.u = paramsG.u(1:sysG.nrOfInputs,:); end
    end
    if isprop(sysG,'nrOfDisturbances') && sysG.nrOfDisturbances>0 && ~isempty(W_pred)
        paramsG.W = coerceWToSys(sysG, W_pred);
    end
    % Rg = reach(sysG, paramsG, struct()); Rg = Rg.timePoint.set;
    Rg = reach_with_options(sysG, paramsG, Kred);

    for k=1:nk
        Yg = safe_linearMap(Rg{k}, sysG.C);
        [loG,hiG] = bounds_interval(Yg);
        [loT,hiT] = bounds_interval(Ytrue{b}{k});
        denom = (hiT - loT) + 1e-12;
        waxis_gray_k_i(k,:) = waxis_gray_k_i(k,:) + ((hiG - loG)./denom).';
    end
end
waxis_gray_k_i = waxis_gray_k_i / B;   % average across blocks

% DDRA widths
waxis_ddra_k_i = zeros(nk,ny);
for b=1:B
    X = VAL.R0 + VAL.x0{b};
    for k=1:nk
        u_k = VAL.u{b}(k,:)';
        Yd = sysT.C * X + sysT.D * u_k;
        [loD,hiD] = bounds_interval(Yd);
        [loT,hiT] = bounds_interval(Ytrue{b}{k});
        denom = (hiT - loT) + 1e-12;
        waxis_ddra_k_i(k,:) = waxis_ddra_k_i(k,:) + ((hiD - loD)./denom).';

        U_pt = zonotope(u_k);
        if iscell(W_eff), Wc = W_eff{1}; else, Wc = W_eff; end
        X = M_AB * cartProd(X, U_pt) + Wc;
        X = reduce(X,'girard',Kred);
    end
end
waxis_ddra_k_i = waxis_ddra_k_i / B;
end

function stats = summarize_axes(waxis_k_i)
% Return [median, p90] per axis across steps
ny = size(waxis_k_i,2);
stats = nan(ny,2);
for j=1:ny
    v = waxis_k_i(:,j);
    stats(j,1) = median(v,'omitnan');
    stats(j,2) = prctile(v,90);
end
end

function dirs = sample_unit_dirs(n, L)
if n==1
    dirs = [1; -1]; return;
end
X = randn(n,L);
dirs = X ./ vecnorm(X);
end

function [supp_med, supp_p90, hmean, hp90] = support_and_hausdorff(sysG, sysT, VAL, W_pred, dirs, Kred)
% Gray vs True directional support ratios & outer Hausdorff (sampled).
B = numel(VAL.x0); nk = size(VAL.y{1},1); L = size(dirs,2);
eps_all = []; h_all = [];

for b=1:B
    % True reach
    paramsT = struct('R0', VAL.R0, 'U', VAL.U, 'u', VAL.u{b}', 'tFinal', sysT.dt*(nk-1));
    %Rt = reach(sysT, paramsT, struct()); Rt = Rt.timePoint.set;
    Rt = reach_with_options(sysT, paramsT, Kred);

    % Gray reach
    paramsG = struct('R0', zero_center(VAL.R0)+VAL.x0{b}, 'u', VAL.u{b}', 'tFinal', sysG.dt*(nk-1));
    if isprop(sysG,'nrOfInputs')
        du = sysG.nrOfInputs - size(paramsG.u,1);
        if du>0, paramsG.u = [paramsG.u; zeros(du,size(paramsG.u,2))]; end
        if du<0, paramsG.u = paramsG.u(1:sysG.nrOfInputs,:); end
    end
    if isprop(sysG,'nrOfDisturbances') && sysG.nrOfDisturbances>0 && ~isempty(W_pred)
        paramsG.W = coerceWToSys(sysG, W_pred);
    end
    % Rg = reach(sysG, paramsG, struct()); Rg = Rg.timePoint.set;
    Rg = reach_with_options(sysG, paramsG, Kred);


    for k=1:nk
        Yt = safe_linearMap(Rt{k}, sysT.C);
        Yg = safe_linearMap(Rg{k}, sysG.C);

        st = support_batch(Yt, dirs);
        sg = support_batch(Yg, dirs);
        eps = sg ./ (st + 1e-12);
        eps_all = [eps_all; eps(:)']; 

        hd = max(sg - st, 0);    % outer half only
        h_all = [h_all; hd(:)']; 
    end
end
% Aggregate across (k,b,dirs)
supp_med = median(eps_all,'all','omitnan');
supp_p90 = prctile(eps_all(:),90);
hmean    = mean(h_all,'all','omitnan');
hp90     = prctile(h_all(:),90);
end

function [supp_med, supp_p90, hmean, hp90] = support_and_hausdorff_DDRA(sysT, VAL, M_AB, W_eff, dirs, Kred)
% DDRA vs True directional supports & outer Hausdorff (sampled)
B = numel(VAL.x0); nk = size(VAL.y{1},1);
eps_all = []; h_all = [];

for b=1:B
    % True reach
    paramsT = struct('R0', VAL.R0, 'U', VAL.U, 'u', VAL.u{b}', 'tFinal', sysT.dt*(nk-1));
    % Rt = reach(sysT, paramsT, struct()); Rt = Rt.timePoint.set;
    Rt = reach_with_options(sysT, paramsT, Kred);

    % DDRA recursion
    X = VAL.R0 + VAL.x0{b};
    Yd = cell(nk,1);
    for k=1:nk
        u_k = VAL.u{b}(k,:)';
        Yd{k} = sysT.C * X + sysT.D * u_k;
        U_pt = zonotope(u_k);
        if iscell(W_eff), Wc = W_eff{1}; else, Wc = W_eff; end
        X = M_AB * cartProd(X, U_pt) + Wc;
        X = reduce(X,'girard',Kred);
    end

    for k=1:nk
        Yt = safe_linearMap(Rt{k}, sysT.C);
        Yx = Yd{k};
        st = support_batch(Yt, dirs);
        sx = support_batch(Yx, dirs);
        eps = sx ./ (st + 1e-12);
        eps_all = [eps_all; eps(:)']; 

        hd = max(sx - st, 0);
        h_all = [h_all; hd(:)']; 
    end
end
supp_med = median(eps_all,'all','omitnan');
supp_p90 = prctile(eps_all(:),90);
hmean    = mean(h_all,'all','omitnan');
hp90     = prctile(h_all(:),90);
end

function s = support_batch(Z, dirs)
% Zonotope support for multiple directions
try
    c = center(Z); G = generators(Z);
    s = (dirs.'*c) + sum(abs(G.'*dirs),1)';   % L×1
catch
    % Fallback: interval hull upper bound along dirs
    L = size(dirs,2); s = zeros(L,1);
    for l=1:L
        I = interval(linearMap(Z, dirs(:,l)')); % 1D set
        s(l) = supremum(I);
    end
end
end

function [lo,hi] = bounds_interval(Z)
I = interval(Z);
try
    lo = infimum(I); hi = supremum(I);
catch
    lo = I.inf; hi = I.sup;   % for very old CORA
end
lo = lo(:); hi = hi(:);
end

function Zc = zero_center(Z)
% Return zonotope with zero center but same generators
Zc = zonotope([zeros(size(center(Z),1),1), generators(Z)]);
end

function Y = safe_linearMap(X, C)
try
    Y = linearMap(X, C);
catch
    % If X is already in output space, return as-is
    Y = X;
end
end

function [enom_mean, enom_p90] = nominal_error(sysT, sysG, VAL)
% E6: l2 error between nominal trajectories (centers, no noise)
B = numel(VAL.x0); nk = size(VAL.y{1},1);
errs = zeros(B,nk);
for b=1:B
    xT = center(VAL.R0) + VAL.x0{b};
    xG = center(VAL.R0) + VAL.x0{b};     % same nominal IC
    yT = zeros(nk,size(VAL.y{b},2));
    yG = zeros(nk,size(VAL.y{b},2));
    for k=1:nk
        u = VAL.u{b}(k,:)';
        yT(k,:) = (sysT.C*xT + sysT.D*u).';
        yG(k,:) = (sysG.C*xG + sysG.D*u).';
        xT = sysT.A*xT + sysT.B*u;
        xG = sysG.A*xG + sysG.B*u;
    end
    errs(b,:) = vecnorm(yT - yG,2,2).';
end
enom_mean = mean(errs,'all');
enom_p90  = prctile(errs(:),90);
end

function [ridge_abs, ridge_mult] = ridge_widening(sysT, W_eff, art)
% E7: Absolute widening at output from W_eff; multiplier only if base saved.
try
    Yw = sysT.C * coerceWToSys(sysT, W_eff);
    [lo,hi] = bounds_interval(Yw);
    ridge_abs = sum(hi-lo);     % axis-sum width attributable to W_eff
catch
    ridge_abs = NaN;
end
if isfield(art,'W_base')
    try
        Yb = sysT.C * coerceWToSys(sysT, art.W_base);
        [lob,hib] = bounds_interval(Yb);
        ridge_mult = (sum(hi-lo)+1e-12)/(sum(hib-lob)+1e-12);
    catch
        ridge_mult = NaN;
    end
else
    ridge_mult = NaN; % not available unless you persist W_base
end
end

function T = upsert_perstep(T, row_id, kvec, names, data)
% Merge-add per-step columns for a given row/k
% names: cellstr of column names; data: numel(k) x numel(names)
if isempty(T)
    T = table(repmat(row_id,numel(kvec),1), kvec, 'VariableNames', {'row','k'});
end
% Ensure row/k key exists
[~,ia] = ismember([repmat(row_id,numel(kvec),1), kvec], [T.row, T.k], 'rows');
missing = find(ia==0);
for idx = missing(:)'
    % T = [T; table(row_id, kvec(idx), 'VariableNames', {'row','k'})]; 
    % --- BEGIN PATCH ---
    % Assume inputs: T (table already loaded), row_id (scalar double),
%                kvec (vector of ks for this row), and that we already
%                computed the set of missing ks.

% Normalize kvec to a column vector early
    kvec = kvec(:);
    
    % Existing ks for this row
    mask_row = ismember(T.row, row_id);                 % logical mask for this row
    k_exist  = T.k(mask_row);
    
    % Which ks are missing?
    k_missing = setdiff(kvec, k_exist, 'stable');       % column vector
    nnew = numel(k_missing);
    
    if nnew > 0
        % Create a typed empty table with identical schema as T
        vars   = T.Properties.VariableNames;
        vtypes = varfun(@class, T, 'OutputFormat', 'cell');
        newRows = table('Size', [nnew, numel(vars)], ...
                        'VariableTypes', vtypes, ...
                        'VariableNames', vars);
    
        % Initialize each variable with a reasonable default
        for j = 1:numel(vars)
            nm  = vars{j};
            cls = vtypes{j};
            switch cls
                case {'double','single'}
                    newRows.(nm)(:) = NaN;
                case 'logical'
                    newRows.(nm)(:) = false;
                case {'string','char'}
                    newRows.(nm)(:) = string(missing);
                case {'datetime','duration','calendarDuration'}
                    newRows.(nm)(:) = NaT;
                case 'categorical'
                    newRows.(nm) = categorical(repmat(missing, nnew, 1));
                otherwise
                    % Fallbacks: try missing, else empty cells
                    try
                        newRows.(nm)(:) = missing;
                    catch
                        newRows.(nm)(:) = {[]};
                    end
            end
        end
    
        % Assign row and k as columns (sizes must match height(newRows))
        newRows.row = repmat(row_id, nnew, 1);          % nnew×1
        newRows.k   = k_missing(:);                     % nnew×1
    
        % Append
        T = [T; newRows];
    end


    % --- END PATCH ---

end
% Fill columns
for j=1:numel(names)
    name = names{j};
    if ~ismember(name, T.Properties.VariableNames)
        T.(name) = nan(height(T),1);
    end
    % write values
    for t=1:numel(kvec)
        rid = (T.row==row_id) & (T.k==kvec(t));
        T.(name)(rid) = data(t,j);
    end
end
end

function Rsets = reach_with_options(sys, params, zonOrd)
% Wrapper to ensure CORA gets a valid options struct
if nargin < 3 || isempty(zonOrd), zonOrd = 100; end
opt = struct();
opt.zonotopeOrder = zonOrd;
% (optional, but common) opt.reductionTechnique = 'girard';
Robj = reach(sys, params, opt);
Rsets = Robj.timePoint.set;
end


