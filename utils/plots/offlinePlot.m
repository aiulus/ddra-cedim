function H = offlinePlot(csv_path, varargin)
% PLOT_SUMMARY_CSV  Flexible trend plots from a summary CSV.
%
% H = plot_summary_csv(csv_path, 'Name', Value, ...)
%
% Required
%   csv_path   : path to CSV (e.g., summary_agg.csv)
%
% Common options (Name/Value):
%   'x'        : string, x-axis column. Default: auto (alpha_w>n_m>n_k>n_s)
%   'y'        : string or cellstr of y columns to plot. Default:
%                tries {'cval_ddra','cval_gray'} then {'sizeI_ddra','sizeI_gray'}.
%   'group'    : string, column to split series by (e.g., 'ddra_variant','D').
%                You can combine multiple y's *and* a group; the legend will show both.
%   'filter'   : struct of equality/in-list filters, e.g., struct('D',2,'dyn',"k-Mass-SD")
%                Use numeric or string/scalar/cell arrays of allowed values.
%   'agg'      : 'median' or 'mean' aggregation over duplicate (x,group) rows. Default 'median'.
%   'err'      : [] (none), 'iqr' (25–75%), or 'p90' (5–95%). Default [].
%   'smooth'   : integer window (>=2) for moving average per series. Default [].
%   'xscale'   : 'linear' or 'log'. Default 'linear'.
%   'yscale'   : 'linear' or 'log'. Default 'linear'.
%   'title'    : plot title (string). Default auto.
%   'legendLoc': legend location. Default 'best'.
%   'save_dir' : folder to save. If empty, no saving. Default ''.
%   'filename' : base filename (no extension). Default auto from x/y/group.
%   'formats'  : cellstr of formats, e.g., {'png','pdf'}. Default {'png','pdf'}.
%   'dpi'      : resolution for raster formats. Default 200.
%   'linewidth': numeric. Default 1.8
%   'markers'  : true/false to show markers. Default true.
%
% Returns
%   H.fig, H.ax, H.lines (struct array), H.data (table used), H.series (table of series defs)
%
% Example
%   plot_summary_csv('summary_agg.csv', ...
%     'x','alpha_w', ...
%     'y',{'cval_ddra','cval_gray'}, ...
%     'group','ddra_variant', ...
%     'filter',struct('D',2,'dyn',"k-Mass-SD"), ...
%     'agg','median', 'err','iqr', ...
%     'save_dir','plots');
%

% ---------- parse inputs ----------
P = inputParser;
addParameter(P,'x',"",@(s)ischar(s)||isstring(s));
addParameter(P,'y',{},@(y)ischar(y)||isstring(y)||iscellstr(y)||iscellstr(string(y)));
addParameter(P,'group',"",@(s)ischar(s)||isstring(s));
addParameter(P,'filter',struct(),@isstruct);
addParameter(P,'agg','median',@(s)any(strcmpi(s,{'median','mean'})));
addParameter(P,'err',[],@(s)isempty(s)||any(strcmpi(s,{'iqr','p90'})));
addParameter(P,'smooth',[],@(x)isempty(x)||(isscalar(x) && x>=2));
addParameter(P,'xscale','linear',@(s)any(strcmpi(s,{'linear','log'})));
addParameter(P,'yscale','linear',@(s)any(strcmpi(s,{'linear','log'})));
addParameter(P,'title',"",@(s)ischar(s)||isstring(s));
addParameter(P,'legendLoc','best',@(s)ischar(s)||isstring(s));
addParameter(P,'save_dir','',@(s)ischar(s)||isstring(s));
addParameter(P,'filename','',@(s)ischar(s)||isstring(s));
addParameter(P,'formats',{'png','pdf'},@(c)iscellstr(c)||iscellstr(string(c)));
addParameter(P,'dpi',200,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(P,'linewidth',1.8,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(P,'markers',true,@(x)islogical(x)||ismember(x,[0,1]));
parse(P,varargin{:});
prm = P.Results;

% ---------- load & basic checks ----------
T = readtable(csv_path);
vars = string(T.Properties.VariableNames);

% utility lambdas
hascol = @(name) any(strcmpi(vars,string(name)));
col    = @(name) vars(strcmpi(vars,string(name)));

% auto-pick x if not provided
if strlength(prm.x)==0
    candidates = ["alpha_w","n_m","n_k","n_s","D"];
    prm.x = "";
    for c = candidates
        if hascol(c) && numel(unique(T.(col(c))))>1, prm.x = c; break; end
    end
    if strlength(prm.x)==0
        error('Could not auto-detect x-axis column. Please pass ''x'', e.g., ''alpha_w''.');
    end
end
xv = col(prm.x);

% auto-pick y if not provided
if isempty(prm.y)
    yTry = {["cval_ddra","cval_gray"], ["sizeI_ddra","sizeI_gray"], ...
            ["t_ddra_total","t_gray_total"], ["cov_auc_ddra","cov_auc_gray"]};
    ok = false;
    for k = 1:numel(yTry)
        yk = yTry{k};
        if all(arrayfun(@(s) hascol(s), yk))
            prm.y = cellstr(yk);
            ok = true; break;
        end
    end
    if ~ok
        % fallback: pick first numeric columns except x
        numCols = vars(varfun(@(x)isnumeric(x)&&~isscalar(x), T, 'OutputFormat','uniform'));
        numCols = setdiff(numCols, xv, 'stable');
        prm.y = cellstr(numCols(1:min(2,numel(numCols))));
    end
end
if ischar(prm.y) || isstring(prm.y), prm.y = cellstr(string(prm.y)); end

% ensure y columns exist
missingY = prm.y(~arrayfun(@(s) hascol(s), string(prm.y)));
if ~isempty(missingY)
    error('Missing y columns: %s', strjoin(missingY, ', '));
end

% apply simple filters
T = apply_filters(T, prm.filter);

% ensure x is numeric (coerce if needed)
T.(xv) = coerce_numeric(T.(xv));

% derive grouping key(s)
gcol = "";
if strlength(prm.group)>0
    if ~hascol(prm.group)
        error('Group column "%s" not found in CSV.', string(prm.group));
    end
    gcol = col(prm.group);
end

% collect long-form series: one series per (group value, y-name)
S = table();
si = 0;
if strlength(gcol)>0
    gvals = unique(T.(gcol));
else
    gvals = {[]}';  % single pseudo-group
end
for iy = 1:numel(prm.y)
    yname = col(prm.y{iy});
    for ig = 1:numel(gvals)
        si = si + 1;
        S.si(si,1) = si;
        S.yname(si,1) = string(yname);
        if strlength(gcol)>0
            gv = gvals(ig);
            % normalize group value to string for labeling
            if iscell(gv), gv = gv{1}; end
            try
                S.gval(si,1) = string(gv);
            catch
                S.gval(si,1) = string(num2str(gv));
            end
        else
            S.gval(si,1) = "";
        end
    end
end

% aggregate per series
aggfun = lower(string(prm.agg));
switch aggfun
    case "median", fagg = @(v) median(v,'omitnan');
    case "mean"  , fagg = @(v) mean(v,'omitnan');
end
switch lower(string(prm.err))
    case 'iqr', qlo = 0.25; qhi = 0.75;
    case 'p90', qlo = 0.05; qhi = 0.95;
    otherwise , qlo = [];   qhi = [];
end

% container for plotted data
Pdat = struct('series',[]);

figure('Color','w'); ax = axes; hold(ax,'on'); grid(ax,'on');

L = [];  % line handles
for iser = 1:height(S)
    yname = S.yname(iser);
    if strlength(gcol)>0
        mask = equal_match(T.(gcol), S.gval(iser));
    else
        mask = true(height(T),1);
    end
    Td = T(mask,:);
    if isempty(Td), continue; end

    % group by x (duplicates across seeds etc.)
    [G, xkeys] = findgroups(Td.(xv));
    yvals = Td.(yname);

    yAgg  = splitapply(@(v) fagg(v), yvals, G);
    if ~isempty(qlo)
        yLo = splitapply(@(v) qtile(v,qlo), yvals, G);
        yHi = splitapply(@(v) qtile(v,qhi), yvals, G);
    else
        yLo = []; yHi = [];
    end

    % sort by x
    [xSorted, ord] = sort(xkeys);
    yAgg = yAgg(ord);
    if ~isempty(yLo), yLo = yLo(ord); yHi = yHi(ord); end

    % optional smoothing
    if ~isempty(prm.smooth) && numel(xSorted)>=prm.smooth
        yAgg = movmean(yAgg, prm.smooth, 'omitnan');
        if ~isempty(yLo)
            yLo = movmean(yLo, prm.smooth, 'omitnan');
            yHi = movmean(yHi, prm.smooth, 'omitnan');
        end
    end

    % plot error band first (not in legend)
    if ~isempty(yLo)
        hh = fill([xSorted; flipud(xSorted)], [yLo; flipud(yHi)], [0 0 0], ...
            'FaceAlpha', 0.12, 'EdgeColor','none', 'Parent', ax);
        uistack(hh,'bottom');
    end

    % line + marker
    if prm.markers
        mk = pick_marker(iser);
    else
        mk = 'none';
    end
    L(iser) = plot(ax, xSorted, yAgg, ...
        'LineWidth', prm.linewidth, ...
        'Marker', mk, ...
        'DisplayName', legend_label(S, iser, gcol));

    % store data
    Pdat(iser).x = xSorted;
    Pdat(iser).y = yAgg;
    Pdat(iser).ylo = yLo;
    Pdat(iser).yhi = yHi;
    Pdat(iser).yname = yname;
    Pdat(iser).gval  = S.gval(iser);
end

% axes formatting
ax.XScale = lower(string(prm.xscale));
ax.YScale = lower(string(prm.yscale));
xlabel(ax, string(xv), 'Interpreter','none');
ylabel(ax, join(string(prm.y), ' / '), 'Interpreter','none');
if strlength(prm.title)>0
    title(ax, string(prm.title));
else
    ttl = sprintf('%s vs %s', strjoin(cellstr(prm.y), ', '), string(xv));
    if strlength(prm.group)>0, ttl = sprintf('%s   grouped by %s', ttl, string(gcol)); end
    title(ax, ttl, 'Interpreter','none');
end
legend(ax, 'Location', prm.legendLoc);

H.fig = gcf; H.ax = ax; H.lines = L; H.data = T; H.series = S;

% save if requested
if strlength(prm.save_dir)>0
    if ~exist(prm.save_dir,'dir'), mkdir(prm.save_dir); end
    if strlength(prm.filename)==0
        yslug = regexprep(strjoin(cellstr(prm.y),'_'), '[^\w-]','');
        if strlength(gcol)>0
            base = sprintf('trends_x_%s_y_%s_by_%s', string(xv), yslug, string(gcol));
        else
            base = sprintf('trends_x_%s_y_%s', string(xv), yslug);
        end
    else
        base = string(prm.filename);
    end
    for f = 1:numel(prm.formats)
        ext = lower(prm.formats{f});
        fpath = fullfile(prm.save_dir, sprintf('%s.%s', base, ext));
        switch ext
            case {'png','tif','tiff','jpg','jpeg','bmp'}
                print(H.fig, fpath, ['-d' ext], sprintf('-r%d', prm.dpi));
            case 'pdf'
                set(H.fig,'PaperPositionMode','auto');
                print(H.fig, fpath, '-dpdf', '-painters');
            case 'eps'
                print(H.fig, fpath, '-depsc', '-painters');
            otherwise
                warning('Unknown format "%s" — skipping.', ext);
        end
    end
end
end

% ================= helpers =================
function T2 = apply_filters(T, F)
    if isempty(fieldnames(F)), T2 = T; return; end
    mask = true(height(T),1);
    fn = string(fieldnames(F));
    for i = 1:numel(fn)
        f = fn(i);
        if ~ismember(f, string(T.Properties.VariableNames))
            warning('Filter field "%s" not in table; ignoring.', f); continue;
        end
        val = F.(f);
        col = T.(f);
        if iscell(val) || (isstring(val) && numel(val)>1) || (isnumeric(val) && numel(val)>1)
            m = false(height(T),1);
            for j = 1:numel(val)
                m = m | equal_match(col, val(j));
            end
            mask = mask & m;
        else
            mask = mask & equal_match(col, val);
        end
    end
    T2 = T(mask,:);
end

function tf = equal_match(col, v)
    if iscell(col), col = string(col); end
    if ischar(col), col = string(col); end
    if ischar(v) || isstring(v)
        tf = strcmpi(string(col), string(v));
    else
        if isdatetime(col)
            tf = col == v;
        else
            tf = (col == v);
        end
    end
end

function x = coerce_numeric(x)
    if isnumeric(x), return; end
    if iscell(x) || isstring(x) || ischar(x)
        try
            x = double(string(x));
        catch
            % fallback: categorical → double codes
            x = double(categorical(x));
        end
    end
end

function q = qtile(v, p)
    v = v(:); v = v(~isnan(v));
    if isempty(v), q = NaN; return; end
    try
        q = quantile(v, p);
    catch
        q = prctile(v, 100*p);
    end
end

function mk = pick_marker(k)
    mks = {'o','s','d','^','v','>','<','p','h','x','+'};
    mk = mks{ 1 + mod(k-1, numel(mks)) };
end

function s = legend_label(S, i, gcol)
    yn = string(S.yname(i));
    gv = string(S.gval(i));
    if strlength(gv)>0 && strlength(gcol)>0
        s = sprintf('%s | %s = %s', yn, string(gcol), gv);
    else
        s = yn;
    end
end
