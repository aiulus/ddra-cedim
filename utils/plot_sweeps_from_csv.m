function plot_sweeps_from_csv(csv_path, varargin)
% plot_sweeps_from_csv  Make runtime + fidelity/conservatism plots directly from summary.csv
%
% Usage:
%   plot_sweeps_from_csv('.../summary.csv');                          % on-screen only
%   plot_sweeps_from_csv('.../summary.csv','SaveDir','plots');        % export PNG/PDF
%   plot_sweeps_from_csv('.../summary.csv','TikZ',true);              % also export TikZ (if available)
%   plot_sweeps_from_csv('.../summary.csv','Which',{'runtime','fc'}); % select panels
%
% Options:
%   'SaveDir'  : output folder (default = [])
%   'TikZ'     : true/false (default = false); requires matlab2tikz on path if true
%   'Which'    : cellstr subset of {'runtime','fc'} (default = both)
%   'Formats'  : cellstr, e.g., {'png','pdf'} (default = {'png','pdf'})
%
% Auto x-axis:
%   prioritizes n_m > D > pe_order based on which varies in the CSV.

p = inputParser;
addParameter(p,'SaveDir',[]);
addParameter(p,'TikZ',false);
addParameter(p,'Which',{'runtime','fc'});
addParameter(p,'Formats',{'png','pdf'});
parse(p,varargin{:});
opt = p.Results;

T = readtable(csv_path);
C = colorscheme('tum');

% ---- coerce helpers
numify = @(v) double(string(v));  % works for numeric/char/string
hascol = @(n) ismember(n, T.Properties.VariableNames);

% Compute totals if missing
if hascol('t_ddra_learn') && hascol('t_ddra_check') && hascol('t_ddra_infer') && ~hascol('t_ddra_total')
    T.t_ddra_total = T.t_ddra_learn + T.t_ddra_check + T.t_ddra_infer;
end
if hascol('t_gray_learn') && hascol('t_gray_val')   && hascol('t_gray_infer') && ~hascol('t_gray_total')
    T.t_gray_total = T.t_gray_learn + T.t_gray_val  + T.t_gray_infer;
end

% ---- choose x-axis
cands = {'n_m','D','pe_order'};
vary = false(size(cands));
for i=1:numel(cands)
    if hascol(cands{i})
        xvals = T.(cands{i});
        vary(i) = numel(unique(xvals(~isnan(xvals)))) > 1;
    end
end
ix = find(vary,1,'first');
if isempty(ix) % fallback to n_m if present, else first numeric column
    if hascol('n_m'), xname='n_m'; else
        numcols = varfun(@isnumeric, T,'OutputFormat','uniform');
        xname = T.Properties.VariableNames{find(numcols,1,'first')};
    end
else
    xname = cands{ix};
end

% Sort by x for clean lines
[xs, idx] = sort(numify(T.(xname)));
T = T(idx,:);

% Labels
rcsi_name = 'RCSI/Gray';

% ===================== RUNTIME PANELS =====================
if any(strcmp(opt.Which,'runtime'))
    f = figure('Color','w','Name','Runtime panels'); 
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    plot(xs, T.t_ddra_total,'-o','Color',C.ddra,'LineWidth',1.6,'DisplayName','DDRA');
    plot(xs, T.t_gray_total,'-s','Color',C.rcsi,'LineWidth',1.6,'DisplayName',rcsi_name);
    xlabel(xname); ylabel('Seconds'); title('Total runtime'); grid on; legend('Location','best');

    nexttile; hold on;
    plot(xs, T.t_ddra_learn,'-o','Color',C.ddra,'LineWidth',1.6,'DisplayName','DDRA learn');
    plot(xs, T.t_gray_learn,'-s','Color',C.rcsi,'LineWidth',1.6,'DisplayName',[rcsi_name ' learn']);
    xlabel(xname); ylabel('Seconds'); title('Learning'); grid on; legend('Location','best');

    nexttile; hold on;
    plot(xs, T.t_ddra_check,'-o','Color',C.ddra,'LineWidth',1.6,'DisplayName','DDRA check');
    plot(xs, T.t_gray_val,  '-s','Color',C.rcsi,'LineWidth',1.6,'DisplayName',[rcsi_name ' validate']);
    xlabel(xname); ylabel('Seconds'); title('Validation/Check'); grid on; legend('Location','best');

    nexttile; hold on;
    plot(xs, T.t_ddra_infer,'-o','Color',C.ddra,'LineWidth',1.6,'DisplayName','DDRA infer');
    plot(xs, T.t_gray_infer,'-s','Color',C.rcsi,'LineWidth',1.6,'DisplayName',[rcsi_name ' infer']);
    xlabel(xname); ylabel('Seconds'); title('Inference'); grid on; legend('Location','best');

    if ~isempty(opt.SaveDir)
        export_figure(f, opt.SaveDir, ['runtime_panels_by_' xname], opt.Formats, opt.TikZ);
    end
end

% ===================== FIDELITY / CONSERVATISM =====================
if any(strcmp(opt.Which,'fc'))
    f2 = figure('Color','w','Name','Fidelity & Conservatism');
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % Fidelity
    nexttile; hold on;
    if hascol('cval_ddra')
        plot(xs, numify(T.cval_ddra),'-o','Color',C.ddra,'LineWidth',1.6,'DisplayName','DDRA');
    end
    if hascol('cval_gray')
        plot(xs, numify(T.cval_gray),'-s','Color',C.rcsi,'LineWidth',1.6,'DisplayName',rcsi_name);
    end
    xlabel(xname); ylabel('Containment on validation (%)');
    title('Fidelity'); grid on; legend('Location','best');

    % Conservatism
    nexttile; hold on;
    if hascol('sizeI_ddra')
        plot(xs, numify(T.sizeI_ddra),'-o','Color',C.ddra,'LineWidth',1.6,'DisplayName','DDRA');
    end
    if hascol('sizeI_gray')
        plot(xs, numify(T.sizeI_gray),'-s','Color',C.rcsi,'LineWidth',1.6,'DisplayName',rcsi_name);
    end
    xlabel(xname); ylabel('Aggregated interval size (proxy)');
    title('Conservatism'); grid on; legend('Location','best');

    if ~isempty(opt.SaveDir)
        export_figure(f2, opt.SaveDir, ['fidelity_conservatism_by_' xname], opt.Formats, opt.TikZ);
    end
end
end
