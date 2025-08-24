function outfiles = save_plot(figOrAx, plots_dir, filename_base, varargin)
% SAVE_PLOT  Robust plot exporter that works with tiledlayout/figures.
% Usage:
%   save_plot(gcf, plots_dir, 'runtime_panels_vs_nm');
%   save_plot(gca, plots_dir, 'myplot', 'Formats', {'png','pdf'}, 'Resolution', 200);

    p = inputParser;
    addParameter(p, 'Formats', {'png'}, @(x) iscellstr(x) || isstring(x));
    addParameter(p, 'Resolution', 200, @(x) isnumeric(x) && isscalar(x));
    parse(p, varargin{:});
    fmts = string(p.Results.Formats);
    dpi  = p.Results.Resolution;

    if ~exist(plots_dir,'dir'), mkdir(plots_dir); end
    if ishghandle(figOrAx) && strcmp(get(figOrAx,'Type'),'axes')
        ax = figOrAx; fig = ancestor(ax,'figure');
    else
        fig = figOrAx; ax = [];
    end

    drawnow; pause(0.01);  % give tiledlayout a beat to render
    outfiles = strings(0);

    for fmt = fmts
        outfile = fullfile(plots_dir, filename_base + "." + fmt);
        try
            if verLessThan('matlab','9.8') % R2020a; exportgraphics arrived in R2020a
                switch fmt
                    case "png"
                        print(fig, outfile, '-dpng', ['-r' num2str(dpi)]);
                    case "pdf"
                        print(fig, outfile, '-dpdf', '-painters');
                    case "eps"
                        print(fig, outfile, '-depsc', '-painters');
                    otherwise
                        print(fig, outfile, '-dpng', ['-r' num2str(dpi)]);
                end
            else
                if isempty(ax)
                    exportgraphics(fig, outfile, 'Resolution', dpi, 'BackgroundColor','white');
                else
                    exportgraphics(ax,  outfile, 'Resolution', dpi, 'BackgroundColor','white');
                end
            end
            outfiles(end+1) = string(outfile); %#ok<AGROW>
        catch ME
            warning('save_plot:exportFailed','Failed to save %s: %s. Trying PRINT fallback.', outfile, ME.message);
            try
                print(fig, outfile, '-dpng', ['-r' num2str(dpi)]);
                outfiles(end+1) = string(outfile); %#ok<AGROW>
            catch ME2
                warning('save_plot:fallbackFailed','Fallback failed for %s: %s', outfile, ME2.message);
            end
        end
    end
end
