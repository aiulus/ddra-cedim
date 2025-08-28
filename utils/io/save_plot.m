function outfiles = save_plot(figOrAx, plots_dir, filename_base, varargin)
% SAVE_PLOT  Robust plot exporter for figures/axes/tiledlayout.
% Usage:
%   save_plot(gcf, plots_dir, 'runtime_panels_vs_nm');
%   save_plot(gca, plots_dir, 'myplot', 'Formats', {'png','pdf'}, 'Resolution', 200);

    % ---- args ----
    p = inputParser;
    addParameter(p, 'Formats', {'png'}, @(x) iscellstr(x) || isstring(x));
    addParameter(p, 'Resolution', 200, @(x) isnumeric(x) && isscalar(x));
    parse(p, varargin{:});
    fmts = string(p.Results.Formats);
    dpi  = p.Results.Resolution;

    if ~exist(plots_dir,'dir'), mkdir(plots_dir); end

    % ---- normalize handle (figure/axes/tiledlayout/numeric) ----
    target = [];
    if isnumeric(figOrAx)
        if isgraphics(figOrAx), target = figOrAx; end
    elseif isgraphics(figOrAx)
        target = figOrAx;
    end
    if isempty(target)
        warning('save_plot:badHandle','Invalid handle; using gcf().');
        target = gcf;
    end

    % Prefer exporting the container itself; keep a figure handle for print()
    if strcmp(get(target,'Type'),'axes')
        fig = ancestor(target,'figure');
    else
        fig = ancestor(target,'figure');
        if isempty(fig) && strcmp(get(target,'Type'),'figure')
            fig = target;
        elseif isempty(fig)
            fig = gcf;
        end
    end

    % ---- render & save ----
    drawnow; pause(0.01);  % give graphics a beat
    outfiles = strings(0);

    for fmt = fmts
        fmtc = lower(char(fmt));
        base = char(filename_base);
        outfile = fullfile(plots_dir, [base '.' fmtc]);

        try
            if verLessThan('matlab','9.8') % < R2020a: no exportgraphics
                % Device mapping for print()
                switch fmtc
                    case 'pdf', dev = '-dpdf'; opts = {'-painters'};
                    case 'eps', dev = '-depsc'; opts = {'-painters'};
                    otherwise,   dev = '-dpng';  opts = {['-r' num2str(dpi)]};
                end
                print(fig, outfile, dev, opts{:});
            else
                % exportgraphics path (container-aware)
                switch fmtc
                    case 'pdf'
                        exportgraphics(target, outfile, 'ContentType','vector', 'BackgroundColor','white');
                    case 'eps'
                        % exportgraphics -> EPS not supported; use print()
                        print(fig, outfile, '-depsc', '-painters');
                    otherwise % png, jpg, etc.
                        exportgraphics(target, outfile, 'Resolution', dpi, 'BackgroundColor','white');
                end
            end
            outfiles(end+1) = string(outfile); 
        catch ME
            warning('save_plot:exportFailed','Failed to save %s: %s. Trying print fallback.', outfile, ME.message);
            try
                % Format-consistent fallback
                switch fmtc
                    case 'pdf', print(fig, outfile, '-dpdf','-painters');
                    case 'eps', print(fig, outfile, '-depsc','-painters');
                    otherwise,   print(fig, outfile, '-dpng', ['-r' num2str(dpi)]);
                end
                outfiles(end+1) = string(outfile);
            catch ME2
                warning('save_plot:fallbackFailed','Fallback failed for %s: %s', outfile, ME2.message);
            end
        end
    end
end
