function export_figure(fig, outdir, basename, formats, tikz)
% export_figure  Save figure to PNG/PDF and optionally TikZ (if matlab2tikz available).
%
% export_figure(gcf,'plots','runtime_panels','png',true)

if ~exist(outdir,'dir'), mkdir(outdir); end

if nargin < 4 || isempty(formats), formats = {'png','pdf'}; end
if nargin < 5, tikz = false; end

% Raster/vector
for i=1:numel(formats)
    ext = lower(formats{i});
    fname = fullfile(outdir, [basename '.' ext]);
    switch ext
        case 'png'
            exportgraphics(fig, fname, 'Resolution', 200);
        case 'pdf'
            exportgraphics(fig, fname, 'ContentType','vector');
        otherwise
            warning('Unknown format: %s (skipped)', ext);
    end
end

% TikZ export if available
if tikz
    if exist('matlab2tikz','file') == 2
        fname = fullfile(outdir, [basename '.tex']);
        matlab2tikz(fname, 'standalone',false, 'height','\figureheight', 'width','\figurewidth');
    else
        warning('matlab2tikz not found on path; TikZ export skipped.');
    end
end
end
