function C = colorscheme(style)
% colorscheme  Return a consistent color struct for plots.
%
% TUM blue (corporate): RGB(0,101,189) ~ #0065BD
% Complementary orange (balanced to TUM blue): ~ #E37222
% Grays: light/mid/dark + black

arguments
    style (1,:) char = 'tum'
end

switch lower(style)
    case 'tum'
        C.blue    = [0 101 189] / 255;   % TUM blue
        C.orange  = [227 114 34] / 255;  % warm orange complement
        C.grayL   = [230 230 230] / 255;
        C.grayM   = [160 160 160] / 255;
        C.grayD   = [80 80 80]   / 255;
        C.black   = [0 0 0];

        % Assign roles
        C.ddra = C.blue;    % main accent (your LaTeX tumblue)
        C.rcsi = C.orange;  % RCSI/Gray in orange (red-scale counterpart)
    otherwise
        error('Unknown scheme');
end
end
