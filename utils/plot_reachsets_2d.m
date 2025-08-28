function plot_reachsets_2d(sys_gray, sys_ddra, VAL, W_eff, b, dims, varargin)
% plot_reachsets_2d  Side-by-side reachable OUTPUT sets for one VAL block.
%
% Usage:
%   plot_reachsets_2d(configs{2}.sys, sys_ddra, VAL, W_eff, 1, [1 2], ...
%                     'Colors', colorscheme('tum'), 'SaveDir','plots', 'TikZ', true)
%
% Inputs:
%   sys_gray : CORA linearSysDT (identified RCSI/Gray system)
%   sys_ddra : struct with fields A,B,C,D,dt (true or same model used for data)
%   VAL      : struct with fields x0{b}, u{b}, (optional) y{b}
%   W_eff    : zonotope disturbance used by DDRA (use zero zonotope if none)
%   b        : which VAL block to visualize (integer)
%   dims     : [i j] output dimensions to plot
%
% Options:
%   'Colors'   : struct from colorscheme('tum')
%   'SaveDir'  : folder (default = [])
%   'Name'     : base filename (default = 'reachsets_bX')
%   'TikZ'     : true/false (default = false)
%   'OptionsGrayReach' : struct reach options (default = sys defaults)

p = inputParser;
addParameter(p,'Colors', colorscheme('tum'));
addParameter(p,'SaveDir', []);
addParameter(p,'Name', sprintf('reachsets_b%d', b));
addParameter(p,'TikZ', false);
addParameter(p,'OptionsGrayReach', struct());
parse(p,varargin{:});
opt = p.Results;

C = opt.Colors;

% --- Extract VAL block
x0 = VAL.x0{b};               % (n_x × 1)
Ublk = VAL.u{b};              % (n_k × n_u)
hasY = isfield(VAL,'y') && numel(VAL.y)>=b && ~isempty(VAL.y{b});
if hasY
    Yblk = VAL.y{b};          % (n_k × n_y)
end

% --- Matrices
A = sys_ddra.A; B = sys_ddra.B;
if isfield(VAL,'C'), Cmat = VAL.C; else, Cmat = sys_ddra.C; end
if isfield(VAL,'D'), Dmat = VAL.D; else, Dmat = sys_ddra.D; end
dt = sys_ddra.dt;

[n_k, n_u] = size(Ublk);
n_y = size(Cmat,1);

% --- Left figure: DDRA forward propagation in output space
f1 = figure('Color','w','Name',sprintf('DDRA reachsets (block %d)',b));
hold on; box on;

% state set initial
Xk = VAL.R0_if_present; % not standard; fallback:
if ~isfield(VAL,'R0') || isempty(VAL.R0)
    % If you know your global R0, pass it via VAL.R0; otherwise assume small ball:
    R0 = zonotope(zeros(size(x0)), 1e-9*eye(numel(x0)));
else
    R0 = VAL.R0;
end
Xk = R0 + x0;  % shift initial set

for k = 1:n_k
    uk = Ublk(k,:)';                      % (n_u × 1)
    % Output reachable set at time k
    Yk = affineMap(Xk, Cmat) + (Dmat * uk);    % '+' → affine shift

    % Plot as filled patch (project dims)
    plot(Yk, dims, 'FaceColor',C.ddra, 'FaceAlpha',0.12, 'EdgeColor','none');

    % Next state set
    if k < n_k
        Xk = affineMap(Xk, A) + (B*uk) + W_eff;
    end
end

% Nominal trajectory (prefer measured Y if available)
if hasY
    plot(Yblk(:,dims(1)), Yblk(:,dims(2)), '-o', 'Color', C.grayD, 'LineWidth',1.0, ...
        'MarkerSize',3, 'DisplayName','nominal (measured)');
else
    % simulate nominal with w=0
    x = x0; pts = zeros(n_k,n_y);
    for k=1:n_k
        yk = Cmat*x + Dmat*Ublk(k,:)';
        pts(k,:) = yk';
        if k<n_k, x = A*x + B*Ublk(k,:)'; end
    end
    plot(pts(:,dims(1)), pts(:,dims(2)), '-o', 'Color', C.grayD, 'LineWidth',1.0, ...
        'MarkerSize',3, 'DisplayName','nominal (sim)');
end

xlabel(sprintf('y_{%d}',dims(1))); ylabel(sprintf('y_{%d}',dims(2)));
title('DDRA reachable outputs'); axis equal; grid on;

if ~isempty(opt.SaveDir)
    export_figure(f1, opt.SaveDir, [opt.Name '_ddra'], {'png','pdf'}, opt.TikZ);
end

% --- Right figure: RCSI/Gray reachable outputs (CORA reach)
f2 = figure('Color','w','Name',sprintf('RCSI/Gray reachsets (block %d)',b));
hold on; box on;

params = struct();
params.R0 = zonotope(x0);                 % exact initial
params.u  = Ublk';                        % (n_u × n_k)
params.tFinal = dt * (n_k-1);

opts = opt.OptionsGrayReach;
R = reach(sys_gray, params, opts);

for k = 1:min(n_k, numel(R.timePoint.set))
    Yk = R.timePoint.set{k};
    plot(Yk, dims, 'FaceColor', C.rcsi, 'FaceAlpha', 0.12, 'EdgeColor','none');
end

% Same nominal overlay for comparison
if hasY
    plot(Yblk(:,dims(1)), Yblk(:,dims(2)), '-o', 'Color', C.grayD, 'LineWidth',1.0, 'MarkerSize',3);
else
    % same simulated nominal as above (reuse)
    plot(pts(:,dims(1)), pts(:,dims(2)), '-o', 'Color', C.grayD, 'LineWidth',1.0, 'MarkerSize',3);
end

xlabel(sprintf('y_{%d}',dims(1))); ylabel(sprintf('y_{%d}',dims(2)));
title('RCSI/Gray reachable outputs'); axis equal; grid on;

if ~isempty(opt.SaveDir)
    export_figure(f2, opt.SaveDir, [opt.Name '_gray'], {'png','pdf'}, opt.TikZ);
end
end
