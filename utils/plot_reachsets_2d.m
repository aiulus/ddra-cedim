function plot_reachsets_2d(sys_gray, sys_ddra, VAL, W_eff, b, dims, varargin)
% plot_reachsets_2d  Side-by-side reachable OUTPUT sets for one VAL block.
    % See earlier version; only the R0 handling differs.
    
    p = inputParser;
    addParameter(p,'Colors', colorscheme('tum'));
    addParameter(p,'SaveDir', []);
    addParameter(p,'Name', sprintf('reachsets_b%d', b));
    addParameter(p,'TikZ', false);
    addParameter(p,'OptionsGrayReach', struct());
    parse(p,varargin{:});
    opt = p.Results;
    C = opt.Colors;
    
    x0   = VAL.x0{b};
    Ublk = VAL.u{b};      % (n_k × n_u)
    hasY = isfield(VAL,'y') && numel(VAL.y)>=b && ~isempty(VAL.y{b});
    if hasY, Yblk = VAL.y{b}; end
    
    A = sys_ddra.A; B = sys_ddra.B; 
    if isfield(VAL, 'C') && ~isempty(VAL.C)
        Cmat = VAL.C;
    else
        Cmat = sys_ddra.C;
    end
    
    if isfield(VAL, 'D') && ~isempty(VAL.D)
        Dmat = VAL.D;
    else
        Dmat = sys_ddra.D;
    end
    dt = sys_ddra.dt;
    
    [n_k, ~] = size(Ublk);
    n_y = size(Cmat,1);
    
    % ---------- DDRA side ----------
    f1 = figure('Color','w','Name',sprintf('DDRA reachsets (block %d)',b)); hold on; box on;
    
    % initial set for plotting
    if isfield(VAL,'R0') && ~isempty(VAL.R0)
        R0 = VAL.R0;
    else
        R0 = zonotope(zeros(size(x0)), 1e-9*eye(numel(x0))); % tiny ball, viz only
    end
    Xk = R0 + x0;
    
    for k = 1:n_k
        uk = Ublk(k,:)';
        Yk = affineMap(Xk, Cmat) + (Dmat * uk);
        plot(Yk, dims, 'FaceColor',C.ddra, 'FaceAlpha',0.12, 'EdgeColor','none');
        if k < n_k
            Xk = affineMap(Xk, A) + (B*uk) + W_eff;
        end
    end
    
    % nominal overlay
    if hasY
        plot(Yblk(:,dims(1)), Yblk(:,dims(2)), '-o','Color',C.grayD,'LineWidth',1.0,'MarkerSize',3);
    else
        x = x0; pts = zeros(n_k,n_y);
        for k=1:n_k
            pts(k,:) = (Cmat*x + Dmat*Ublk(k,:)')';
            if k<n_k, x = A*x + B*Ublk(k,:)'; end
        end
        plot(pts(:,dims(1)), pts(:,dims(2)), '-o','Color',C.grayD,'LineWidth',1.0,'MarkerSize',3);
    end
    xlabel(sprintf('y_{%d}',dims(1))); ylabel(sprintf('y_{%d}',dims(2)));
    title('DDRA reachable outputs'); axis equal; grid on;
    
    if ~isempty(opt.SaveDir)
        export_figure(f1, opt.SaveDir, [opt.Name '_ddra'], {'png','pdf'}, opt.TikZ);
    end
    
    % ---------- RCSI/Gray side ----------
    f2 = figure('Color','w','Name',sprintf('RCSI/Gray reachsets (block %d)',b)); hold on; box on;
    
    params = struct();
    params.R0 = zonotope(x0);           % exact initial
    params.u  = Ublk';                   % (n_u × n_k)
    params.tFinal = dt * (n_k-1);
    R = reach(sys_gray, params, opt.OptionsGrayReach);
    
    for k = 1:min(n_k, numel(R.timePoint.set))
        plot(R.timePoint.set{k}, dims, 'FaceColor',C.rcsi,'FaceAlpha',0.12,'EdgeColor','none');
    end
    
    if hasY
        plot(Yblk(:,dims(1)), Yblk(:,dims(2)), '-o','Color',C.grayD,'LineWidth',1.0,'MarkerSize',3);
    else
        plot(pts(:,dims(1)), pts(:,dims(2)), '-o','Color',C.grayD,'LineWidth',1.0,'MarkerSize',3);
    end
    xlabel(sprintf('y_{%d}',dims(1))); ylabel(sprintf('y_{%d}',dims(2)));
    title('RCSI/Gray reachable outputs'); axis equal; grid on;
    
    if ~isempty(opt.SaveDir)
        export_figure(f2, opt.SaveDir, [opt.Name '_gray'], {'png','pdf'}, opt.TikZ);
    end
end
