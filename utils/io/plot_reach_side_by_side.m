function plot_reach_side_by_side(sys_ddra, R0, U, W, Xsets_ddra, configs_gray, dims, ttl, outpath)
% dims: [i j] output/state dims to show
% Xsets_ddra: cell {k=1..T} sets from DDRA inference (data-propagated)
% configs_gray: same struct list you already use with validateReach (true + gray)

if nargin<8, outpath=""; end

% Compute “true” LTI propagation (same recurrence used in DDRA’s X_model)
T = numel(Xsets_ddra);
X_true = cell(T,1); X_true{1} = R0;
for k=1:T-1
    X_true{k} = reduce(X_true{k}, 'girard', 400);
    X_true{k+1} = sys_ddra.A*X_true{k} + sys_ddra.B*U + W;
end

figure('Color','w'); tiledlayout(1,T, 'TileSpacing','tight','Padding','compact');

for k=1:T
    nexttile; hold on; grid on;
    if k<=numel(X_true),  plot(X_true{k}, dims, 'FaceColor',[.75 .75 .75], 'EdgeColor',[.5 .5 .5], 'DisplayName','True (LTI)'); end
    if k<=numel(Xsets_ddra), plot(Xsets_ddra{k}, dims, 'r', 'DisplayName','DDRA'); end

    % For Gray: show “true” config (idx=1) and the gray one (idx=2)
    if ~isempty(configs_gray) && iscell(configs_gray)
        % Recompute R at just this k via your validateReach (fast subset)
        % Or, if you already cached R, pass it in; here we just show the gray set if available:
        % (Assuming you called validateReach earlier and cached R{i}{k})
        % If not cached, you can skip this to avoid recomputation.
    end

    xlabel(sprintf('x_%d', dims(1))); ylabel(sprintf('x_%d', dims(2)));
    title(sprintf('%s | k=%d', ttl, k));
    legend('Location','best'); axis tight;
end

if ~isempty(outpath)
    exportgraphics(gcf, outpath, 'ContentType','vector');
end
end
