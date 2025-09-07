function stats = paired_compare(T, metric_stub)
    % Compare Gray–DDRA for any metric pair whose names include metric_stub
    % e.g., metric_stub = 'E1_AUC' (will search for *_gray and *_ddra)
    V = T.Properties.VariableNames;
    g = V(contains(V, metric_stub) & contains(V,'gray', 'IgnoreCase',true));
    d = V(contains(V, metric_stub) & contains(V,'ddra', 'IgnoreCase',true));
    assert(~isempty(g)&&~isempty(d),'Could not find matching Gray/DDRA columns for %s', metric_stub);

    x = T.(g{1}); y = T.(d{1});
    ok = ~isnan(x) & ~isnan(y);
    [p,~,st] = signrank(x(ok), y(ok));           % Wilcoxon signed-rank
    delta = cliffsDelta(x(ok), y(ok));           % custom below
    stats = struct('metric',metric_stub,'n',sum(ok),'median_delta',median(x(ok)-y(ok)), ...
                   'p_wilcoxon',p, 'cliffs_delta',delta);
    disp(stats);
end

function d = cliffsDelta(x,y)
    % Simple Cliff's delta (effect size)
    x = x(:); y = y(:); nx = numel(x); ny = numel(y);
    comp = 0;
    for i=1:nx, comp = comp + sum(x(i)>y) - sum(x(i)<y); end
    d = comp/(nx*ny);
end
