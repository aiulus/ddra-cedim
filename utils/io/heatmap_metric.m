function heatmap_metric(summary_csv, xfield, yfield, metric)
S = readtable(summary_csv);
[X,~,ix] = unique(S.(xfield));
[Y,~,iy] = unique(S.(yfield));
M = nan(numel(Y), numel(X));
for r = 1:height(S), M(iy(r), ix(r)) = S.(metric)(r); end
figure('Color','w'); h = heatmap(X, Y, M);
xlabel(xfield); ylabel(yfield); title(metric); colormap('parula');
end
