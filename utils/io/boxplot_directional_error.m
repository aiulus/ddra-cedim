function boxplot_directional_error(T)
% BOX PLOT of directional support-ratio summaries (E2/E3)
% Works with either legacy names 'dir_eps_*' or new names 'supp_*'.

% -- Fetch variables (first available among candidates)
gm = pickvar(T, {'dir_eps_med_gray','supp_med_gray'});    % Gray median
dm = pickvar(T, {'dir_eps_med_ddra','supp_med_ddra'});    % DDRA median
gp = pickvar(T, {'dir_eps_p90_gray','supp_p90_gray'});    % Gray p90
dp = pickvar(T, {'dir_eps_p90_ddra','supp_p90_ddra'});    % DDRA p90

if isempty(gm) || isempty(dm)
    vnames = T.Properties.VariableNames;
    error(['E2 median columns for Gray/DDRA not found. Looked for any of:\n' ...
           ' Gray:  dir_eps_med_gray OR supp_med_gray\n' ...
           ' DDRA:  dir_eps_med_ddra OR supp_med_ddra\n' ...
           'Available columns:\n  %s'], strjoin(vnames, ', '));
end

% -- Coerce to numeric column vectors
gm = double(gm(:)); dm = double(dm(:));

% -- Pairwise mask: same rows for both series
mask = isfinite(gm) & isfinite(dm);
Ymed = [gm(mask), dm(mask)];

if isempty(Ymed)
    warning('No finite paired data for medians; skipping median boxplot.'); 
else
    figure('Name','Directional support ratios (median)');
    boxplot(Ymed);                       % do NOT pass 'Labels' here
    set(gca, 'XTick', 1:2, ...
             'XTickLabel', {'Gray (median)','DDRA (median)'});
    ylabel('support(\hat{Y}) / support(Y_{true})');
    grid on
end

% p90 plot if available
if ~isempty(gp) && ~isempty(dp)
    gp = double(gp(:)); dp = double(dp(:));
    maskp = isfinite(gp) & isfinite(dp);
    Yp90 = [gp(maskp), dp(maskp)];

    if isempty(Yp90)
        warning('No finite paired data for p90; skipping p90 boxplot.');
    else
        figure('Name','Directional support ratios (p90)');
        boxplot(Yp90);                   % again, no 'Labels'
        set(gca, 'XTick', 1:2, ...
                 'XTickLabel', {'Gray (p90)','DDRA (p90)'});
        ylabel('support(\hat{Y}) / support(Y_{true})');
        grid on
    end
end
end

function x = pickvar(T, candidates)
% Return T.(name) for the first existing name in candidates; [] if none.
x = [];
for i = 1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        x = T.(candidates{i});
        return
    end
end
end
