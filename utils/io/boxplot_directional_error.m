function boxplot_directional_error(T)
% BOX PLOT of directional support-ratio summaries (E2/E3)
% Works with either legacy names 'dir_eps_*' or new names 'supp_*'.

vnames = T.Properties.VariableNames;

% Pick first available among candidates
gm = pickvar(T, {'dir_eps_med_gray','supp_med_gray'});    % Gray median
dm = pickvar(T, {'dir_eps_med_ddra','supp_med_ddra'});    % DDRA median
gp = pickvar(T, {'dir_eps_p90_gray','supp_p90_gray'});    % Gray p90
dp = pickvar(T, {'dir_eps_p90_ddra','supp_p90_ddra'});    % DDRA p90

if isempty(gm) || isempty(dm)
    error(['E2 median columns for Gray/DDRA not found. Looked for any of:\n' ...
           ' Gray:  dir_eps_med_gray OR supp_med_gray\n' ...
           ' DDRA:  dir_eps_med_ddra OR supp_med_ddra\n' ...
           'Available columns:\n  %s'], strjoin(vnames, ', '));
end

% Convert to column vectors, drop NaNs
gm = gm(:); dm = dm(:); gp = gp(:); dp = dp(:);
gm = gm(~isnan(gm)); dm = dm(~isnan(dm));
gp = gp(~isnan(gp)); dp = dp(~isnan(dp));

figure('Name','Directional support ratios (median)');
boxplot([gm, dm], 'Labels', {'Gray (median)','DDRA (median)'});
ylabel('support(\hat Y)/support(Y_{true})');

if ~isempty(gp) && ~isempty(dp)
    figure('Name','Directional support ratios (p90)');
    boxplot([gp, dp], 'Labels', {'Gray (p90)','DDRA (p90)'});
    ylabel('support(\hat Y)/support(Y_{true})');
end
end

function x = pickvar(T, candidates)
x = [];
for i = 1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        x = T.(candidates{i});
        return
    end
end
end
