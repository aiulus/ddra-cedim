function plot_fidelity_conservatism(summary_csv, perstep_csv)
    S = readtable(summary_csv);
    P = readtable(perstep_csv);
    
    % per-row conservatism (Gray vs DDRA width); per-row coverage AUCs
    G = groupsummary(P,'row','mean',{'wid_gray','wid_ddra','cov_gray','cov_ddra'});
    G.conserv_gray_over_ddra = G.mean_wid_gray ./ max(G.mean_wid_ddra, eps);
    
    J = outerjoin(S, G, 'Keys','row','MergeKeys',true);
    
    figure('Color','w'); tiledlayout(1,2,'TileSpacing','compact');
    
    nexttile; grid on; box on; hold on;
    scatter(J.cval_gray, J.conserv_gray_over_ddra, 24, 'filled');
    xlabel('Fidelity (Gray % on VAL)'); ylabel('Conservatism (Gray/DDRA width)');
    title('Fidelity vs Conservatism'); set(gca,'FontSize',11);
    
    nexttile; grid on; box on; hold on;
    scatter(J.mean_cov_ddra, J.mean_cov_gray, 24, 'filled');
    xlabel('DDRA coverage AUC_k'); ylabel('Gray coverage AUC_k');
    title('Coverage (AUC_k)'); axis([0 1 0 1]); refline(1,0); axis square; set(gca,'FontSize',11);
end
