function [X_model, X_data] = ddra_linear_inference(sys, lookup, M_AB)
    %% Expects:
    %   * sys = c2d(.)
    %   * lookup = struct(...) 
    %   * M_AB : learned matrix zonotope of [A B]
    %
    % Returns:
    %   * X_model, X_data (exactly as the original function)

    %------------- BEGIN CODE --------------

    %% system dimensions
    dim_x = lookup.dim_x; 
    dim_u = lookup.dim_u; 
    dim_y = lookup.dim_y; 

    n_k_val = lookup.n_k_val; % length of the evaluation trajectory

    %% initial-state uncertainty X0
    eta_x = lookup.eta_x; 
    alpha_x = lookup.alpha_x; 
    G_x0  = diag(ones(dim_x, eta_x));
    c_x = lookup.c_x;
    c_delta_x = lookup.c_delta_x;
    c_x = c_x + c_delta_x;
    X0 = zonotope(c_x, alpha_x * G_x0);

    %% input uncertainty U
    eta_u = lookup.eta_u; 
    alpha_u = lookup.alpha_u; 
    G_u = diag(ones(dim_u, eta_u));
    c_u = lookup.c_u;
    c_delta_u = lookup.c_delta_u;
    c_u = c_u + c_delta_u;
    U = zonotope(c_u, alpha_u * G_u);

    %% process noise W
    c_w = lookup.c_w;
    alpha_w = lookup.alpha_w;
    eta_w = lookup.eta_w; 
    G_w = ones(dim_x, eta_w);
    W = zonotope(c_w, alpha_w * G_w);

    %% propagate sets (unchanged)
    X_model = cell(n_k_val+1,1);
    X_data = cell(n_k_val+1,1);
    X_model{1} = X0; 
    X_data{1} = X0;

    for i = 1:n_k_val
        X_model{i,1} = reduce(X_model{i,1}, 'girard', 400);
        X_model{i+1,1}= sys.A * X_model{i} + sys.B * U + W;

        X_data{i,1} = reduce(X_data{i,1}, 'girard', 400);
        X_data{i+1,1} = M_AB * (cartProd(X_data{i}, U)) + W;
    end

    %% visualization (unchanged)
    projectedDims = {[1 2],[3 4],[4 5]};
    aux_visualize_original(X0, X_model, X_data, projectedDims);

end

%--------------------------- HELPER FUNCTIONS ----------------------------

function aux_visualize_original(X0, X_model, X_data, projectedDims)
    axx{1} = [0.75,1.5,0.5,4]; 
    axx{2} = [0.75,3,0.8,2.2];
    axx{3} = [0.75,2.3,0.75,2.8];
    index=1;
    numberofplots = 5;%length(X_model)
    for plotRun=1:length(projectedDims)
        figure('Renderer', 'painters', 'Position', [10 10 700 900])
             
        index=index+1;
        % plot initial set
        handleX0 = plot(X0,projectedDims{plotRun},'k-','LineWidth',2);
        hold on;
       
        % plot reachable sets starting from index 2, since index 1 = X0
        
        % plot reachable sets from model
        for iSet=2:numberofplots
            handleModel=  plot(X_model{iSet},projectedDims{plotRun},'b','Filled',true,'FaceColor',[.8 .8 .8],'EdgeColor','b');
        end
        
        % plot reachable sets from data
        for iSet=2:numberofplots
            handleData=   plot(X_data{iSet},projectedDims{plotRun},'r');
        end
        
        % label plot
        xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
        ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
        %axis(axx{plotRun});
        % skip warning for extra legend entries
        warOrig = warning; warning('off','all');
        legend([handleX0,handleModel,handleData],...
            'Initial Set','Set from Model','Set from Data','Location','northwest');
        warning(warOrig);
        ax = gca;
        ax.FontSize = 22;
        %set(gcf, 'Position',  [50, 50, 800, 400])
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset;
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
        %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    end
end

function MZ = aux_buildMatrixZonotope(W, dim_x, totalsamples)
    index=1; GW = cell(1, size(W.G,2) * totalsamples);
    for i=1:size(W.G,2)
        Z = [W.center, W.generators];
        vec = Z(:,i+1);
        GW{index} = [vec, zeros(dim_x, totalsamples-1)];
        for j=1:totalsamples-1
            GW{j+index} = [GW{index+j-1}(:,2:end) GW{index+j-1}(:,1)];
        end
        index = j + index + 1;
    end
    MZ = matZonotope(zeros(dim_x, totalsamples), GW);
end
