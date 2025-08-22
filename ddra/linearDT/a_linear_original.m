function [X_model, X_data] = a_linear_original(sys, lookup)
    %% Expects 
    %   * sys = c2d(.)
    %   * lookup = struct(
    %
    %
    %% Original 
    % Author:       Amr Alanwar
    % Written:      28-October-2020
    % Last update:  
    % Last revision:---
    %
    %% Refactoring
    % Author:   Aybüke Ulusarslan
    % Written:  21-August-2025
    
    %------------- BEGIN CODE --------------
    
    rand('seed',1);
    
    %% system dynamics
    dim_x = lookup.dim_x; dim_u = lookup.dim_u; dim_y = lookup.dim_y;

    n_s = lookup.n_s; % Number of samples per unique trajectory
    n_m = lookup.n_m; % Number of unique test cases (x0, u_{1:n_k})
    n_k = lookup.n_k; % Trajectory length per test case
    totalsamples = n_m*n_k;
    
    
    %% Compute \calX_0
    eta_x = lookup.eta_x; % Number of generator vectors for \calX_0
    alpha_x = lookup.alpha_x; % Scaling factor for G_x0
    G_x0 = diag(ones(dim_x,eta_x)); % Generator template for \calX_0
    
    c_x =lookup.c_x;
    c_delta_x = lookup.c_delta_x;
    c_x = c_x + c_delta_x;
    
    X0 = zonotope(c_x, alpha_x*G_x0);
    
    %% Compute \calU_i
    eta_u = lookup.eta_u; % Number of generator vectors for \calU_i
    alpha_u = lookup.alpha_u; % Scaling factor for G_u
    G_u = diag(ones(dim_u, eta_u)); % Generator template for \calU_i
    
    c_u = lookup.c_u;
    c_delta_u = lookup.c_delta_u;
    c_u = c_u + c_delta_u;
    
    U = zonotope(c_u, alpha_u*G_u);
    
    %% Compute W
    c_w = lookup.c_w;
    alpha_w = lookup.alpha_w;
    eta_w = lookup.eta_w;
    G_w = ones(dim_x, eta_w); % Generator template for W
    W = zonotope(c_w, alpha_w*G_w);
    
    %Construct matrix zonotpe \mathcal{M}_w
    index=1;
    for i=1:size(W.G,2)
        Z = [W.center, W.generators];
        vec = Z(:,i+1);
        GW{index}= [ vec,zeros(dim_x,totalsamples-1)];
        for j=1:totalsamples-1
            GW{j+index}= [GW{index+j-1}(:,2:end) GW{index+j-1}(:,1)];
        end
        index = j+index+1;
    end
    W_mz = matZonotope(zeros(dim_x,totalsamples), GW);
    
    
    % randomly choose constant inputs for each step / sampling time
    u = zeros(totalsamples);
    for i=1:totalsamples
        u(i) = randPoint(U);
    end
    
    
    %simulate the system to get the data
    x0 = X0.center;
    x(:,1) = x0;
    utraj = zeros(dim_u, n_k);
    index=1;
    for j=1:dim_x:(n_m*dim_x)
        x(j:j+dim_x-1,1) = randPoint(X0);
        for i=1:n_k
            utraj(j,i) = u(index);
            x(j:j+dim_x-1,i+1) = sys.A*x(j:j+dim_x-1,i) + sys.B*u(index) + randPoint(W);      
            index=index+1;
        end
    end
    
    
    % concatenate the data trajectories 
    index_0 = 1;
    index_1 = 1;
    X_plus = zeros(dim_x, n_k);
    X_minus = zeros(dim_x, n_k);
    U_minus = zeros(dim_u, n_k);
    for j=1:dim_x:n_m*dim_x
        for i=2:n_k+1
            X_plus(:,index_1) = x(j:j+dim_x-1,i);
            index_1 = index_1 +1;
        end
        for i=1:n_k
            U_minus(:,index_0) = utraj(j,i);
            X_minus(:,index_0) = x(j:j+dim_x-1,i);
            index_0 = index_0 +1;
        end
    end
    
    U_full = U_minus(:,1:totalsamples);
    X_0T = X_minus(:,1:totalsamples);
    X_1T = X_plus(:,1:totalsamples);
    
    
    % plot simulated trajectory
    figure;
    subplot(1,2,1); hold on; box on; plot(x(1,:),x(2,:),'b'); xlabel('x_1'); ylabel('x_2');
    subplot(1,2,2); hold on; box on; plot(x(3,:),x(4,:),'b'); xlabel('x_3'); ylabel('x_4');
    close;
    
    X1W_cen =  X_1T - W_mz.center;
    X1W = matZonotope(X1W_cen,W_mz.generator);
    
    % set of A and B
    M_AB = X1W * pinv([X_0T;U_full]);
    
    % validate that A and B are within AB
    intAB11 = intervalMatrix(M_AB);
    intAB1 = intAB11.int;
    intAB1.sup >= [sys.A,sys.B]
    intAB1.inf <= [sys.A,sys.B]
    
    
    
    %% compute next step sets from model / data
    
    % set number of steps in analysis
    totalsteps = 5;
    X_model = cell(totalsteps+1,1);
    X_data = cell(totalsteps+1,1);
    % init sets for loop
    X_model{1} = X0; X_data{1} = X0;
    
    for i=1:totalsteps
        % 1) model-based computation
        X_model{i,1}=reduce(X_model{i,1},'girard',400);
        X_model{i+1,1} = sys.A * X_model{i} + sys.B * U + W;
        % 2) Data Driven approach
        X_data{i,1}=reduce(X_data{i,1},'girard',400);
        X_data{i+1,1} = M_AB * (cartProd(X_data{i},U)) + W; 
    end
    
    
    
    
    %% visualization
    projectedDims = {[1 2],[3 4],[4 5]};
    axx{1} = [0.75,1.5,0.5,4]; axx{2} = [0.75,3,0.8,2.2];axx{3} = [0.75,2.3,0.75,2.8];
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
    
    
    
    %------------- END OF CODE --------------


end

