
function [R_true, R_data, out] = ddra_nonlinear(sys, lookup)
% ddra_nonlinear  Data-Driven Reachability for nonlinear DT systems
% Implements Algorithm 6 (Lipschitz Reachability) from Alanwar et al.
% Returns:
%   R_data : reachSet (CORA) with timePoint.set{k} = R'_k
%   R_true : (optional) model-based reach set; left empty unless requested
%   out    : struct with artifacts (data matrices, bounds, etc.)

% ---------- 0) Unpack  ----------
rand('seed',1);
dim_x = lookup.dim_x; dim_u = lookup.dim_u;

n_s = lookup.n_s; n_m = lookup.n_m; n_k = lookup.n_k; n_k_val = lookup.n_k_val;
totalsamples = n_m * n_s * n_k;

% X0
Gx = diag(ones(dim_x, lookup.eta_x));
X0 = zonotope(lookup.c_x + lookup.c_delta_x, lookup.alpha_x * Gx);

% U
Gu = diag(ones(dim_u, lookup.eta_u));
U  = zonotope(lookup.c_u + lookup.c_delta_u, lookup.alpha_u * Gu);

% W and its matrix zonotope Mw
Gw = ones(dim_x, lookup.eta_w);
W  = zonotope(lookup.c_w, lookup.alpha_w * Gw);
Mw = aux_buildMatrixZonotope(W, dim_x, totalsamples);

% dynamics handle (x_{k+1} = f(x_k,u_k))
if isfield(lookup,'fun') && ~isempty(lookup.fun)
    fun = lookup.fun;
else
    try, fun = sys.mFile; catch, error('Provide lookup.fun: x_next = fun(x,u).'); end
end

% ---------- 1) Generate dataset D = (U_-, X) ----------
u_all = zeros(dim_u, totalsamples);
for i=1:totalsamples, u_all(:,i) = randPoint(U); end

n_blocks = n_m * n_s;
x_all  = zeros(dim_x*n_blocks, n_k+1);
x_free = zeros(dim_x*n_blocks, n_k+1);
for b=1:n_blocks
    r = (b-1)*dim_x + 1;
    x0b = randPoint(X0);
    x_all(r:r+dim_x-1,1)  = x0b;
    x_free(r:r+dim_x-1,1) = x0b;
end

idx = 1;
for b=1:n_blocks
    rx = (b-1)*dim_x + 1;
    for t=1:n_k
        u_bt   = u_all(:,idx);
        x_curr = x_free(rx:rx+dim_x-1,t);
        x_nf   = fun(x_curr, u_bt);
        x_free(rx:rx+dim_x-1,t+1) = x_nf;
        x_all (rx:rx+dim_x-1,t+1) = x_nf + randPoint(W);
        idx = idx + 1;
    end
end

% Flatten blocks to (dim × totalsamples)
X_minus = zeros(dim_x, totalsamples);
X_plus  = zeros(dim_x, totalsamples);
U_minus = zeros(dim_u, totalsamples);
col = 1;
for b=1:n_blocks
    rx = (b-1)*dim_x + 1;
    for t=1:n_k
        X_minus(:,col) = x_all(rx:rx+dim_x-1, t);
        X_plus(:, col) = x_all(rx:rx+dim_x-1, t+1);
        U_minus(:,col) = u_all(:, col);
        col = col + 1;
    end
end

% ---------- 2) Lipschitz inflation Z_eps (L and gamma via helper) ----------
stepsLip      = getf(lookup,'stepsLip',1);
initpointsLip = getf(lookup,'initpointsLip',50);
normType      = getf(lookup,'normType',2);
addZeps       = getf(lookup,'addZeps',true);

[gamma,L] = aux_compLipConst(fun, U, X0, stepsLip, initpointsLip, dim_x, normType);
eps_vec   = (L(:) .* gamma(:)) / 2;               % per-coord L*delta/2
Z_eps     = zonotope([zeros(dim_x,1), diag(eps_vec)]);

% ---------- 3) Algorithm 6 propagation ----------
% As a reachSet object
R_cells = cell(n_k_val,1);
R_cells{1} = X0;

oneRow = ones(1, totalsamples);
C_Mw   = Mw.center;                   % center matrix of Mw

for k = 1:n_k_val-1
    % x*, u*  (linearization point)
    x_star = center(R_cells{k});
    u_star = center(U);

    % local LS model around (x*,u*)
    Phi = [ oneRow;
            X_minus - x_star .* ones(dim_x, totalsamples);
            U_minus - u_star .* ones(dim_u, totalsamples) ];
    Mprime = (X_plus - C_Mw) * pinv(Phi);

    % residuals across data columns
    Res = X_plus - Mprime * Phi;      % each col is r_j

    % interval hull across columns, then subtract process noise W
    lb = min(Res, [], 2);             % per-state lower
    ub = max(Res, [], 2);             % per-state upper
    %ZL = zonotope(interval(lb, ub)) - W;
    ZL = minkDiff(zonotope(interval(lb, ub)), W);

    % propagate: R_{k+1}’ = M’( {1}×(R_k’ - x*)×(U - u*) ) + W + ZL + Z_eps
    RK_shift = R_cells{k} + (-x_star);
    U_shift  = U + (-u_star);
    prodSet  = cartProd(1, cartProd(RK_shift, U_shift));  % {1}×...
    R_next   = Mprime * prodSet + W + ZL;
    if addZeps
        R_next = R_next + Z_eps;
    end

    % (Optional) light reduction to keep generator count sane
    R_cells{k+1} = reduce(R_next, 'girard', 200);
end

% wrap as CORA reachSet
tp.set  = R_cells(2:end);
tp.time = num2cell((sys.dt:sys.dt:sys.dt*(n_k_val-1))');
R_data  = reachSet(tp);

% ---------- 4) Model-based (optional) ----------
R_true = [];
if isfield(lookup,'compute_model_reach') && lookup.compute_model_reach
    % Try a standard model-based reach as a sanity check (best-effort)
    try
        params.R0 = X0; params.U = U; params.tFinal = sys.dt * n_k_val - sys.dt;
        opt.dim_x = dim_x; opt.zonotopeOrder = getf(lookup,{'reach','zonotopeOrder'},100);
        opt.tensorOrder = getf(lookup,{'reach','tensorOrder'},2); opt.errorOrder = getf(lookup,{'reach','errorOrder'},5);
        opt.W = W; opt.tStart = 0; opt.uTrans = 0;
        [R_true, ~] = reach_DT(sys, params, opt); 
    catch
        % leave R_true = []
    end
end

% ---------- 5) (Optional) plots ----------
if isfield(lookup,'plot_settings') && isfield(lookup.plot_settings,'projectedDims')
    aux_visualize( ...
        X0, R_true, R_data, ...
        lookup.plot_settings.projectedDims, ...
        getf(lookup.plot_settings,'numberofplots',5) );
end


% ---------- 6) Artifacts ----------
out = struct();
out.X_0T = X_minus; out.X_1T = X_plus; out.U_full = U_minus;
out.W = W; out.Mw = Mw; out.eps = eps_vec; out.L = L(:); out.gamma = gamma(:);
out.Mprime_last = exist('Mprime','var')*Mprime;

end

% ----------------- helpers -----------------
function aux_visualize(X0, R_true, R_data, projectedDims, numberofplots)
    % Make nonlinear visualization match the linear script’s style.
    % Plots X0 once, then first `numberofplots-1` time steps for model (filled gray)
    % and data (red outline) for each projection in projectedDims.

    % Extract stepwise sets from reachSets (if present)
    sets_true = {};
    if ~isempty(R_true) && isprop(R_true,'timePoint') && isfield(R_true.timePoint,'set')
        sets_true = R_true.timePoint.set;
    end
    sets_data = {};
    if ~isempty(R_data) && isprop(R_data,'timePoint') && isfield(R_data.timePoint,'set')
        sets_data = R_data.timePoint.set;
    end

    % How many steps to draw (cap by available)
    nT_true = numel(sets_true);
    nT_data = numel(sets_data);
    n_to_plot_model = min(numberofplots-1, nT_true);
    n_to_plot_data  = min(numberofplots-1, nT_data);

    % Iterate projections (like the linear helper)
    for p = 1:numel(projectedDims)
        dims = projectedDims{p};

        figure('Renderer','painters','Position',[10 10 700 900]);
        hold on; box on;

        % Initial set
        hX0 = plot(X0, dims, 'k-','LineWidth',2);

        % Model-based reachable sets (filled gray), if available
        hModel = [];
        for iSet = 1:n_to_plot_model
            hModel = plot(sets_true{iSet}, dims, 'b', ...
                          'Filled', true, 'FaceColor', [.8 .8 .8], 'EdgeColor', 'b');
        end

        % Data-driven reachable sets (red outline)
        hData = [];
        for iSet = 1:n_to_plot_data
            hData = plot(sets_data{iSet}, dims, 'r');
        end

        % Labels + legend
        xlabel(sprintf('x_{%d}',dims(1)));
        ylabel(sprintf('x_{%d}',dims(2)));

        % Build legend
        legHandles = hX0;
        legNames   = {'Initial Set'};
        if ~isempty(hModel), legHandles = [legHandles, hModel]; legNames{end+1} = 'Set from Model'; end
        if ~isempty(hData),  legHandles = [legHandles, hData];  legNames{end+1} = 'Set from Data';  end

        warOrig = warning; warning('off','all');
        legend(legHandles, legNames, 'Location','northwest');
        warning(warOrig);

        % Aesthetics 
        ax = gca; ax.FontSize = 22;
        outerpos = ax.OuterPosition; ti = ax.TightInset;
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width  = outerpos(3) - ti(1) - ti(3) - 0.01;
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
    end
end

function v = getf(S, key, def)
    if ischar(key) || isstring(key)
        v = def; if isfield(S,key), v = S.(key); end; return;
    end
    cur = S; v = def;
    for i=1:numel(key)
        if ~isfield(cur, key{i}), return; end
        cur = cur.(key{i});
    end
    v = cur;
end

function MZ = aux_buildMatrixZonotope(W, dim_x, totalsamples)
    % identical layout to linear case (cell-of-generators; OK in CORA ≥2024.2 with warnings)
    GW = cell(1, size(W.G,2) * totalsamples);
    index = 1;
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

% ============================ HELPERS =============================

function val = getfieldwithdef(S, path, def)
    % path can be 'name' or {'a','b','c'} for nested
    if ischar(path) || isstring(path)
        if isfield(S, path), val = S.(path); else, val = def; end
        return;
    end
    % nested struct path
    cur = S;
    for i=1:numel(path)
        k = path{i};
        if ~isfield(cur, k)
            val = def; return;
        end
        cur = cur.(k);
    end
    val = cur;
end

function [gamma, L] = aux_compLipConst(fun, U, X0, steps, initpoints, dim_x, normType)
    if nargin < 7 || isempty(normType), normType = 2; end
    totalsamples = initpoints * steps;

    % sample inputs
    u = zeros(length(center(U)), totalsamples);
    for i=1:totalsamples, u(:,i) = randPoint(U); end

    % roll out noise-free trajectories
    x_free = zeros(dim_x*initpoints, steps+1);
    idx = 1;
    for j=1:dim_x:initpoints*dim_x
        x_free(j:j+dim_x-1,1) = randPoint(X0);
        for i=1:steps
            x_free(j:j+dim_x-1, i+1) = fun(x_free(j:j+dim_x-1, i), u(:,idx));
            idx = idx + 1;
        end
    end

    % vectorize
    x0 = zeros(dim_x, totalsamples);
    x1 = zeros(dim_x, totalsamples);
    p0=1; p1=1;
    for j=1:dim_x:initpoints*dim_x
        for i=2:steps+1, x1(:,p1) = x_free(j:j+dim_x-1,i); p1=p1+1; end
        for i=1:steps,   x0(:,p0) = x_free(j:j+dim_x-1,i); p0=p0+1; end
    end

    % coordinate-wise Lipschitz bound
    L = zeros(dim_x,1); gamma = zeros(dim_x,1);
    for d = 1:dim_x
        bestL=0; bestG=0;
        for i=1:totalsamples
            z1 = [x0(d,i); u(:,i)]; f1 = x1(d,i);
            for j=1:totalsamples
                z2 = [x0(d,j); u(:,j)]; f2 = x1(d,j);
                if any(z1 ~= z2)
                    newL = norm(f1-f2, normType) / norm(z1-z2, normType);
                    newG = norm(z1-z2, normType);
                    if newL>bestL, bestL=newL; end
                    if newG>bestG, bestG=newG; end
                end
            end
        end
        L(d)=bestL; gamma(d)=bestG;
    end
end

function [R ,R_data]= reach_DT(obj,params,options,varargin)
% reach - computes the reachable sets of the discrete time system
%
% Syntax:  
%    R = reach_DT(obj,params,options)
%    [R,res] = reach_DT(obj,params,options,spec)
%
% Inputs:
%    obj - nonlinearSysDT object
%    params - parameter defining the reachability problem
%    options - options for the computation of the reachable set
%    spec - object of class specification 
%
% Outputs:
%    R - object of class reachSet storing the reachable set
%    res  - 1 if specifications are satisfied, 0 if not
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSysDT

% Author:       Matthias Althoff, Niklas Kochdumper, Amr Alanwar
% Written:      21-August-2012
% Last update:  29-January-2018
% Last revision:---

%------------- BEGIN CODE --------------

    % options preprocessing
    options = params2options(params,options);
    %options = checkOptionsReach(obj,options,0);
    
    spec = [];
    if nargin >= 4
       spec = varargin{1}; 
    end

    % compute symbolic derivatives
    derivatives(obj,options);

    % initialize cell array that stores the reachable sets
    t = options.tStart:obj.dt:options.tFinal;

    steps = length(t)-1;
    R = cell(steps+1,1);
    R_data = cell(steps+1,1);
    R{1} = params.R0;
    R_data{1} = params.R0;
    % loop over all reachablity steps
    for i = 1:steps

        % if a trajectory should be tracked
        if isfield(options,'uTransVec')
            options.uTrans = options.uTransVec(:,i);
        end  
        %reduce
        R{i} = reduce(R{i},'girard',20);
        R_data{i} = reduce(R_data{i},'girard',100);
        % compute next reachable set
        [R{i+1},R_data{i+1}] = linReach_DT(obj,R{i},R_data{i},options);

        if isfield(options,'verbose') && options.verbose 
            disp(t(i));
        end
        
        % check specification
        if ~isempty(spec)
           if ~check(spec,R{i+1})
               timePoint.set = R(2:i+1);
               timePoint.time = num2cell(t(2:i+1)');
               R = reachSet(timePoint);
               return;
           end
        end
    end

    % create reachable set object
    timePoint.set = R(2:end);
    timePoint.time = num2cell(t(2:end)');
    
    timePoint_data.set = R_data(2:end);
    timePoint_data.time = num2cell(t(2:end)');
    
    R = reachSet(timePoint);
    R_data = reachSet(timePoint_data);
end

%------------- END OF CODE --------------
function [Rtp ,Rtp_data] = linReach_DT(obj,Rinit,R_data,options)
% linReach - computes the reachable set after linearization
%
% Syntax:  
%    [Rtp] = linReach_DT(obj,Rinit,options)
%
% Inputs:
%    obj - nonlinearSysDT system object
%    Rinit - initial reachable set
%    options - options struct
%
% Outputs:
%    Rtp - resulting reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Amr Alanwar, Matthias Althoff, Niklas Kochdumper 
% Written:      21-August-2012
% Last update:  29-January-2018 (NK)
%               29-October-2020 (Amr) add data driven reachability
% Last revision:---

%------------- BEGIN CODE --------------

% linearize nonlinear system
[obj,A_lin,U] = linearize_DT(obj,Rinit,options); 

%translate Rinit by linearization point
Rdelta = Rinit + (-obj.linError.p.x);

% compute reachable set of linearized system
Rtp = A_lin*Rdelta + U;

% obtain linearization error
if options.tensorOrder > 2
    Verror = linError_thirdOrder(obj, options, Rdelta); 
else
    Verror = linError_mixed_noInt_DT(obj, options, Rdelta);   
end


%add interval of actual error
Rtp=Rtp+Verror+options.W;
% %%%%%%%%-------------------Data driven reachability-----------
options.Uorig= options.U +  options.uTrans;
 xStar = R_data.center;
 uStar =options.Uorig.center;
 xStarMat = repmat(xStar,1,size(options.X_0T,2));
 uStarMat = repmat(uStar,1,size(options.U_full,2));
 oneMat = repmat([1],1,size(options.U_full,2));
 IAB = (options.X_1T )*pinv([oneMat; options.X_0T+(-1*xStarMat);options.U_full+-1*uStarMat]);

V =  options.X_1T + -1*(IAB*[oneMat; options.X_0T+(-1*xStarMat);options.U_full+-1*uStarMat] + options.Wmatzono);
 VInt = intervalMatrix(V);
 leftLimit = VInt.Inf;
 rightLimit = VInt.Sup;
 
 V_one= zonotope(interval(min(leftLimit')',max(rightLimit')'));
 
 
 
 Rtp_data = IAB*cartProd([1],cartProd(R_data+(-1*xStar),options.Uorig+(-1*uStar))) +V_one+ options.W ;


end
%------------- END OF CODE --------------

function [obj,A_lin,U] = linearize_DT(obj,R,options)
% linearize - linearizes the nonlinearSysDT object
%
% Syntax:  
%    [obj,A_lin,U] = linearize(obj,R,options)
%
% Inputs:
%    obj - nonlinearSysDT system object
%    R - initial reachable set
%    options - options struct
%
% Outputs:
%    obj - nonlinearSysDT system object with additional properties
%    A_lin - system matrix of the linearized system
%    U - reachable set due to the inputs
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-August-2012
% Last update:  29-January-2018
% Last revision:---

%------------- BEGIN CODE --------------

%linearization point p.u of the input is the center of the input u
p.u = center(options.U) + options.uTrans;

%linearization point p.x and p.y
x0 = center(R);
p.x = x0;

%substitute p into the system equation in order to obtain the constant
%input
f0 = obj.mFile(p.x, p.u);

%get jacobian matrices
[A_lin,B_lin] = obj.jacobian(p.x, p.u);


uTrans = f0; %B*Ucenter from linOptions.U not added as the system is linearized around center(U)
Udelta = B_lin*(options.U+(-center(options.U)));
U = Udelta + uTrans;

%save linearization point
%% Commented out later
%obj.linError.p=p;

%------------- END OF CODE --------------
end