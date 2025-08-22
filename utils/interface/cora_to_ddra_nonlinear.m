function [sysDisc, lookup] = cora_to_ddra_nonlinear(dynamics, type, overrides)
% cora_to_ddra_nonlinear
% Convert nonlinear dynamics from custom_loadDynamics() to the inputs that
% ddra_nonlinear(sysDisc, lookup) expects.
%
% Usage:
%   [sysDisc, lookup] = cora_to_ddra_nonlinear("lorenz","standard");
%   [R, R_data] = ddra_nonlinear(sysDisc, lookup);
%
% Notes:
%   - Works for CORA nonlinear classes (nonlinearSysDT / nonlinearDT / nonlinearARX).
%   - For linear systems, this will error and nudge you to use ddra_linear instead.
%   - You can pass overrides: any fields in 'overrides' are merged into 'lookup'.

    if nargin < 2 || isempty(type), type = "standard"; end
    if nargin < 3, overrides = struct(); end

    % 1) Load via your helper (as you already do elsewhere)
    [sys, R0, U, ~] = custom_loadDynamics(dynamics, type);

    % 2) Guard against linear models (use ddra_linear for those)
    sysClass = lower(class(sys));
    if ~contains(sysClass,'nonlinear') && ~contains(sysClass,'arx')
        error('cora_to_ddra_nonlinear:ExpectedNonlinear', ...
              ['System "%s" is not nonlinear (class: %s). ', ...
               'Use ddra_linear(...) for linear models.'], dynamics, class(sys));
    end

    % 3) Start from ddra_nonlinear defaults (if present), else inline defaults
    if exist('ddra_nonlinear_defaults','file')
        lookup = ddra_nonlinear_defaults();
    else
        lookup = struct();
        lookup.seed = 1; lookup.dt = 0.015; lookup.NN = 5;
        lookup.dim_x = []; lookup.dim_u = [];
        lookup.initpoints = 30; lookup.steps = 20;
        lookup.wfac = 1e-4; lookup.stepsLip = 1; lookup.initpointsLip = 50;
        lookup.normType = 2; lookup.addZeps = true;
        lookup.params.R0 = []; lookup.params.U = [];
        lookup.reach.zonotopeOrder = 100; lookup.reach.tensorOrder = 2; lookup.reach.errorOrder = 5;
        lookup.plot.do_plot = true; lookup.plot.dims = [1 2];
    end

    % 4) Fill required fields from loaded system & sets
    lookup.params.R0 = R0;
    lookup.params.U  = U;

    % Dimensions (robustly):
    %   - dim_x from R0
    %   - dim_u by sampling U once (works for cartProd as well)
    lookup.dim_x = length(center(R0));
    try
        u_probe = randPoint(U);
        lookup.dim_u = length(u_probe);
    catch
        % fallback (rare): assume U has center
        lookup.dim_u = length(center(U));
    end

    % Time step & fun handle (if available on the object)
    try lookup.dt = sys.dt;      catch, end
    try lookup.fun = sys.f;      catch, end  % ddra_nonlinear will ask for lookup.fun if needed

    % 5) Merge user overrides last
    lookup = local_merge(lookup, overrides);

    % 6) Return the system object as-is (reach_DT accepts CORA nonlinear types)
    sysDisc = sys;
end

% --------------------------- helpers -------------------------------------
function out = local_merge(a, b)
    out = a;
    if ~isstruct(b), return; end
    f = fieldnames(b);
    for i = 1:numel(f)
        k = f{i};
        if isfield(a,k) && isstruct(a.(k)) && isstruct(b.(k))
            out.(k) = local_merge(a.(k), b.(k));
        else
            out.(k) = b.(k);
        end
    end
end
