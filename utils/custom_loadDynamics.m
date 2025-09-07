function [sys, R0, U, p_true] = custom_loadDynamics(dynamics, type, p)
% loadDynamics - load system dynamics and uncertainty sets
%
% Syntax:
%    [sys, R0, U] = loadDynamics(dynamics)
%
% Inputs:
%    dynamics - string specifying the dynamical system
%    type - string "standard" (default) for the normal uncertainty sets or
%               "diag" for diagonal generator matrices with random elements
%               "rand" for non-diagonal generator matrices
%    p - [optional] model parameters 
%
% Outputs:
%    sys - dynamical system
%    R0 - default initial state set
%    U - default input set
%    p_true - true model parameters
%
% References:
%   [1] L. Luetzow, M. Althoff, "Reachability analysis of ARMAX models," in 
%       Proc. of the 62nd IEEE Conference on Decision and Control, pp. 
%       7021–7028, 2023.
%   [2] E. N. Lorenz, “Deterministic nonperiodic flow,” Journal of
%       Atmospheric Sciences, vol. 20, no. 2, pp. 130 – 141, 1963.
%   [3] A. Kroll and H. Schulte, “Benchmark problems for nonlinear system
%       identification and control using soft computing methods: Need and
%       overview," Applied Soft Computing, vol. 25, pp. 496–513, 2014.
%   [4] J. M. Bravo. "Robust MPC of constrained discrete-time nonlinear 
%       systems based on approximated reachable sets", Automatica, 2006.
%   [5] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%       uncertain parameters using conservative linearization", in Proc. 
%       of the 62nd IEEE Conference on Decision and Control, pp.
%       4042-4048, 2008.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: testCase

% Authors:       Laura Luetzow
% Written:       01-March-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin == 1
    type = "standard";
end

switch dynamics
    case {"kLipMSD","k-Lip-MSD"}
        % Mass–Spring–Damper chain with onsite saturating (tanh) springs.
        % Nonlinear globally Lipschitz; scalable in D with tridiagonal coupling.
        %
        % Continuous-time model (about the linear equilibrium xbar):
        %   qdot = v
        %   M vdot = -Kc q - Cc v - Knl * tanh(Gamma * q) + u
        %
        % Here Kc is the usual chain stiffness (fixed ends), Cc ground damping.
        % The onsite term Knl*tanh(Gamma*q) is odd, saturating, and has bounded slope:
        %   d/dq [Knl*tanh(Gamma*q)] = Knl * diag( Gamma .* sech^2(Gamma.*q) )  <=  Knl*Gamma.
        % Hence the vector field is globally Lipschitz (linear part + bounded-slope nonlinearity).
        %
        % p: scalar D or struct with fields:
        %   p.k     : number of masses D (default 3)
        %   p.m     : masses [scalar or D-vector] (default 1)
        %   p.k_s   : linear springs k_0..k_D [scalar or (D+1)-vector] (default 1)
        %   p.c_d   : ground damping b_1..b_D [scalar or D-vector] (default 0.1)
        %   p.l     : natural lengths l_0..l_D [scalar or (D+1)-vector] (default 0)
        %   p.L     : right-wall position X_{D+1}=L (default 0)
        %   p.k_nl  : onsite sat-spring gains (diag entries of Knl) [scalar or D-vector] (default 0.5)
        %   p.gamma : sat slopes (diag entries of Gamma)             [scalar or D-vector] (default 1.0)
        %   p.dt    : sampling time (default 0.01)
        %
        % States: z = [q; v] in R^{2D}, with q = x - xbar (deflection about equilibrium)
        % Inputs: u in R^{D} (one force per mass)
        % Output: full state [q; v]

        if nargin < 3 || isempty(p), p = 3; end
        if isnumeric(p)
            D = round(p);
            p = struct('k',D,'m',1,'k_s',1,'c_d',0.1,'l',0,'L',0,...
                       'k_nl',0.5,'gamma',1.0,'dt',0.01);
        elseif isstruct(p)
            if ~isfield(p,'k'),     p.k = 3;      end
            if ~isfield(p,'m'),     p.m = 1;      end
            if ~isfield(p,'k_s'),   p.k_s = 1;    end
            if ~isfield(p,'c_d'),   p.c_d = 0.1;  end
            if ~isfield(p,'l'),     p.l = 0;      end
            if ~isfield(p,'L'),     p.L = 0;      end
            if ~isfield(p,'k_nl'),  p.k_nl = 0.5; end
            if ~isfield(p,'gamma'), p.gamma = 1.0;end
            if ~isfield(p,'dt'),    p.dt = 0.01;  end
        else
            error('Parameter p must be a scalar (k) or a struct.');
        end

        D   = p.k;   dt = p.dt;

        % --- broadcast scalars to vectors ---
        masses   = p.m;    if isscalar(masses),    masses    = masses   * ones(D,1);    else, masses = masses(:);    end
        cds = p.c_d;  if isscalar(cds),  cds  = cds * ones(D,1);    else, cds = cds(:);end
        ks  = p.k_s;  if isscalar(ks),   ks   = ks  * ones(D+1,1);  else, ks = ks(:);  end
        l   = p.l;    if isscalar(l),    l    = l   * ones(D+1,1);  else, l = l(:);    end
        Lw  = p.L;
        knl = p.k_nl; if isscalar(knl),  knl  = knl * ones(D,1);    else, knl = knl(:);end
        gam = p.gamma;if isscalar(gam),  gam  = gam * ones(D,1);    else, gam = gam(:);end

        % --- linear chain stiffness (fixed ends) ---
        main = ks(1:D) + ks(2:D+1);
        off  = -ks(2:D);
        Kc   = diag(main) + diag(off,1) + diag(off,-1);

        Cc = diag(cds);
        m_vec = p.m; if isscalar(m_vec), m_vec = m_vec*ones(D,1); else, m_vec = m_vec(:); end
        M = diag(m_vec);

        % --- preload from natural lengths and right wall; equilibrium shift q = x - xbar ---
        s = zeros(D,1);
        s(1) = ks(1)*l(1) - ks(2)*l(2);
        if D >= 3
            idx = 2:(D-1);
            s(idx) = ks(idx).*l(idx) - ks(idx+1).*l(idx+1);
        end
        s(D) = ks(D)*l(D) - ks(D+1)*l(D+1) - ks(D+1)*Lw;
        xbar = Kc \ s;

        % --- linear part (ZOH discretized) ---
        ZD = zeros(D); ID = eye(D);
        A_c = [ZD, ID; -(M\Kc), -(M\Cc)];
        B_c = [ZD;    M\ID];
        sysc = ss(A_c, B_c, eye(2*D), zeros(2*D,D));
        sysd = c2d(sysc, dt);

                % ---------- control augmentation (optional) ----------
        ctrl = struct(); if isfield(p,'ctrl'), ctrl = p.ctrl; end
        ID = eye(D);

        if isfield(ctrl,'S')
            S = ctrl.S;
        elseif isfield(ctrl,'idx')
            S = ID(:, unique(ctrl.idx(:)'));
        elseif isfield(ctrl,'arch') && strcmpi(ctrl.arch,'boundary') && D>=2
            S = ID(:, [1 D]);
        elseif isfield(ctrl,'arch') && strcmpi(ctrl.arch,'alternate')
            S = ID(:, 1:2:D);
        else
            S = ID;
        end

        if isfield(ctrl,'K'), K = ctrl.K; else, K = zeros(D,2*D); end

        Bk = sysd.B * K;     % (2D x 2D)
        Bv = sysd.B * S;     % (2D x m)

        % one-step map: (A + B K) z + B S v + saturated onsite term
        Knl   = diag(knl);
        Gamma = diag(gam);
        fun = @(z,v) (sysd.A + Bk)*z + Bv*v ...
                     + [zeros(D,1); dt * ( M \ ( -Knl * tanh(Gamma * z(1:D)) ) )];

        m_u = size(S,2);    
        dim_x = 2*D; dim_u = m_u;
        sys = nonlinearSysDT('kLipMSD', fun, dt, dim_x, dim_u);

        % ---------- sets (dim_u becomes m) ----------
        switch type
            case "standard"
                R0 = zonotope([zeros(dim_x,1), 0.1*eye(dim_x)]);
                U  = zonotope([zeros(dim_u,1), 0.2*eye(dim_u)]);
            case "diag"
                R0 = zonotope([0.05*randn(dim_x,1), 0.05*eye(dim_x)]);
                U  = zonotope([0.05*randn(dim_u,1), 0.1*eye(dim_u)]);
            case "rand"
                R0 = zonotope([0.1*randn(dim_x,1), 0.05*randn(dim_x, max(2, ceil(dim_x/6)))]);
                U  = zonotope([0.1*randn(dim_u,1), 0.05*randn(dim_u, max(2, ceil(dim_u/3)))]);
        end

        p_true = struct('m',m_vec,'k_s',ks,'c_d',cds,'l',l,'L',Lw,...
                        'k_nl',knl,'gamma',gam,'xbar',xbar,'dt',dt);
        p_true.ctrl = struct('K', K, 'S', S);
        p_true.m_act = size(S,2);

    case {"kDuffingMSD", "k- Duffing-MSD"}
        % MSD chain with onsite cubic springs: locally Lipschitz polynomial.
        % States z=[q;v] in R^{2D}, q = x - xbar. Inputs u: one force per mass.
        %
        % p: scalar D or struct:
        %   p.k     : #masses D (default 3)
        %   p.m     : masses [scalar or D-vector] (default 1)
        %   p.k_s   : linear springs k_0..k_D [scalar or (D+1)-vector] (default 1)
        %   p.c_d   : ground damping b_1..b_D [scalar or D-vector] (default 0.1)
        %   p.l     : natural lengths l_0..l_D (default 0)
        %   p.L     : right-wall position (default 0)
        %   p.gamma : onsite cubic coeffs Gamma_i >=0 [scalar or D-vector] (default 0.2)
        %   p.dt    : sampling time (default 0.01)

        if nargin < 3 || isempty(p), p = 3; end
        if isnumeric(p)
            D = round(p);
            p = struct('k',D,'m',1,'k_s',1,'c_d',0.1,'l',0,'L',0,'gamma',0.2,'dt',0.01);
        elseif isstruct(p)
            if ~isfield(p,'k'),     p.k = 3;     end
            if ~isfield(p,'m'),     p.m = 1;     end
            if ~isfield(p,'k_s'),   p.k_s = 1;   end
            if ~isfield(p,'c_d'),   p.c_d = 0.1; end
            if ~isfield(p,'l'),     p.l = 0;     end
            if ~isfield(p,'L'),     p.L = 0;     end
            if ~isfield(p,'gamma'), p.gamma = 0.2; end
            if ~isfield(p,'dt'),    p.dt = 0.01; end
        else
            error('Parameter p must be a scalar (k) or a struct.');
        end

        D  = p.k;          dt = p.dt;

        % Broadcast scalars to vectors
        masses   = p.m;    if isscalar(masses),    masses    = masses   * ones(D,1);   else, masses = masses(:);    end
        cds = p.c_d;  if isscalar(cds),  cds  = cds * ones(D,1);   else, cds = cds(:);end
        ks  = p.k_s;  if isscalar(ks),   ks   = ks  * ones(D+1,1); else, ks = ks(:);  end
        l   = p.l;    if isscalar(l),    l    = l   * ones(D+1,1); else, l = l(:);    end
        Lw  = p.L;
        gam = p.gamma; if isscalar(gam), gam  = gam * ones(D,1);   else, gam = gam(:);end

        % Linear chain stiffness (as in kMSD)
        main = ks(1:D) + ks(2:D+1);
        off  = -ks(2:D);
        Kc   = diag(main) + diag(off,1) + diag(off,-1);

        Cc = diag(cds);
        m_vec = p.m; if isscalar(m_vec), m_vec = m_vec*ones(D,1); else, m_vec = m_vec(:); end
        M = diag(m_vec);

        % Affine preload -> equilibrium shift
        s = zeros(D,1);
        s(1) = ks(1)*l(1) - ks(2)*l(2);
        if D >= 3
            idx = 2:(D-1);
            s(idx) = ks(idx).*l(idx) - ks(idx+1).*l(idx+1);
        end
        s(D) = ks(D)*l(D) - ks(D+1)*l(D+1) - ks(D+1)*Lw;
        xbar = Kc \ s;

        % Nonlinear Euler step
        % z = [q; v];  qdot = v;  vdot = M^{-1}(-Kc q - Cc v - Gamma q.^3 + u)
        Gamma = diag(gam);
                % ---------- control augmentation (optional) ----------
        ctrl = struct(); if isfield(p,'ctrl'), ctrl = p.ctrl; end
        ID = eye(D);

        if isfield(ctrl,'S')
            S = ctrl.S;
        elseif isfield(ctrl,'idx')
            S = ID(:, unique(ctrl.idx(:)'));
        elseif isfield(ctrl,'arch') && strcmpi(ctrl.arch,'boundary') && D>=2
            S = ID(:, [1 D]);
        elseif isfield(ctrl,'arch') && strcmpi(ctrl.arch,'alternate')
            S = ID(:, 1:2:D);
        else
            S = ID;
        end

        if isfield(ctrl,'K'), K = ctrl.K; else, K = zeros(D,2*D); end

        % z = [q; v];  qdot = v;
        % vdot = M^{-1}(-Kc q - Cc v - Gamma q.^3 + K z + S v_exo)
        fun = @(z,v_exo) [ ...
            z(1:D) + dt * z(D+1:end) ; ...
            z(D+1:end) + dt * ( M \ ( -Kc*z(1:D) - Cc*z(D+1:end) ...
                                     - Gamma*(z(1:D).^3) + K*z + S*v_exo ) ) ];
        
        m_u = size(S,2);    
        dim_x = 2*D; dim_u = m_u;
        sys = nonlinearSysDT('kDuffingMSD', fun, dt, dim_x, dim_u);

        % ---------- sets (dim_u becomes m) ----------
        switch type
            case "standard"
                R0 = zonotope([zeros(dim_x,1), 0.1*eye(dim_x)]);
                U  = zonotope([zeros(dim_u,1), 0.2*eye(dim_u)]);
            case "diag"
                R0 = zonotope([0.05*randn(dim_x,1), 0.05*eye(dim_x)]);
                U  = zonotope([0.05*randn(dim_u,1), 0.1*eye(dim_u)]);
            case "rand"
                R0 = zonotope([0.1*randn(dim_x,1), 0.05*randn(dim_x, max(2, ceil(dim_x/6)))]);
                U  = zonotope([0.1*randn(dim_u,1), 0.05*randn(dim_u, max(2, ceil(dim_u/3)))]);
        end

        p_true = struct('m',m_vec,'k_s',ks,'c_d',cds,'l',l,'L',Lw,'gamma',gam,'xbar',xbar,'dt',dt);
        p_true.ctrl = struct('K', K, 'S', S);
        p_true.m_act = size(S,2);

    % ===== NEW: globally Lipschitz 2D system (state-space) =====
    case "lipschitz2D"
        % x+ = x + dt*(A*x + g(x) + B*u) ,  g(x) = [tanh(k1*x1); tanh(k2*x2)]
        % globally Lipschitz since |tanh'| <= 1
        if nargin < 3 || ~isstruct(p), p = struct(); end
        dt = getfieldwithdef(p,'dt',0.05);
        k1 = getfielddef(p,'k1',0.8);
        k2 = getfieldwithdef(p,'k2',0.6);
        A  = getfieldwithdef(p,'A', [ -0.6  0.2; -0.1 -0.5 ]);
        B  = getfieldwithdef(p,'B', eye(2));

        dim_x = 2; dim_u = 2; dim_y = 2;
        fun = @(x,u) x + dt*(A*x + [tanh(k1*x(1)); tanh(k2*x(2))] + B*u);
        sys = nonlinearSysDT('lipschitz2D', fun, dt, dim_x, dim_u);

        % sets
        switch type
            case "rand"
                R0 = zonotope([ [0;0], 0.2*eye(dim_x) ]);
                U  = zonotope([ [0;0], 0.2*eye(dim_u) ]);
            case "diag"
                R0 = zonotope([ [0.1; -0.1], 0.1*eye(dim_x) ]);
                U  = zonotope([ [0.0;  0.0], 0.15*eye(dim_u) ]);
            otherwise % "standard"
                R0 = zonotope([ [0;0], 0.1*eye(dim_x) ]);
                U  = zonotope([ [0;0], 0.1*eye(dim_u) ]);
        end
        p_true = struct('dt',dt,'k1',k1,'k2',k2,'A',A,'B',B);

    % ===== NEW: globally Lipschitz 2D system (NARX) for black-box RCSI =====
    case "lipschitz2D_ARX"
        % y+ = A*y + [tanh(k1*y1); tanh(k2*y2)] + B*u   (p = 1)
        if nargin < 3 || ~isstruct(p), p = struct(); end
        dt = getfieldwithdef(p,'dt',0.05);
        k1 = getfieldwithdef(p,'k1',0.8);
        k2 = getfieldwithdef(p,'k2',0.6);
        A  = getfieldwithdef(p,'A', [ -0.6  0.2; -0.1 -0.5 ]);
        B  = getfieldwithdef(p,'B', eye(2));
        dim_y = 2; dim_u = 2; p_dim = 1;
        f = @(y,u) A*y(1:dim_y,1) + [tanh(k1*y(1,1)); tanh(k2*y(2,1))] + B*u(1:dim_u,1);
        sys = nonlinearARX('lipschitz2D_ARX', f, dt, dim_y, dim_u, p_dim);

        switch type
            case "rand"
                R0 = zonotope(zeros(dim_y*p_dim,1));                % ARX “state” is history
                U  = zonotope([ [0;0], 0.2*eye(dim_u) ]);
            case "diag"
                R0 = zonotope(zeros(dim_y*p_dim,1));
                U  = zonotope([ [0;0], 0.15*eye(dim_u) ]);
            otherwise
                R0 = zonotope(zeros(dim_y*p_dim,1));
                U  = zonotope([ [0;0], 0.1*eye(dim_u) ]);
        end
        p_true = struct('dt',dt,'k1',k1,'k2',k2,'A',A,'B',B);

    % ===== NEW: polynomial 2D system (state-space) =====
    case "poly2D"
        % x+ = x + dt*(A*x + [alpha*x1^2; beta*x1*x2] + B*u)
        if nargin < 3 || ~isstruct(p), p = struct(); end
        dt    = getfieldwithdef(p,'dt',0.05);
        alpha = getfieldwithdef(p,'alpha', 0.2);
        beta  = getfieldwithdef(p,'beta', -0.15);
        A     = getfieldwithdef(p,'A', [ -0.4  0.1; -0.2 -0.3 ]);
        B     = getfieldwithdef(p,'B', eye(2));
        dim_x = 2; dim_u = 2; dim_y = 2;
        fun = @(x,u) x + dt*( A*x + [alpha*x(1)^2; beta*x(1)*x(2)] + B*u );
        sys = nonlinearSysDT('poly2D', fun, dt, dim_x, dim_u);

        switch type
            case "rand"
                R0 = zonotope([ [0;0], 0.15*eye(dim_x) ]);
                U  = zonotope([ [0;0], 0.2*eye(dim_u) ]);
            case "diag"
                R0 = zonotope([ [0.1;-0.1], 0.1*eye(dim_x) ]);
                U  = zonotope([ [0;0], 0.15*eye(dim_u) ]);
            otherwise
                R0 = zonotope([ [0;0], 0.1*eye(dim_x) ]);
                U  = zonotope([ [0;0], 0.1*eye(dim_u) ]);
        end
        p_true = struct('dt',dt,'alpha',alpha,'beta',beta,'A',A,'B',B);

    % ===== NEW: polynomial 2D system (NARX) for black-box RCSI =====
    case "poly2D_ARX"
        % y+ = A*y + [alpha*y1^2; beta*y1*y2] + B*u   (p = 1)
        if nargin < 3 || ~isstruct(p), p = struct(); end
        dt    = getfieldwithdef(p,'dt',0.05);
        alpha = getfieldwithdef(p,'alpha', 0.2);
        beta  = getfieldwithdef(p,'beta', -0.15);
        A     = getfieldwithdef(p,'A', [ -0.4  0.1; -0.2 -0.3 ]);
        B     = getfieldwithdef(p,'B', eye(2));
        dim_y = 2; dim_u = 2; p_dim = 1;
        f = @(y,u) A*y(1:dim_y,1) + [alpha*y(1,1)^2; beta*y(1,1)*y(2,1)] + B*u(1:dim_u,1);
        sys = nonlinearARX('poly2D_ARX', f, dt, dim_y, dim_u, p_dim);

        switch type
            case "rand"
                R0 = zonotope(zeros(dim_y*p_dim,1));
                U  = zonotope([ [0;0], 0.2*eye(dim_u) ]);
            case "diag"
                R0 = zonotope(zeros(dim_y*p_dim,1));
                U  = zonotope([ [0;0], 0.15*eye(dim_u) ]);
            otherwise
                R0 = zonotope(zeros(dim_y*p_dim,1));
                U  = zonotope([ [0;0], 0.1*eye(dim_u) ]);
        end
        p_true = struct('dt',dt,'alpha',alpha,'beta',beta,'A',A,'B',B);

    case {"k-Mass-SD", "kMSD"}
        % Mass–Spring–Damper chain with D masses and fixed ends
        % Model matches the definition: ground damping b_i, natural lengths l_i,
        % left wall X_0 = 0, right wall X_{D+1} = L.
        %
        % p can be a scalar D or a struct with fields:
        %   p.k     : number of masses D (default 3)
        %   p.m     : mass(es) [scalar or D-vector] (default 1)
        %   p.k_s   : spring constants k_0..k_D [scalar or (D+1)-vector] (default 1)
        %   p.c_d   : ground damping b_1..b_D [scalar or D-vector] (default 0.1)
        %   p.l     : natural lengths l_0..l_D [scalar or (D+1)-vector] (default 0)
        %   p.L     : right-wall position X_{D+1}=L (default 0)
        %   p.dt    : sampling time (default 0.01)
        %
        % States: z = [q; v] in R^{2D}, with q = x - xbar (deflection)
        % Inputs: u in R^{D} (one force per mass)
        % Output: full state [q; v]

        % --- parse parameters ---
        if nargin < 3 || isempty(p), p = 3; end
        if isnumeric(p)
            D = round(p);
            p = struct('k', D, 'm', 1, 'k_s', 1, 'c_d', 0.1, 'l', 0, 'L', 0, 'dt', 0.01);
        elseif isstruct(p)
            if ~isfield(p,'k'),   p.k   = 3;    end
            if ~isfield(p,'m'),   p.m   = 1;    end
            if ~isfield(p,'k_s'), p.k_s = 1;    end
            if ~isfield(p,'c_d'), p.c_d = 0.1;  end
            if ~isfield(p,'l'),   p.l   = 0;    end
            if ~isfield(p,'L'),   p.L   = 0;    end
            if ~isfield(p,'dt'),  p.dt  = 0.01; end
        else
            error('Parameter p must be a scalar (k) or a struct.');
        end

        D   = p.k;
        dt  = p.dt;

        % Broadcast scalars to vectors
        masses   = p.m;   if isscalar(masses),   masses   = masses   * ones(D,1);     else, masses = masses(:);   end
        cds = p.c_d; if isscalar(cds), cds = cds * ones(D,1);     else, cds = cds(:); end
        ks  = p.k_s; if isscalar(ks),  ks  = ks  * ones(D+1,1);   else, ks  = ks(:);  end
        l   = p.l;   if isscalar(l),   l   = l   * ones(D+1,1);   else, l   = l(:);   end
        Lw  = p.L;

        % --- stiffness (fixed ends) ---
        % K = tridiag with diagonal ks(1:D) + ks(2:D+1), off-diagonals -ks(2:D)
        main = ks(1:D) + ks(2:D+1);
        off  = -ks(2:D);
        Kc   = diag(main) + diag(off,1) + diag(off,-1);

        % --- ground damping (diagonal) ---
        Cc = diag(cds);

        % --- mass matrix ---
        m_vec = p.m; if isscalar(m_vec), m_vec = m_vec*ones(D,1); else, m_vec = m_vec(:); end
        M = diag(m_vec);


        % --- affine load from natural lengths & right wall X_{D+1}=L ---
        s = zeros(D,1);
        % s_1 = k_0 l_0 - k_1 l_1
        s(1) = ks(1)*l(1) - ks(2)*l(2);
        % s_i = k_{i-1}l_{i-1} - k_i l_i, 2..D-1
        if D >= 3
            idx = 2:(D-1);
            s(idx) = ks(idx).*l(idx) - ks(idx+1).*l(idx+1);
        end
        % s_D = k_{D-1}l_{D-1} - k_D l_D - k_D * L
        s(D) = ks(D)*l(D) - ks(D+1)*l(D+1) - ks(D+1)*Lw;

        % --- equilibrium shift: Kc * xbar = s ---
        % (Kc is SPD for ks>0 with fixed ends)
        xbar = Kc \ s;

        % --- first-order LTI around q = x - xbar ---
        ZD = zeros(D); ID = eye(D);
        A_c = [ZD, ID; - (M \ Kc), - (M \ Cc)];
        B_c = [ZD;    M \ ID];
        C_c = eye(2*D);
        D_c = zeros(2*D, D);

        % --- discretize (ZOH) ---
        sysc = ss(A_c, B_c, C_c, D_c);
        sysd = c2d(sysc, dt);

        % ---------- control augmentation (optional) ----------
        % Build K (D x 2D) and S (D x m)
        ctrl = struct(); if isfield(p,'ctrl'), ctrl = p.ctrl; end
        ID = eye(D);

        % S selection
        if isfield(ctrl,'S')
            S = ctrl.S;
        elseif isfield(ctrl,'idx')
            S = ID(:, unique(ctrl.idx(:)'));
        elseif isfield(ctrl,'arch') && strcmpi(ctrl.arch,'boundary') && D>=2
            S = ID(:, [1 D]);
        elseif isfield(ctrl,'arch') && strcmpi(ctrl.arch,'alternate')
            S = ID(:, 1:2:D);
        else
            S = ID;  % full actuation
        end

        % K gain
        if isfield(ctrl,'K')
            K = ctrl.K;  % (D x 2D)
        else
            K = zeros(D, 2*D);
        end

        % Closed-loop A, exogenous B
        m_u = size(S,2);
        A_cl = sysd.A + sysd.B*K;
        B_ex = sysd.B*S;
        C_d  = eye(2*D);
        D_d  = zeros(2*D, m_u);
        sys  = linearSysDT(A_cl, B_ex, [], C_d, D_d, dt);
        % Make the disturbance map E equal to the input map B_ex
        %sys  = linearSysDT(A_cl, B_ex, B_ex, C_d, D_d, dt);


        % ---------- sets (dim_u becomes m) ----------   
        dim_x = 2*D; dim_u = m_u;
        switch type
            case "standard"
                R0 = zonotope([zeros(dim_x,1), 0.1*eye(dim_x)]);
                U  = zonotope([zeros(dim_u,1), 0.2*eye(dim_u)]);
            case "diag"
                R0 = zonotope([0.05*randn(dim_x,1), 0.05*eye(dim_x)]);
                U  = zonotope([0.05*randn(dim_u,1), 0.1*eye(dim_u)]);
            case "rand"
                R0 = zonotope([0.1*randn(dim_x,1), 0.05*randn(dim_x, max(2, ceil(dim_x/6)))]);
                U  = zonotope([0.1*randn(dim_u,1), 0.05*randn(dim_u, max(2, ceil(dim_u/3)))]);
        end

        % Expose true parameters (and shift) for downstream use
        p_true = struct('m',m_vec,'k_s',ks,'c_d',cds,'l',l,'L',Lw,'xbar',xbar);
        p_true.ctrl = struct('K', K, 'S', S);
        p_true.m_act = size(S,2);


    case "ddra5"
        A = [-1 -4  0  0  0;
              4 -1  0  0  0;
              0  0 -3  1  0;
              0  0 -1 -3  0;
              0  0  0  0 -2];
        B = ones(5,1);
        C = [1 0 0 0 0];
        D = 0;
        dt = 0.05;
    
        % Discretize via MATLAB c2d, then wrap as CORA linearSysDT
        sysc = ss(A,B,C,D);
        sysd = c2d(sysc, dt);
        sys = linearSysDT(sysd.A, sysd.B, [], sysd.C, sysd.D, dt);
    
        % Uncertainty sets 
        dim_x = 5; dim_u = 1;
        c_R0 = zeros(dim_x,1);                
        G_R0 = 0.1*eye(dim_x);                
        c_U = 1.0*ones(dim_u,1);              
        G_U = 0.25*eye(dim_u);                
        R0 = zonotope([c_R0, G_R0]);
        U = zonotope([c_U,  G_U]);
    
        p_true = [];  % not used here

        V = zonotope(zeros(size(C, 1),1));  % empty/degenerate V

    case "ddra5_ident"
        A = [-1 -4  0  0  0;
              4 -1  0  0  0;
              0  0 -3  1  0;
              0  0 -1 -3  0;
              0  0  0  0 -2];
        B = ones(5,1);
        C = ones(5,5);
        D = zeros(5,1);
        dt = 0.05;
    
        % Discretize via MATLAB c2d, then wrap as CORA linearSysDT
        sysc = ss(A,B,C,D);
        sysd = c2d(sysc, dt);
        sys = linearSysDT(sysd.A, sysd.B, [], sysd.C, sysd.D, dt);
    
        % Uncertainty sets 
        dim_x = 5; dim_u = 1;
        c_R0 = zeros(dim_x,1);                
        G_R0 = 0.1*eye(dim_x);                
        c_U = 1.0*ones(dim_u,1);              
        G_U = 0.25*eye(dim_u);                
        R0 = zonotope([c_R0, G_R0]);
        U = zonotope([c_U,  G_U]);
    
        p_true = [];  % not used here

        V = zonotope(zeros(size(C, 1),1));  % empty/degenerate V

    case "ddra5_2p"
        A = [-1 -4  0  0  0;
              4 -1  0  0  0;
              0  0 -3  1  0;
              0  0 -1 -3  0;
              0  0  0  0 -2];
        B = ones(5,1);
        C = [1 0 0 0 0;
             0 0 0 0 0];
        D = zeros(2,1);
        dt = 0.05;
    
        % Discretize via MATLAB c2d, then wrap as CORA linearSysDT
        sysc = ss(A,B,C,D);
        sysd = c2d(sysc, dt);
        sys = linearSysDT(sysd.A, sysd.B, [], sysd.C, sysd.D, dt);
    
        % Uncertainty sets 
        dim_x = 5; dim_u = 1;
        c_R0 = zeros(dim_x,1);                
        G_R0 = 0.1*eye(dim_x);                
        c_U = 1.0*ones(dim_u,1);              
        G_U = 0.25*eye(dim_u);                
        R0 = zonotope([c_R0, G_R0]);
        U = zonotope([c_U,  G_U]);
    
        p_true = [];  % not used here

        V = zonotope(zeros(size(C, 1),1));  % empty/degenerate V


    case "platoon"
        % Platoon dynamics (linear time-varying) as in CORA example
        % p can be a struct with fields:
        %   p.N_v : number of vehicles (default 2)
        %   p.N_k : horizon length used by platoonN (default 20)
        %
        % Note: 'type' is ignored for this case; no parameter vector p_true.
        if nargin < 3 || ~isstruct(p)
            p = struct();
        end
        if ~isfield(p,'n_n'); p.n_n = 2; end
        if ~isfield(p,'N_k'); p.N_k = 20; end

        dt = 0.5;              % discretization step
        n_n = p.n_n;            % #vehicles
        N_k = p.N_k;            % #time steps used to build LTV sys
        N_u = n_n;              % inputs = vehicles
        N_n = 3*n_n;            % states = 3 per vehicle

        % Build CORA LTV system
        sys = platoonN(dt, n_n, N_k);

        c_R0 = randn(N_n,1);
        for i = 0:n_n-1
            % enforce positive initial spacing for each vehicle’s spacing state
            c_R0(i*n_n + 2) = 3*abs(c_R0(i*n_n + 2));
        end
        alpha_R0 = 2*rand(N_n,1);
        c_U  = randn(N_u,1);
        alpha_U = rand(N_u,1);

        R0 = zonotope([c_R0, diag(alpha_R0)]);
        U  = zonotope([c_U,  diag(alpha_U)]);

        p_true = [];   % no gray-box parameter vector
    
    case "pedestrian"
        % pedestrian model as a state-space model [1]
        p_true = [1 0.01 5e-5 0.01]';
        if nargin < 3
            p = p_true;
        end
        A = [p(1)	0	    p(2)	0
            0	    p(1)	0	    p(2)
            0	    0	    p(1)	0
            0	    0	    0	    p(1)];
        B =[p(3)    0       0       0
            0	    p(3)    0       0
            p(4)    0       0       0
            0	    p(4)    0       0];
        C =[1	    0	    0	    0
            0	    1	    0	    0];
        D =[0	    0       1       0
            0	    0       0       1];

        dt = 0.01;
        sys = linearSysDT(A,B,[],C,D, dt);

        % create uncertainty sets
        dim_x = length(sys.A);
        dim_u = 2;
        dim_v = 2;
        switch type
            case "rand" % ---
                c_R0 = [-0.76; -9.68; 0.21; -5.42];
                c_U = [-0.16; -8.93];
                c_V = [1.48; -7.06];
                G_R0 = [-0.02 0.13 0.10 0.06
                    0.30 -0.24  0.21 -0.16
                    0.28  0.14  0.15  0.18
                    0.28  0.33 -0.06 -0.23];
                G_U = [0.07   -0.25
                    -0.28   -0.11];
                G_V = [-0.08    0.01
                    -0.00   -0.03];
            case "diag" % ---
                c_R0 = 0.1*[-0.76; -9.68; 0.21; -5.42];
                c_U = 0.1*[-0.16; -8.93];
                c_V = 0.1*[1.48; -7.06];
                G_R0 = diag([0.22 0.13 0.10 0.06]);
                G_U = diag([0.07 0.25]);
                G_V = diag([0.08 0.01]);
            case "standard" % ---
                c_R0 = zeros(dim_x,1);
                G_R0 = [];
                c_U = 0.1+zeros(dim_u,1);
                G_U = 0.2*diag(ones(dim_u,1));
                c_V =  -0.05+zeros(dim_v,1);
                G_V = 0.1*[diag(ones(dim_v,1)) ones(dim_v,1)];
        end
        R0 = zonotope([c_R0,G_R0]);
        V = zonotope([c_V,G_V]);
        U = cartProd(zonotope([c_U,G_U]), V);

    case "pedestrianARX"
        % pedestrian model as an ARX model [1]
        p_true = [2 -1 5e-5 -2]';
        if nargin < 3
            p = p_true;
        end
        A{1,1} = [  p(1)	0	    
                    0	    p(1)];
        A{2,1} = [  p(2)	0	    
                    0	    p(2)];
        B{1,1} = [  0	    0       1       0
                    0	    0       0       1];
        B{2,1} = [  p(3)    0       p(4)    0
                    0	    p(3)    0       p(4)];
        B{3,1} = [  p(3)    0       1       0
                    0	    p(3)    0       1];
        dt = 0.01;
        sys = linearARX(A, B, dt);

        % create uncertainty sets
        dim_x = 4;
        dim_u = 2;
        dim_v = 2;
        switch type
            case "rand" % ---
                c_U = [-0.16; -8.93];
                c_V = [1.48; -7.06];
                G_U = [0.07   -0.25
                    -0.28   -0.11];
                G_V = [-0.08    0.01
                    -0.00   -0.03];
            case "diag" % ---
                c_U = 0.1*[-0.16; -8.93];
                c_V = 0.1*[1.48; -7.06];
                G_U = diag([0.07 0.25]);
                G_V = diag([0.08 0.01]);
            case "standard" % ---
                c_U = 0.1+zeros(dim_u,1);
                G_U = 0.2*diag(ones(dim_u,1));
                c_V =  -0.05+zeros(dim_v,1);
                G_V = 0.1*[diag(ones(dim_v,1)) ones(dim_v,1)];
        end
        R0 = zonotope(zeros(dim_x,1));
        V = zonotope([c_V,G_V]);
        U = cartProd(zonotope([c_U,G_U]), V);            

    case "lorenz"
        % Lorenz system [2]
        p_true = [10 28 8/3]';
        if nargin < 3
            p = p_true;
        end
        dt = 0.01;
        fun = @(x,u) aux_dynLorenz(x,u,dt,p);
        dim_x = 3;
        dim_u = 3;
        dim_y = 2;
        out_fun = @(x,u) x(1:dim_y);
        sys = nonlinearSysDT('lorenz', fun, dt, dim_x, dim_u, out_fun, dim_y);

        switch type 
            case "rand" % ---
                c_R0 = [6.01; 9.36; -3.73];
                c_U = [7.56; -8.03; -1.57];
                G_R0 = [-0.11   -0.11    0.04
                    -0.07    0.08    0.14
                    0.02     0       0.02];
                G_U = [0.06   -0.02    0.04
                    -0.04    0.08   -0.09
                    -0.07    0.20   -0.10];
            case "diag" % ---
                c_R0 = 0.1*[6.01; 9.36; -3.73];
                c_U = 0.1*[7.56; -8.03; -1.57];
                G_R0 = 0.03*diag([0.11   0.11    0.24]);
                G_U = diag([0.06   -0.02    0.04]);
            case "standard" % ---
                c_U = [0.5;0.1;-0.2];
                c_R0 = [2; -1; 4];
                G_U =  diag([0.1;2;0.2]);
                G_R0 = 0.2*eye(dim_x);
        end
        R0 = zonotope([c_R0,G_R0]);
        W = zonotope([c_U,G_U]);
        V = zonotope([]);
        U = cartProd(W, V);
    
    case "lorenz_2D"
        % first two dimensions of the Lorenz system [2]
        p_true = [10 28]';
        if nargin < 3
            p = p_true;
        end
        dt = 0.01;
        fun = @(x,u) aux_dynLorenz2D(x,u,dt,p);
        dim_x = 2;
        dim_u = 2;
        dim_y = 2;
        out_fun = @(x,u) x(1:dim_y);
        sys = nonlinearSysDT('lorenz', fun, dt, dim_x, dim_u, out_fun, dim_y);

        switch type 
            case "rand" % ---
                c_R0 = [6.01; 9.36];
                c_U = [7.56; -8.03];
                G_R0 = [-0.11   -0.11
                    -0.07    0.08];
                G_U = [0.06   -0.02
                    -0.04    0.08];
            case "diag" % ---
                c_R0 = 0.1*[6.01; 9.36];
                c_U = 0.1*[7.56; -8.03];
                G_R0 = 0.03*diag([0.11   0.11]);
                G_U = diag([0.06   -0.02]);
            case "standard" % ---
                c_U = [0.5;0.1];
                c_R0 = [2; -1];
                G_U =  diag([0.1;2]);
                G_R0 = 0.2*eye(dim_x);
        end
        R0 = zonotope([c_R0,G_R0]);
        W = zonotope([c_U,G_U]);
        V = zonotope([]);
        U = cartProd(W, V);

    case "NARX" 
        % artificial NARX model, adapted from [3]
        p_true = [0.8 1.2]';
        if nargin < 3
            p = p_true;
        end

        f = @(y,u) [y(1,1)/(1+y(2,1)^2) + p(1)*u(3,1); ...
            (y(1,1) * y(2,1))/(1+y(2,1)^2)+ p(2)*u(6,1)];
        dt = 0.1;
        dim_y = 2;
        dim_u = 2;
        p_dim = 2;
        sys = nonlinearARX(dynamics,f,dt,dim_y, dim_u, p_dim);

        % initilization
        R0 = zonotope(zeros(dim_y*p_dim,1));

        % input
        switch type
            case "rand" % ---
                c_U = [-1.66; 4.41];
                G_U = [-0.1   0.13
                    0.25   -0.09];
            case "diag" % ---
                c_U = 0.1*[-1.66; 4.41];
                G_U = 0.7*diag([0.1   0.13]);
            case "standard" % ---
                c_U = [0;0.05];
                G_U =  0.2*eye(2);
        end
        U = zonotope(c_U, G_U);

    case "Square" 
        % artificial, simple NARX model
        p_true = [0.8 1.2]';
        if nargin < 3
            p = p_true;
        end

        f = @(y,u) [y(1,1)^2 + p(1)*u(3,1); ...
            y(2,1)^2+p(2)*u(2,1)];
        dt = 0.1;
        dim_y = 2;
        dim_u = 2;
        p_dim = 1;
        sys = nonlinearARX(dynamics,f,dt,dim_y, dim_u, p_dim);

        % initilization
        R0 = zonotope(zeros(dim_y*p_dim,1));

        % input
        switch type
            case "rand" % ---
                c_U = [-1.66; 4.41];
                G_U = [-0.1   0.13
                    0.25   -0.09];
            case "diag" % ---
                c_U = 0.1*[-1.66; 4.41];
                G_U = 0.7*diag([0.1   0.13]);
            case "standard" % ---
                c_U = [0;0.05];
                G_U =  0.2*eye(2);
        end
        U = zonotope(c_U, G_U);

    case "bicycle"
        % bicycle dynamics (see DOTBicycleDynamics_SRX_velEq.m)
        dt = 0.001;
        p_true = []; % no parameters defined
        fun = @(x,u) x + dt*DOTBicycleDynamics_SRX_velEq(x,u);
        dim_x = 6;
        dim_u = 8;
        out_fun = @(x,u) [x(4); x(5)] + [u(7); u(8)];
        dim_y = 2;
        sys = nonlinearSysDT('bicycle', fun, dt, dim_x, dim_u, out_fun, dim_y);

        if nargin > 1 && type == "rand" % ---
            c_R0 = [1.86; 3.46; 3.97; 5.39; 4.19; 6.85];
            c_W = [2.04; 8.78; 0.27; 6.70; 4.17; 5.59];
            c_V = [-0.02; 0.06];
            G_R0 = 0.01*[1.78 -1.37  0.79  0.60 -1.17 -1.4800
                1.77 -0.29  0.93 -0.54 -0.69  0.26
                -1.87  1.27 -0.49 -0.16  0.93 -2.02
                -1.05  0.07  1.80  0.61 -1.48  0.20
                -0.42  0.45  0.59 -1.04 -0.56  0.43
                1.40 -0.32 -0.64 -0.35 -0.03 -1.27];
            G_W = 0.01*diag([0.55 0.17 -0.19 0.58 -0.85 0.81]);
            G_W(1,2) = 0.2; G_W(2,5) = 1; G_W(6,1) = -0.5;
            G_V = 0.002*eye(dim_y);
        else % ---
            c_R0 = [1.2;0.5; 0; 0; 0; 0];
            c_W = zeros(dim_x,1);
            c_V = zeros(dim_y,1);
            G_R0 = eye(6);
            G_W = eye(dim_x);
            G_V = eye(dim_y);
        end
        R0 = zonotope([c_R0,G_R0]);
        W = zonotope([c_W,G_W]);
        V = zonotope([c_V,G_V]);
        U = cartProd(W, V);

    case "bicycleHO"
        % higher-order bicycle dynamics (see highorderBicycleDynamics.m)
        dt = 0.001;
        p_true = []; % no parameters defined
        fun = @(x,u) x + dt*highorderBicycleDynamics(x,u);
        dim_x = 18;
        dim_u = 4;
        out_fun = @(x,u) [x(5); x(6)] + [u(3); u(4)];
        dim_y = 2;
        sys = nonlinearSysDT('bicycleHO', fun, dt, dim_x, dim_u, out_fun, dim_y);

        if type ~= "standard"
            throw(CORAerror('CORA:specialError',"Only standard uncertainty sets defined."))
        end
        R0 = zonotope([[1.2;0.5; 0; 5; zeros(14,1)],0.01*eye(dim_x)]);
        W = zonotope([zeros(2,1),0.004*eye(2)]);
        V = zonotope([zeros(dim_y,1),0.002*eye(dim_y)]);
        U = cartProd(W, V);

    case "cstrDiscr"
        % stirred-tank reactor system [5]        
        dt = 0.015;
        p_true = []; % no parameters defined
        fun = @(x,u) cstrDiscr(x,u,dt);
        dim_x = 2;
        dim_u = 4;
        out_fun = @(x,u) [x(1); x(2)] + [u(3); u(4)];
        dim_y = 2;
        sys = nonlinearSysDT('cstrDiscr', fun, dt, dim_x, dim_u, out_fun, dim_y);
        
        if type ~= "standard"
            throw(CORAerror('CORA:specialError',"Only standard uncertainty sets defined."))
        end
        c_R0 = [-0.15;-45];
        c_W = zeros(dim_x,1);
        c_V = zeros(dim_y,1);
        G_R0 = diag([0.005;3]);
        G_W = diag([0.1;2]);
        G_V = 0.002*eye(dim_y);
        R0 = zonotope([c_R0,G_R0]);
        W = zonotope([c_W,G_W]);
        V = zonotope([c_V,G_V]);
        U = cartProd(W, V);

    case "tank"
        % stirred-tank reactor system with 6 dimensions [5]
        dt = 0.5;
        p_true = []; % no parameters defined
        fun = @(x,u) tank6EqDT(x,u,dt);
        dim_x = 6;
        dim_u = 4;
        out_fun = @(x,u) [x(1); x(2)] + [u(3); u(4)];
        dim_y = 2;
        sys = nonlinearSysDT('tank6', fun, dt, dim_x, dim_u, out_fun, dim_y);
        
        if type ~= "standard"
            throw(CORAerror('CORA:specialError',"Only standard uncertainty sets defined."))
        end
        c_R0 = [2; 4; 4; 2; 10; 4];
        c_W = zeros(2,1);
        c_V = zeros(dim_y,1);
        G_R0 = 0.2*eye(6);
        G_W = diag([0.1;2]);
        G_V = 0.002*eye(dim_y);
        R0 = zonotope([c_R0,G_R0]);
        W = zonotope([c_W,G_W]);
        V = zonotope([c_V,G_V]);
        U = cartProd(W, V);

    case "tank30"
        % stirred-tank reactor system with 30 dimensions [5]
        dt = 0.5;
        p_true = []; % no parameters defined
        fun = @(x,u) tank30EqDT_inflow15(x,u,dt);
        dim_x = 30;
        dim_u = 15;
        dim_y = 6;
        out_fun = @(x,u) x(1:dim_y) + u(1:dim_y);
        sys = nonlinearSysDT('tank30', fun, dt, dim_x, dim_u, out_fun, dim_y);
        
        if type ~= "standard"
            throw(CORAerror('CORA:specialError',"Only standard uncertainty sets defined."))
        end
        R0 = zonotope([12*rand(dim_x,1),0.2*eye(dim_x)]);
        W = zonotope([zeros(dim_y-dim_y,1),0.01*eye(dim_u-dim_y)]);
        V = zonotope([zeros(dim_y,1),0.002*eye(dim_y)]);
        U = cartProd(W, V);

    case "tank60"
        % stirred-tank reactor system with 60 dimensions [4]
        dt = 0.5;
        p_true = []; % no parameters defined
        fun = @(x,u) tank60EqDT_inflow30(x,u,dt);
        dim_x = 60;
        dim_u = 30;
        dim_y = 2;
        out_fun = @(x,u) x(1:dim_y) + u(1:dim_y);
        sys = nonlinearSysDT('tank60', fun, dt, dim_x, dim_u, out_fun, dim_y);
        
        if type ~= "standard"
            throw(CORAerror('CORA:specialError',"Only standard uncertainty sets defined."))
        end
        R0 = zonotope([12*rand(dim_x,1),0.2*eye(dim_x)]);
        W = zonotope([zeros(dim_y-dim_y,1),0.01*eye(dim_u-dim_y)]);
        V = zonotope([zeros(dim_y,1),0.002*eye(dim_y)]);
        U = cartProd(W, V);

end
end


% Auxiliary functions -----------------------------------------------------

function xnew = aux_dynLorenz2D(x,u,dt,p)
% dynamics of the Lorenz system

xdot = [(p(1)+u(1))*(x(2)-x(1)); ...
    (p(2)+u(2))*x(1)-x(2)-x(1)];
xnew = x + dt*xdot; 
end

function xnew = aux_dynLorenz(x,u,dt,p)
% dynamics of the Lorenz system

xdot = [(p(1)+u(1))*(x(2)-x(1)); ...
    (p(2)+u(2))*x(1)-x(2)-x(1)*x(3); ...
    x(1)*x(2)-(p(3)+u(3))*x(3)];
xnew = x + dt*xdot; 
end

% ------------------------------ END OF CODE ------------------------------
