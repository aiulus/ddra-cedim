function sysd = cora_to_matlab_ss(sys_cora)
% cora_to_matlab_ss  Adapt CORA DT systems to MATLAB-friendly forms.
%
% Returns:
%   - ss(A,B,C,D,Ts)             for linearSysDT
%   - ss(A,B,C,D,Ts)             for linearARX (lifted; stacked-input form)
%   - struct (nonlinear wrapper) for nonlinearSysDT / nonlinearARX
%   - struct (LTV wrapper)       for linearSysDT with cell A,B[,C,D]

    % Already MATLAB ss?
    if isa(sys_cora,'ss')
        sysd = sys_cora; 
        return;
    end

    % Sampling time (best-effort)
    Ts = get_prop(sys_cora,'dt',[]);

    % ---------- linearSysDT (LTI) ----------
    if isa(sys_cora,'linearSysDT') && ~iscell(sys_cora.A)
        A = sys_cora.A; B = sys_cora.B;
        C = fallback(sys_cora,'C',eye(size(A,1)));
        D = fallback(sys_cora,'D',zeros(size(C,1), size(B,2)));
        assert(~isempty(Ts),'cora_to_matlab_ss: missing dt for linearSysDT.');
        sysd = ss(A,B,C,D,Ts);
        return;
    end

    % ---------- linearSysDT (LTV via cell arrays) ----------
    if isa(sys_cora,'linearSysDT') && iscell(sys_cora.A) && iscell(sys_cora.B)
        Ak = sys_cora.A; Bk = sys_cora.B;
        Ck = fallback(sys_cora,'C',[]);
        Dk = fallback(sys_cora,'D',[]);
        assert(~isempty(Ts),'LTV system missing dt.');
        nx = size(Ak{1},1); nu = size(Bk{1},2);
        ny = ~isempty(Ck) * size(Ck{1},1) + isempty(Ck) * nx;
        sysd = struct();
        sysd.type = 'LTV';
        sysd.dt   = Ts;
        sysd.nx   = nx;
        sysd.nu   = nu;
        sysd.ny   = ny;
        sysd.Ak   = Ak; sysd.Bk = Bk; sysd.Ck = Ck; sysd.Dk = Dk;
        sysd.step = @(x,u,k) Ak{clamp_k(k,Ak)}*x + Bk{clamp_k(k,Bk)}*u;
        if ~isempty(Ck)
            sysd.output = @(x,u,k) Ck{clamp_k(k,Ck)}*x + Dk{clamp_k(k,Dk)}*u;
        else
            sysd.output = @(x,~,~) x;
        end
        sysd.meta = keep_some(sys_cora, {'name','dt'});
        return;
    end

    % ---------- linearARX -> lifted ss ----------
    if isa(sys_cora,'linearARX')
        Abar = sys_cora.A_bar;      % cell(1..p)
        Bbar = sys_cora.B_bar;      % cell(1..q)
        assert(~isempty(Ts),'linearARX: missing dt.');
        [Ass,Bss,Css,Dss] = lift_linearARX_to_ss(Abar,Bbar,'stacked-input');
        sysd = ss(Ass,Bss,Css,Dss,Ts);
        return;
    end

    % ---------- nonlinearSysDT (state-space) ----------
    if isa(sys_cora,'nonlinearSysDT')
        f   = sys_cora.mFile;
        out = fallback(sys_cora,'out_mFile',[]);
        nx  = sys_cora.nrOfDims;
        nu  = sys_cora.nrOfInputs;
        ny  = sys_cora.nrOfOutputs;
        assert(~isempty(Ts) && ~isempty(f) && ~isempty(nx) && ~isempty(nu), ...
               'nonlinearSysDT: missing dt/mFile/dims.');
        sysd = struct();
        sysd.type   = 'nonlinear';
        sysd.dt     = Ts;
        sysd.nx     = nx;
        sysd.nu     = nu;
        sysd.ny     = ny;
        sysd.step   = @(x,u,varargin) f(x,u);
        sysd.output = @(x,u) default_output(out,x,u,nx);
        sysd.meta   = keep_some(sys_cora, {'name','dt'});
        return;
    end

    % ---------- nonlinearARX ----------
    if isa(sys_cora,'nonlinearARX')
        f  = sys_cora.mFile;
        p  = sys_cora.n_p;
        ny = sys_cora.nrOfOutputs;
        nu = sys_cora.nrOfInputs;
        assert(~isempty(Ts) && ~isempty(f) && ~isempty(p), 'nonlinearARX: missing dt/mFile/n_p.');
        nx = ny * p;  % stacked output history
        sysd = struct();
        sysd.type   = 'nonlinearARX';
        sysd.dt     = Ts;
        sysd.nx     = nx;
        sysd.nu     = nu;
        sysd.ny     = ny;
        sysd.step   = @(x,u) arx_step_nonlinear(f, x, u, ny, nu, p);
        sysd.output = @(x,~) x(1:ny);
        sysd.meta   = struct('p',p,'name',get_prop(sys_cora,'name','nonlinearARX'));
        return;
    end

    error('cora_to_matlab_ss: unsupported CORA object of class "%s".', class(sys_cora));
end

% ---------------- helpers ----------------
function val = get_prop(obj, name, default)
    try
        val = obj.(name);
    catch
        val = default;
    end
end

function v = fallback(obj, f, d)
    try
        v = obj.(f);
        if isempty(v), v = d; end
    catch
        v = d;
    end
end

function i = clamp_k(k, Kcell)
    if nargin<1 || isempty(k), i = 1; return; end
    n = numel(Kcell); i = max(1, min(n, k));
end

function S = keep_some(obj, fields)
    S = struct();
    for i = 1:numel(fields)
        fn = fields{i};
        try, S.(fn) = obj.(fn); catch, end
    end
end

function y = default_output(out,x,u,nx)
    if ~isempty(out)
        try, y = out(x,u); return; catch, end
    end
    y = x;  % identity
end

function xnext = arx_step_nonlinear(f, xhist, u, ny, nu, p)
% xhist = [y(k); y(k-1); ...; y(k-p+1)]
% Many of your f(y,u) definitions accept stacked y,u with newest first.
    ywin = reshape(xhist, ny, []); 
    ynew = f(xhist, u);             % uses your chosen stacking convention
    xnext = [ynew; xhist(1:ny*(p-1))];
end

function [A,B,C,D] = lift_linearARX_to_ss(Abar,Bbar,mode)
% y(k) = sum_{i=1..p} Abar{i} y(k-i) + sum_{i=1..q} Bbar{i} u(k-i+1)
% state x(k) = [y(k); y(k-1); ...; y(k-p+1)]
    p  = numel(Abar);
    ny = size(Abar{1},1);
    q  = numel(Bbar);
    nu = size(Bbar{1},2);

    % State transition
    A = zeros(ny*p, ny*p);
    A(1:ny, 1:ny*p) = cell2mat(Abar(:).'); % [A1 ... Ap]
    if p>1
        A(ny+1:end, 1:ny*(p-1)) = eye(ny*(p-1));
    end

    C = [eye(ny), zeros(ny, ny*(p-1))];
    D = zeros(ny, nu);

    switch string(mode)
        case "stacked-input"
            % Input is stacked as u_stk(k) = [u(k); u(k-1); ...; u(k-q+1)]
            B = zeros(ny*p, nu*q);
            % top block gets [B1, B2, ..., Bq]
            for i = 1:q
                cols = (i-1)*nu + (1:nu);
                B(1:ny, cols) = Bbar{i};
            end
        otherwise % "current-only" (approximate, ignores B{i>1})
            B = zeros(ny*p, nu);
            if q>=1, B(1:ny,1:nu) = Bbar{1}; end
    end
end
