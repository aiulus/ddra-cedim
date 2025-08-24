function [Xsets, sizeI] = ddra_reach_poly(R0, U, W, D, C)
% Quadratic lifting (like example_poly). U can be a single set or cell{1..N}.

if ~iscell(U), U = repmat({U}, 1, getfielddef(C.shared,'n_k_val', getfielddef(C.shared,'n_k',1))); end
N  = numel(U);
nx = size(center(R0),1);
nu = size(center(U{1}),1);

% Build design matrix from data
X0 = D.Xminus;       % (nx x T)
X1 = D.Xplus;        % (nx x T)
U0 = D.Uminus;       % (nu x T)
Phi = quad_design(X0, U0);   % m x T

% Matrix-zonotope for W across all samples (new CORA API)
Wmz = make_W_matrix_zono(W, size(X1,2));

% Mp as matrix-zonotope: (X1 + (-1)*Wmz) * pinv(Phi)
% Most CORA overloads accept: matZonotope * numeric → matZonotope
ABmz = (X1 + (-1)*Wmz) * pinv(Phi);

% propagate
Xsets = cell(N+1,1); Xsets{1} = R0; sizeI = 0;
for k=1:N
    % interval lift of current (Xk,Uk)
    IX = interval(Xsets{k}); IU = interval(U{k});
    H  = quad_design_interval(IX, IU);   % interval vector → turn to zono
    c  = 0.5*(H.sup + H.inf);
    g  = 0.5*(H.sup - H.inf);
    Hset = zonotope(c, diag(g));
    Xnext = ABmz * Hset + W;             % set-map through matZonotope
    Xnext = reduce(Xnext,'girard', getfielddef(C.lowmem,'zonotopeOrder_cap',100));
    Xsets{k+1} = Xnext;

    Iv = interval(Xnext); sizeI = sizeI + sum(abs(Iv.sup - Iv.inf));
end
end

function Phi = quad_design(X,U)
% Build [1; X; X.^2; x1*x2; U; U.^2; u1*u2; X.*U; x1*u2; x2*u1] for general nx,nu
[nx,T] = size(X); [nu,Tu] = size(U); assert(T==Tu);
Phi = [ ones(1,T); X; X.^2 ];
if nx>=2, Phi = [Phi; X(1,:).*X(2,:)]; end
Phi = [Phi; U; U.^2];
if nu>=2, Phi = [Phi; U(1,:).*U(2,:)]; end
Phi = [Phi; X.*U];
if nx>=2 && nu>=2
    Phi = [Phi; X(1,:).*U(2,:); X(2,:).*U(1,:)]; 
end
end

function H = quad_design_interval(IX,IU)
% same basis on intervals; return struct with .inf, .sup
combine = @(a,b) struct('inf', [a.inf;b.inf], 'sup', [a.sup;b.sup]);
one.inf = 1; one.sup = 1;
X2.inf  = IX.inf.^2; X2.sup  = IX.sup.^2;
U2.inf  = IU.inf.^2; U2.sup  = IU.sup.^2;

% cross terms via interval arithmetic
mul = @(A,B) struct('inf', min(A.inf.*B.inf, A.inf.*B.sup, A.sup.*B.inf, A.sup.*B.sup), ...
                    'sup', max(A.inf.*B.inf, A.inf.*B.sup, A.sup.*B.inf, A.sup.*B.sup));
XU  = mul(IX, IU);

H = combine(struct('inf',one.inf,'sup',one.sup), IX);
H = combine(H, X2);
if size(IX.inf,1)>=2, H = combine(H, mul(IX(1,:), IX(2,:))); end
H = combine(H, IU);
H = combine(H, U2);
if size(IU.inf,1)>=2, H = combine(H, mul(IU(1,:), IU(2,:))); end
H = combine(H, XU);
if size(IX.inf,1)>=2 && size(IU.inf,1)>=2
    H = combine(H, mul(IX(1,:), IU(2,:)));
    H = combine(H, mul(IX(2,:), IU(1,:)));
end
end
