function [Xsets, sizeI] = ddra_reach_lipschitz(R0, U, W, D, C)
% Alg. 6: local LS model + Z_L + optional Z_eps (from Lipschitz).
% R0: initial set (zonotope), U: single input set (zonotope) or cell{1..N}, W: process noise
% D: struct(Xminus,Uminus,Xplus); C.nlip: .ZepsFlag (bool), .Lvec, .gamma (optional)

if ~iscell(U), U = repmat({U}, 1, getfielddef(C.shared,'n_k_val', getfielddef(C.shared,'n_k',1))); end
N  = numel(U);
nx = size(center(R0),1);

% LS linearization about sample means
xS = mean(D.Xminus,2);
uS = mean(D.Uminus,2);
Phi = [ones(1,size(D.Xminus,2)); D.Xminus - xS; D.Uminus - uS];   % [1; x-x*; u-u*]
Mhat = (D.Xplus) * pinv(Phi);

% residual enclosure Z_L (component-wise box)
res = D.Xplus - Mhat*Phi;      % residuals
rmax = max(res,[],2); rmin = min(res,[],2);
ZL = zonotope(0.5*(rmax+rmin), diag(0.5*(rmax-rmin)));

% Z_eps (optional): eps_i = 0.5 * L_i * gamma_i
Zeps = zonotope(zeros(nx,1));
if getfielddef(C.nlip,'ZepsFlag', true)
    Lvec  = getfielddef(C.nlip,'Lvec', []);
    gamma = getfielddef(C.nlip,'gamma',[]);
    if isempty(Lvec) || isempty(gamma)
        % on-demand estimate (lightweight)
        if isfield(C.nlip,'dyn_fun')
            [Lvec, gamma] = estimate_Lipschitz_from_data(C.nlip.dyn_fun, U{1}, R0, 1, 30, nx);
        else
            Lvec = zeros(nx,1); gamma = zeros(nx,1); % safe no-op
        end
    end
    eps = 0.5 * abs(Lvec(:)) .* abs(gamma(:));
    Zeps = zonotope(zeros(nx,1), diag(eps));
end

Xsets  = cell(N+1,1); Xsets{1} = R0;
sizeI  = 0;
for k=1:N
    Uk = U{k};
    Xaff = cartProd(Xsets{k} - xS, Uk - uS);
    Xnext = Mhat * aug1(Xaff) + W + ZL + Zeps;
    Xnext = reduce(Xnext, 'girard', getfielddef(C.lowmem,'zonotopeOrder_cap',100));
    Xsets{k+1} = Xnext;

    Iv = interval(Xnext); sizeI = sizeI + sum(abs(Iv.sup - Iv.inf));
end
end

function A = aug1(S)
c = center(S); G = generators(S);  % S in R^{n+m}
A = zonotope([1; c], blkdiag(zeros(1,size(G,2)), G));
end
