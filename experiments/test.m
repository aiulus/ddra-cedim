% --- quick synthetic dims ---
nx=3; nu=2; ny=2; dt=0.1; nk=5;
A=eye(nx)+0.1*randn(nx); B=randn(nx,nu); C=randn(ny,nx); D=zeros(ny,nu);
sysdt = linearSysDT(A,B,[],C,D,dt);

% sets
R0 = zonotope(zeros(nx,1), 0.2*eye(nx));
U  = zonotope(zeros(nu,1), 0.2*eye(nu));
W  = zonotope(zeros(nx,1), 0.01*eye(nx));

% PE input & VAL
[U_nom, ~] = genPEInput("sinwave", 3, nu, nk, dt, U, struct('strength',1.5));
VAL = struct('x0',{{zeros(nx,1)}}, 'u',{{U_nom'}}, 'y',{{}}, 'R0',R0, 'U',U);

% DDRA "identity" M_AB just to exercise code
M_AB = matZonotope([A B], 0.0*ones(nx, size([A B],2), 0));

% reach (true vs gray both = sysdt here)
plot_reach_all_onepanel(sysdt, sysdt, VAL, 'MAB', M_AB, 'W', W, 'Dims',[1 2], 'ShowSamples', false, 'Save','./plots/reach_smoke');
save_plot(gcf, './plots', 'reach_smoke_copy', 'Formats', {'png','pdf'});
