function [Xminus, Uminus, Xplus, W, Zinfo] = ddra_generate_data(sys, R0, U, C, pe)
    % Simulate blocks to build (X^-, U^-, X^+) and define W
    n_m=C.shared.n_m; n_s=C.shared.n_s; n_k=C.shared.n_k;
    dim_x=size(sys.A,1);
    W = zonotope(zeros(dim_x,1), C.ddra.alpha_w*ones(dim_x, C.ddra.eta_w));
    Ublocks = cell(n_m*n_s,1);
    for b=1:numel(Ublocks), Ublocks{b} = gen_input_block(size(sys.B,2), n_k, U, pe); end
    Xminus=[]; Uminus=[]; Xplus=[];
    for b=1:numel(Ublocks)
        x = randPoint(R0);
        for t=1:n_k
            Xminus(:,end+1) = x; 
            u = Ublocks{b}(:,t); Uminus(:,end+1) = u; 
            w = randPoint(W); x = sys.A*x + sys.B*u + w; Xplus(:,end+1) = x; 
        end
    end
    Z = [Xminus; Uminus];
    Zinfo.rankZ = rank(Z); Zinfo.condZ = cond(Z*Z');
end