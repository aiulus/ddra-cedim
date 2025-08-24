function [Xminus,Uminus,Xplus,W,Zinfo] = ddra_generate_data_mult(sys, R0, U, C, pe)
    % like ddra_generate_data, but draws A_k from A_mz each step
    dim_x = size(sys.A,1);
    W = zonotope(zeros(dim_x,1), C.ddra.alpha_w*ones(dim_x, C.ddra.eta_w));
    A_mz = build_A_mz(sys.A, getfielddef(C,'noise',struct()).mult);

    n_m=C.shared.n_m; n_s=C.shared.n_s; n_k=C.shared.n_k;
    Ublocks = cell(n_m*n_s,1);
    for b=1:numel(Ublocks), Ublocks{b} = gen_input_block(size(sys.B,2), n_k, U, pe); end
    Xminus=[]; Uminus=[]; Xplus=[];

    for b=1:numel(Ublocks)
        x = randPoint(R0);
        for t=1:n_k
            Xminus(:,end+1) = x;
            u = Ublocks{b}(:,t); Uminus(:,end+1) = u;
            % sample one A_k from A_mz
            Ak = A_mz.center; 
            if ~isempty(A_mz.generator)
                % sample uniform \xi in [-1,1] for each generator
                for g=1:numel(A_mz.generator)
                    Ak = Ak + (2*rand-1)*A_mz.generator{g};
                end
            end
            w = randPoint(W);
            x = Ak*x + sys.B*u + w;
            Xplus(:,end+1) = x;
        end
    end
    Z = [Xminus;Uminus];
    Zinfo.rankZ = rank(Z); Zinfo.condZ = cond(Z*Z');
end
