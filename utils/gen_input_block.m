function Ublk = gen_input_block(dim_u, n_k, U, pe)
    mode = getfielddef(pe,'mode','randn');
    switch mode
        case 'randn'
            Ublk = zeros(dim_u,n_k);
            for t=1:n_k, Ublk(:,t)=randPoint(U); end
        case 'multisine'
            r = max(1, round(getfielddef(pe,'order',2)));
            t = (0:n_k-1);
            Ublk = zeros(dim_u,n_k);
            G = generators(U); c = center(U);
            for i=1:dim_u
                sig = zeros(1,n_k);
                for j=1:r
                    f = (j)/(max(8,n_k)); ph = 2*pi*rand();
                    sig = sig + sin(2*pi*f*t + ph);
                end
                sig = sig / max(1e-9,max(abs(sig)));
                scale = 0.8*max(1e-9, max(abs(G(i,:))));
                Ublk(i,:) = c(i) + scale * sig;
            end
        otherwise
            error('Unknown PE mode: %s', mode);
    end
end