function [sizeI_ddra, contain_pct] = ddra_infer_size_streaming(sys, R0, U, W, M_AB, C)
% Streaming DDRA: size proxy + MC containment, memory friendly.
% Unified reduction knob via red_order(C), with a legacy switch:
%   C.lowmem.legacy_ddra_streaming = true  -> match old behavior (no post-reduce)
%                                         = false -> also reduce post-image (recommended)

    % ---- unified reduction order (respects lowmem.cap inside helper) ----
    K = red_order(C);

    % ---- legacy toggle (default = old behavior preserved OFF? keep ON for safety) ----
    LM = getfielddef(C,'lowmem',struct());
    legacy_mode = getfielddef(LM,'legacy_ddra_streaming', true);   

    % ---- MC budget ----
    Nmc = getfielddef(LM,'ddra_mc', getfielddef(C,'mc_batch', 20));  % default 20

    % ---- main loop ----
    n_k_val = C.shared.n_k_val;
    X = R0;
    sizeI_ddra = 0;
    contain = 0; 
    total   = 0;

    for k = 1:n_k_val
        % pre-image reduction (unified K)
        X = reduce(X,'girard',K);

        % affine image with learned tube
        Xnext = M_AB * cartProd(X, U) + W;

        % optional post-image reduction:
        if ~legacy_mode
            % Only attempt reduce if object supports it
            try
                Xnext = reduce(Xnext,'girard',K);
            catch
                % some CORA versions may not reduce matZonotopes; safe to skip
            end
        end

        % size proxy on the (possibly reduced) set
        Iv = interval(Xnext);
        sizeI_ddra = sizeI_ddra + sum(abs(Iv.sup(:) - Iv.inf(:)));

        % on-the-fly MC containment
        if Nmc > 0
            for i = 1:Nmc
                x0 = randPoint(R0); u = randPoint(U); w = randPoint(W);
                x1 = sys.A*x0 + sys.B*u + w;
                if fast_contains_zonotope_point(Xnext, x1, 1e-3, 1e-6)
                    contain = contain + 1;
                end
                total = total + 1;
            end
        end

        X = Xnext;  % step forward
    end

    if total > 0
        contain_pct = 100 * contain / total;
    else
        contain_pct = NaN;
    end
end
