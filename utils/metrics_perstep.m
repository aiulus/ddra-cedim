function [wid_k, cov_k, fv_blocks, wid_axes_k] = perstep_gray_metrics(sys_gray, VAL, Wpred, opts)
% Returns per-step width (scalar L1), per-step coverage (0..1), first-violation per block,
% and axis-wise widths (ny x n_k), computed only from VAL and the identified gray model.
% opts.zonotopeOrder (default 60)

    if ~isfield(opts,'zonotopeOrder'), opts.zonotopeOrder = 60; end
    ny   = size(safe_mat(sys_gray,'C',eye(sys_gray.nrOfDims)),1);
    n_k  = size(VAL.u{1},2);
    Bnum = numel(VAL.x0);

    cov_count = zeros(n_k,1);
    cov_tot   = zeros(n_k,1);
    wid_axes_acc = zeros(ny, n_k);
    fv_blocks = inf(Bnum,1);

    for b = 1:Bnum
        % build params per block
        U = VAL.u{b};
        params = struct('R0', VAL.R0 + VAL.x0{b}, 'u', U, ...
                        'tStart',0, 'tFinal', sys_gray.dt*(n_k-1));
        if sys_gray.nrOfDisturbances > 0 && ~isempty(Wpred)
            params.W = Wpred;
        end
        Rg = reach(sys_gray, params, struct('zonotopeOrder',opts.zonotopeOrder,'verbose',false));
        Yg = map_sets_to_output(Rg.timePoint.set, safe_mat(sys_gray,'C',eye(sys_gray.nrOfDims)), ...
                                safe_mat(sys_gray,'D',zeros(ny,size(U,1))), U);

        % measurements for this block (n_k x ny)
        Ymeas = VAL.y{b}; if size(Ymeas,1)~=n_k, Ymeas = Ymeas.'; end

        for k = 1:n_k
            Ik = interval(toZono(Yg{k}));
            wj = Ik.sup - Ik.inf;              % (ny x 1)
            wid_axes_acc(:,k) = wid_axes_acc(:,k) + wj;
            cov_tot(k) = cov_tot(k) + 1;

            yk = Ymeas(k,:).';
            inside = all(yk <= Ik.sup + 1e-9) && all(yk >= Ik.inf - 1e-9);
            cov_count(k) = cov_count(k) + inside;
            if ~inside && isinf(fv_blocks(b)), fv_blocks(b) = k; end
        end
    end

    wid_axes_k = wid_axes_acc ./ max(1,cov_tot');            % average across blocks
    wid_k      = sum(wid_axes_k, 1).';                       % scalar L1 width per step
    cov_k      = cov_count ./ max(1,cov_tot);                % per-step coverage [0,1]
end


function [wid_k, cov_k, fv_blocks, wid_axes_k] = perstep_ddra_metrics(sys_true_dt, VAL, M_AB, W_used, reduceOrder)
% DDRA propagation using M_AB (same as your plotter), then coverage vs VAL.y.
% reduceOrder default 60.

    if nargin<5 || isempty(reduceOrder), reduceOrder = 60; end
    ny   = size(safe_mat(sys_true_dt,'C',eye(sys_true_dt.nrOfDims)),1);
    n_k  = size(VAL.u{1},2);
    Bnum = numel(VAL.x0);

    cov_count = zeros(n_k,1);
    cov_tot   = zeros(n_k,1);
    wid_axes_acc = zeros(ny, n_k);
    fv_blocks = inf(Bnum,1);

    for b = 1:Bnum
        U = VAL.u{b};
        X = cell(1,n_k); X{1} = VAL.R0 + VAL.x0{b};
        for k = 1:n_k-1
            Zk = cartProd(reduce(toZono(X{k}),'girard',reduceOrder), zonotope(U(:,k)));
            X{k+1} = reduce(M_AB * Zk + W_used, 'girard', reduceOrder);
        end
        Y = map_sets_to_output(X, safe_mat(sys_true_dt,'C',eye(sys_true_dt.nrOfDims)), ...
                               safe_mat(sys_true_dt,'D',zeros(ny,size(U,1))), U);
        Ymeas = VAL.y{b}; if size(Ymeas,1)~=n_k, Ymeas = Ymeas.'; end

        for k = 1:n_k
            Ik = interval(toZono(Y{k}));
            wj = Ik.sup - Ik.inf;
            wid_axes_acc(:,k) = wid_axes_acc(:,k) + wj;
            cov_tot(k) = cov_tot(k) + 1;

            yk = Ymeas(k,:).';
            inside = all(yk <= Ik.sup + 1e-9) && all(yk >= Ik.inf - 1e-9);
            cov_count(k) = cov_count(k) + inside;
            if ~inside && isinf(fv_blocks(b)), fv_blocks(b) = k; end
        end
    end

    wid_axes_k = wid_axes_acc ./ max(1,cov_tot');
    wid_k      = sum(wid_axes_k, 1).';
    cov_k      = cov_count ./ max(1,cov_tot);
end
