function plot_sets_overlay_2D(sys_gray, sys_true_dt, tc, params_gray, W_for_gray, R0, U, W_used, M_AB, Ccfg, dims)
% dims: e.g., [1 2]
    if nargin<11 || isempty(dims), dims = [1 2]; end
    ny = sys_true_dt.nrOfOutputs;
    assert(all(dims>=1 & dims<=ny), 'dims exceed output dimension');

    % --- RCSI reach on this tc
    nk = size(tc.y,1);
    p = struct();
    p.R0 = params_gray.R0 + tc.initialState;
    p.u  = tc.u';
    du = sys_gray.nrOfInputs - size(p.u,1); if du>0, p.u=[p.u; zeros(du, size(p.u,2), size(p.u,3))]; end
    p.tFinal = sys_gray.dt*(nk-1);
    if sys_gray.nrOfDisturbances>0, p.W = W_for_gray; end

    opts = rmfield_if(params_gray, 'options'); % defensive
    opts = rmfield_if(opts, 'cs');
    Rg = reach(sys_gray, p, getfielddef(opts,'',struct())); Rg = Rg.timePoint.set; % state sets

    % Map to outputs (RCSI)
    Yg = cell(nk,1);
    for k=1:nk
        Xk = Rg{k}; if ~isa(Xk,'contSet'), Xk = zonotope(Xk); end
        try, Yg{k} = linearMap(Xk, sys_gray.C); catch, Yg{k} = Xk; end
    end

    % --- DDRA inference on this tc
    VALsingle = struct(); VALsingle.x0 = {tc.initialState};
    VALsingle.u  = {squeeze(tc.u(:,:,1))};
    VALsingle.y  = {squeeze(tc.y(:,:,1))};
    [Xsets_ddra, ~] = ddra_infer(sys_true_dt, R0, U, W_used, M_AB, Ccfg, VALsingle);
    Yd = cellfun(@(X) map_to_output(X, sys_true_dt.C), Xsets_ddra, 'uni', 0);

    % --- figure
    figure('Color','w'); hold on; box on; grid on;

    % Nominal outputs (black line + crosses)
    yn = squeeze(tc.y(:,dims,1));
    plot(yn(:,1), yn(:,2), 'k-', 'LineWidth',1.0, 'DisplayName','nominal y');
    plot(yn(:,1), yn(:,2), 'kx', 'MarkerSize',6, 'HandleVisibility','off');

    % Initial state at output (black dot)
    y0 = sys_true_dt.C * tc.initialState;
    plot(y0(dims(1)), y0(dims(2)), 'ko', 'MarkerFaceColor','k', 'MarkerSize',6, 'DisplayName','initial');

    % Plot sets
    colG = [0.10 0.45 0.85];  % blue
    colD = [0.90 0.40 0.00];  % orange
    for k=1:nk
        if isa(Yg{k},'contSet')
            plot(Yg{k}, dims, 'FaceColor', colG, 'FaceAlpha', 0.12, 'EdgeColor', colG, 'LineWidth',0.5);
        end
        if isa(Yd{k},'contSet')
            plot(Yd{k}, dims, 'FaceColor', colD, 'FaceAlpha', 0.12, 'EdgeColor', colD, 'LineWidth',0.5);
        end
    end
    legend('Location','best');
    xlabel(sprintf('y_{%d}',dims(1))); ylabel(sprintf('y_{%d}',dims(2)));
    title(sprintf('Reachable sets overlay (dims %d-%d)',dims(1),dims(2)));
end

function s = rmfield_if(s, f)
    if isstruct(s) && isfield(s,f), s = rmfield(s,f); end
end
