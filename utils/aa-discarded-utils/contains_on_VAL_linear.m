function contain_pct = contains_on_VAL_linear(sys, W, Xsets, VAL, tol)
% CONTAINS_ON_VAL_LINEAR
%   Compute % of simulated points (using VAL.{x0,u,w}) contained in DDRA sets Xsets.
%   VAL.u{b} is expected to be (n_k × n_u). 

    if nargin < 5 || isempty(tol), tol = 1e-6; end

    B = numel(VAL.x0);
    num_in = 0; num_all = 0;

    % Precompute zero/center disturbance for the "no recorded w" case
    if isa(W,'zonotope')
        w_def = center(W);
    else
        w_def = zeros(size(sys.A,1),1);
    end
    if isempty(w_def), w_def = zeros(size(sys.A,1),1); end

    for b = 1:B
        x   = VAL.x0{b};
        Ubk = VAL.u{b};                  % (n_k × n_u)
        n_k = size(Ubk,1);

        % Optional exact disturbance record
        hasW = isfield(VAL,'w') && numel(VAL.w) >= b && ~isempty(VAL.w{b});

        for k = 1:n_k
            % NOTE: VAL.u{b} is time-major -> take k-th ROW and transpose
            u_k = Ubk(k,:).';            % (n_u × 1)

            if hasW
                w_k = VAL.w{b}(:,k);     % (n_x × 1)
            else
                w_k = w_def;             % center (or zeros)
            end

            % One-step simulate
            x = sys.A*x + sys.B*u_k + w_k;

            % Check containment in the DDRA set for step k
            Zk = Xsets{k};               % expected cell {1..n_k}
            if ~isa(Zk,'contSet'), Zk = zonotope(Zk); end
            if contains_interval(x, Zk, tol)
                num_in = num_in + 1;
            end
            num_all = num_all + 1;
        end
    end

    contain_pct = (num_all > 0) * (100 * num_in / num_all);
    if num_all == 0, contain_pct = NaN; end
end
