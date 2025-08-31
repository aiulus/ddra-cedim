function ok = check_PE_order(TS, L)
%CHECK_PE_ORDER  Verify mL-order PE on the training inputs using a block-Hankel rank test.
% Accepts:
%   • TS = cell array of CORA testCase (preferred)
%   • TS = struct with fields U_blocks (m×n_k×M) or U (m×T or n_k×m)
%   • TS = numeric U directly (m×T or n_k×m)
%
% Returns:
%   ok = (rank(Hankel(U,L)) == m*L)

    % ---- 1) Collect U across all training blocks as m×T ----
    U = ts_concat_U(TS);            % m × T
    [m,T] = size(U);

    % ---- 2) Quick feasibility guard ----
    if L < 1 || ~isfinite(L), ok = false; return; end
    L = floor(L);
    if T < 2*L
        ok = false; return;
    end

    % ---- 3) Build block-Hankel and test rank ----
    cols = T - L + 1;
    H = zeros(m*L, cols);
    for i = 1:L
        H((i-1)*m+1:i*m, :) = U(:, i:i+cols-1);
    end
    ok = (rank(H) == m*L);
end

function Ucat = ts_concat_U(TS)
%TS_CONCAT_U  Normalize and concatenate inputs to m×T.
% Handles:
%   • cell array of testCase: each tc.u can be [n_k×m], [n_k×m×s], or [m×n_k]
%   • struct with U_blocks (m×n_k×M) or U (m×T or n_k×m)
%   • numeric U (m×T or n_k×m)

    % ----- cell array of testCase -----
    if iscell(TS) && ~isempty(TS) && isa(TS{1}, 'testCase')
        Uall = {};
        for i = 1:numel(TS)
            Ui = TS{i}.u;  % usually n_k×m or n_k×m×1
            Ui = squeeze(Ui);      % drop singleton 3rd dim if present
            if size(Ui,1) < size(Ui,2)   % (m×n_k) → transpose
                Ui = Ui.';
            end
            % Now Ui is n_k×m; transpose to m×n_k for concatenation
            Uall{end+1} = Ui.'; 
        end
        Ucat = cat(2, Uall{:});   % m × T
        return;
    end

    % ----- struct with U_blocks / U -----
    if isstruct(TS)
        if isfield(TS, 'U_blocks') && ~isempty(TS.U_blocks)
            % U_blocks: (m × n_k × M) or (n_u × n_k × M)
            Ub = TS.U_blocks;
            if size(Ub,1) > size(Ub,2)    % likely (m×n_k×M) → good
                Ucat = reshape(permute(Ub, [1 3 2]), size(Ub,1), []); % m × (M*n_k)
            else                           % (n_k×m×M) → swap
                Ub = permute(Ub, [2 1 3]); % m×n_k×M
                Ucat = reshape(permute(Ub, [1 3 2]), size(Ub,1), []);
            end
            return;
        elseif isfield(TS,'U') && ~isempty(TS.U)
            U = TS.U;
            if size(U,1) < size(U,2)      % (n_k×m) → m×n_k
                U = U.';
            end
            Ucat = U;                      % m × T
            return;
        end
    end

    % ----- numeric U -----
    if isnumeric(TS) && ~isempty(TS)
        U = TS;
        if size(U,1) < size(U,2)          % (n_k×m) → (m×n_k)
            U = U.';
        end
        Ucat = U;
        return;
    end

    error('check_PE_order: TS missing recognizable inputs (testCase cell, U_blocks/U, or numeric U).');
end
