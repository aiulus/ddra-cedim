function TS = testSuite_fromDDRA(sys_cora, ~, DATASET, n_k_override)
% Build CORA testCase objects from DATASET.
% Accepts BOTH new (x0_list/U_blocks/Y_blocks) and old (x0_blocks/u_blocks) layouts.

    % --- detect layout & define accessors ---
    if isfield(DATASET,'x0_list') && isfield(DATASET,'U_blocks')
        M       = size(DATASET.x0_list, 2);
        get_x0  = @(b) DATASET.x0_list(:, b);              % (n_x×1)
        get_U   = @(b) DATASET.U_blocks(:, :, b);          % (n_u×n_k)
        have_Y  = isfield(DATASET,'Y_blocks');
        if have_Y, get_Y = @(b) DATASET.Y_blocks(:, :, b); % (n_y×n_k)
        end
        n_k_ds  = size(DATASET.U_blocks, 2);
    elseif isfield(DATASET,'x0_blocks') && isfield(DATASET,'u_blocks')
        M       = numel(DATASET.x0_blocks);
        get_x0  = @(b) DATASET.x0_blocks{b};               % (n_x×1)
        get_U   = @(b) DATASET.u_blocks{b};                % (n_u×n_k)
        have_Y  = false;
        n_k_ds  = size(get_U(1), 2);
    else
        error('testSuite_fromDDRA: DATASET missing x0/U fields.');
    end

    % --- pick horizon ---
    if nargin >= 4 && ~isempty(n_k_override)
        n_k = n_k_override;
    else
        n_k = n_k_ds;
    end

    % --- construct testCase objects ---
    ny  = sys_cora.nrOfOutputs;
    TS  = cell(1, M);
    for b = 1:M
        x0   = get_x0(b);         % (n_x×1)
        Ublk = get_U(b);          % (n_u×n_k)

        % match requested horizon by crop/pad (rare)
        if size(Ublk,2) > n_k, Ublk = Ublk(:,1:n_k); end
        if size(Ublk,2) < n_k, Ublk = [Ublk, zeros(size(Ublk,1), n_k-size(Ublk,2))]; end

        % CORA testCase expects time-major u: (a×p×s) = (n_k × n_u × 1)
        u_tc = reshape(Ublk.', n_k, size(Ublk,1), 1);

        % Measured outputs y: (a×q×s) = (n_k × n_y × 1)
        if have_Y
            %Yblk = get_Y(b);            % (n_y×n_k)
            %y_tc = reshape(Yblk.', n_k, ny, 1);
            y_tc = permute(DATASET.Y_blocks(:,:,b), [2 1 3]);
        else
            % Fallback: if DATASET has no Y, hand in NaNs (RCSI can still reach/compare)
            %y_tc = nan(n_k, ny, 1);
            y_tc = nan(n_k, ny, 1);
        end

        % Create the actual CORA testCase object
        TS{b} = testCase(y_tc, u_tc, x0, sys_cora.dt, class(sys_cora));
    end
end