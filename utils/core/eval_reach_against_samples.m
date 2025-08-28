function [contain_pct, sizeI] = eval_reach_against_samples(R_data, testSuite_val, sys)
% Validate containment of samples in reach sets.
% - Prefer state-space validation when dimensions match.
% - Otherwise project state set via linear output matrix C (if available).
% - Fallback: take first n_x rows of y with a warning.

    if nargin < 3, sys = []; end

    % infer set dimension from the first reach set
    assert(~isempty(R_data) && isprop(R_data,'timePoint') && ~isempty(R_data.timePoint.set), ...
        'R_data is empty.');
    S1 = R_data.timePoint.set{1};
    box1 = interval(S1);
    nx = numel(box1.inf);

    % Attempt linear output projection
    use_linear_output = (~isempty(sys) && isprop(sys,'C') && ~isempty(sys.C) && size(sys.C,2) == nx);
    C = []; if use_linear_output, C = sys.C; end

    num_in = 0; num_out = 0; size_acc = 0; steps_tot = 0;

    for m = 1:length(testSuite_val)
        has_x = isfield(testSuite_val{m}, 'x') && ~isempty(testSuite_val{m}.x);
        has_y = isfield(testSuite_val{m}, 'y') && ~isempty(testSuite_val{m}.y);

        if ~has_x && ~has_y
            warning('eval_reach_against_samples: no x or y in testSuite{%d}; skipping.', m);
            continue;
        end

        % pull sets for this suite
        sets = R_data.timePoint.set;
        % figure out horizon alignment from whichever data is available
        if has_x
            X = testSuite_val{m}.x;   % [nx? x (k+1) x ns]
            T = min(numel(sets), size(X,2)-1);
        else
            Y = testSuite_val{m}.y;   % [ny x (k+1) x ns]
            T = min(numel(sets), size(Y,2)-1);
        end

        for k = 1:T
            S = sets{k};
            % size proxy always in set space
            b = interval(S); lbS = b.inf; ubS = b.sup;
            size_acc = size_acc + sum(abs(ubS - lbS));
            steps_tot = steps_tot + 1;

            if has_x && size(testSuite_val{m}.x,1) == nx
                % --- State-space validation (preferred) ---
                X = testSuite_val{m}.x;
                xk_all = squeeze(X(:,k+1,:));  % nx x ns
                inside = all(xk_all >= lbS & xk_all <= ubS, 1);

            elseif use_linear_output && has_y
                % --- Output-space validation via linear map y = Cx ---
                CS = C * S;
                bY = interval(CS); lbY = bY.inf; ubY = bY.sup;   % ny x 1
                Y = testSuite_val{m}.y;
                yk_all = squeeze(Y(:,k+1,:));                   % ny x ns
                inside = all(yk_all >= lbY & yk_all <= ubY, 1);

            elseif has_y && size(testSuite_val{m}.y,1) >= nx
                % --- Fallback: take first nx rows of y ---
                if k == 1
                    warning('eval_reach_against_samples: falling back to first %d rows of y for validation (suite %d).', nx, m);
                end
                Y = testSuite_val{m}.y;
                yk_all = squeeze(Y(1:nx, k+1, :));              % nx x ns
                inside = all(yk_all >= lbS & yk_all <= ubS, 1);

            else
                warning('eval_reach_against_samples: dimension mismatch at suite %d, step %d; skipping.', m, k);
                continue;
            end

            num_in  = num_in  + sum(inside);
            num_out = num_out + sum(~inside);
        end
    end

    contain_pct = 100 * num_in / max(1, (num_in+num_out));
    sizeI = size_acc / max(1, steps_tot);
end
