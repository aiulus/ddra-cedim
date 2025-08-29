function [Y, U] = case_to_YU(tc)
% Y: (ny × T), U: (m × T). Works for 2-D or 3-D CORA arrays.

    Y = [];
    U = [];

    % ----- Y -----
    if ~isempty(tc.y)
        y = tc.y;
        % CORA: (T × ny × s). Use s=1 (first sample) by default.
        if ndims(y) == 3
            y2 = y(:,:,1);            % (T × ny)
        else
            y2 = y;                   % maybe already (T × ny)
        end
        if size(y2,1) < size(y2,2)    % if (ny × T), transpose to (T × ny)
            y2 = y2.';
        end
        Y = permute(y2, [2 1]);       % -> (ny × T)
    end

    % ----- U -----
    if ~isempty(tc.u)
        u = tc.u;
        if ndims(u) == 3
            u2 = u(:,:,1);            % (T × m)
        else
            u2 = u;
        end
        % ensure rows = time
        if size(u2,1) < size(u2,2)    % if (m × T), transpose to (T × m)
            u2 = u2.';
        end
        U = permute(u2, [2 1]);       % -> (m × T)
    end
end
