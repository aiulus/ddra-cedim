function Wd = normalizeWForGray(sys_, W_in)
% Conservative preimage: find Wd in disturbance space so that E*Wd ⊇ W_in.

    % Basic guards
    if isempty(W_in) || ~isa(W_in,'zonotope')
        q = getfielddef(sys_,'nrOfDisturbances',0);
        Wd = iff(q>0, zonotope(zeros(q,1)), []);  return;
    end
    E = getfielddef(sys_,'E',[]);
    if isempty(E)
        % fall back to B if you use E=B in your models
        E = getfielddef(sys_, 'B', []);
    end
    q = size(E,2);  if q==0, Wd = []; return; end

    % If W_in already lives in disturbance space, pass through
    if size(center(W_in),1) == q, Wd = W_in; return; end

    % Interval radius in state space
    Ix = interval(W_in);
    rad_x = 0.5*(supremum(Ix) - infimum(Ix));      % nx×1

    % Row-sum bound: choose scalar t so that |E w|_∞ ≥ rad_x elementwise
    rs = sum(abs(E),2);  rs(rs==0) = eps;
    t  = max(rad_x ./ rs);                           % scalar

    % Build diagonal q-dim box with radius t (center 0)
    Wd = zonotope(zeros(q,1), t*eye(q));
end
