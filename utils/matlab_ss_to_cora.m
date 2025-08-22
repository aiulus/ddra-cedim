function sys_cora = matlab_ss_to_cora(sys_ss, dt)
% If sys_ss is continuous, you must pass dt (sampling time).
% If sys_ss is already discrete, dt is ignored.

    if isct(sys_ss)
        assert(nargin>=2 && dt>0, 'Provide dt for c2d conversion.');
        sysd = c2d(sys_ss, dt);         % MATLAB c2d
        Ts = dt;
    else
        sysd = sys_ss;
        Ts = sysd.Ts;
        if Ts<=0
            error('System is discrete but Ts not set. Give a valid Ts.');
        end
    end

    % Build CORA linearSysDT
    A = sysd.A; B = sysd.B; C = sysd.C; D = sysd.D;
    sys_cora = linearSysDT(A,B,[],C,D,Ts);
end
