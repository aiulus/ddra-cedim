function sysdt = normalize_to_linearSysDT(sys_true, dt_target)
% normalize_to_linearSysDT: return CORA linearSysDT at dt_target
    if isa(sys_true,'linearSysDT')
        if abs(sys_true.dt - dt_target) < 1e-12
            sysdt = sys_true;
        else
            % re-discretize from continuous is not possible here; assume matrices are at dt_target already
            sysdt = linearSysDT(sys_true.A, sys_true.B, [], sys_true.C, sys_true.D, dt_target);
        end
        return
    end
    if isa(sys_true,'ss')
        sysd = d2d(sys_true, dt_target);
        sysdt = linearSysDT(sysd.A, sysd.B, [], sysd.C, sysd.D, dt_target);
        return
    end
    % struct with fields A,B,C,D,(dt) or numeric matrices
    if isstruct(sys_true) && all(isfield(sys_true, {'A','B','C','D'}))
        sysdt = linearSysDT(sys_true.A, sys_true.B, [], sys_true.C, sys_true.D, dt_target);
        return
    end
    % last resort: try to interpret as ss, then c2d
    sysc = ss(sys_true.A, sys_true.B, sys_true.C, sys_true.D);
    sysd = c2d(sysc, dt_target);
    sysdt = linearSysDT(sysd.A, sysd.B, [], sysd.C, sysd.D, dt_target);
end
