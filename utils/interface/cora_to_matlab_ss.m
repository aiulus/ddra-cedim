function sysd_ss = cora_to_matlab_ss(sys_cora)
% Only supports linearSysDT. Returns a discrete-time ss (with Ts).
    assert(isa(sys_cora,'linearSysDT'), ...
        'cora_to_matlab_ss: only linearSysDT supported');
    A = sys_cora.A; B = sys_cora.B; C = sys_cora.C; D = sys_cora.D;
    Ts = sys_cora.dt;
    sysd_ss = ss(A,B,C,D,Ts);
end
