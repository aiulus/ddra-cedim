function sys2 = augment_u_with_v(sys)
% Fallback: linearSysDT only.
% Converts y = C x + D u  into  y = C x + [D  I] [u; v]

if ~isa(sys,'linearSysDT')
    error('augment_u_with_v: only implemented for linearSysDT in this fallback.');
end
A  = sys.A;  B  = sys.B;  C = sys.C;  D = sys.D;  dt = sys.dt;
ny = size(C,1); nx = size(A,1);
Bv = zeros(nx, ny);                % v does NOT affect the state, only outputs
Dv = eye(ny);
sys2 = linearSysDT(A, [B, Bv], [], C, [D, Dv], dt);
end
