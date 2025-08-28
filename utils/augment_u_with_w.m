function sys2 = augment_u_with_w(sys)
% Fallback: linearSysDT only.
% Converts x+ = A x + B u  into  x+ = A x + [B  I] [u; w]
% --> fold W into U via cartProd(U,W).

if ~isa(sys,'linearSysDT')
    error('augment_u_with_w: only implemented for linearSysDT in this fallback.');
end
A  = sys.A;  B  = sys.B;  C = sys.C;  D = sys.D;  dt = sys.dt;
nx = size(A,1);
Bw = eye(nx);                      % additive process noise enters state directly
Dw = zeros(size(C,1), nx);         % no direct feedthrough of w to y
sys2 = linearSysDT(A, [B, Bw], [], C, [D, Dw], dt);
end
