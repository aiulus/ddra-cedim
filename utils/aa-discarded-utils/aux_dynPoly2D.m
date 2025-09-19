function xnew = aux_dynPoly2D(x,u,dt)
    xdot = [ -0.5*x(1) + 0.1*x(1)^3 + 0.2*x(1)*x(2) + 0.7*u(1);
              0.15*x(1)^2 - 0.6*x(2) + 0.1*x(2)^3 + 0.5*u(2) ];
    xnew = x + dt*xdot;
end
