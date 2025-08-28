function xnew = aux_dynLipschitz2D(x,u,dt)
    xdot = [ -0.6*x(1) + 0.4*tanh(x(2)) + 0.8*u(1);
              0.3*sin(x(1)) - 0.7*x(2)   + 0.5*u(2) ];
    xnew = x + dt*xdot;
end
