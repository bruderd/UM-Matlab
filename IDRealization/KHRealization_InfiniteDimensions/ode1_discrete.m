 function yout = ode1_discrete(F,t0,h,tfinal,y0,input_t,input_u,param)
   % ODE1  A simple ODE solver.
   %   yout = ODE1(F,t0,h,tfinal,y0) uses Euler's
   %   method with fixed step size h on the interval
   %      t0 <= t <= tfinal
   %   to solve
   %      dy/dt = F(t,y)
   %   with y(t0) = y0.

   %   Copyright 2014 - 2015 The MathWorks, Inc.
   
      y = y0;
      yout = y;
      for t = t0 : h : tfinal-h
         y = F(t,y,input_t,input_u,param);
         yout = [yout y];
      end