%% Solve for xdot
    % solves the implicit differential equations to solve for xdot at a
    % given point x. This is used to check if there are multiple solutions
    % at a given point.
function [solution] = solve4xdot(x, u, params)

     fun = @(dxdt)vf(x, u, dxdt, params);
     
     xdot0 = [0, 0, 0, 0];
     
     options=optimset('TolFun', 1e-12);
     
     solution = lsqnonlin(fun,xdot0);
%     solution = lsqnonlin(fun,xdot0,[-Inf,-Inf,0,-Inf],[Inf,Inf,Inf,Inf],options);