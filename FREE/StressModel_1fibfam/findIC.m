% findIC.m
%   Finds the "globally" optimal solution of FBstress near a given pressure.
%   Solution can then be used as initial contition to iterative simulation,
%   away from zero.
%
%   NOTE: Must run setParams.m before using this script


function [IC, fval, exitFlag, output] = findIC(u, params, x0)

[P_rest, gama_rest, r_rest, L_rest, phi_rest, T_rest] = deal(params.x_rest(1), params.x_rest(2), params.x_rest(3), params.x_rest(4), params.x_rest(5), params.x_rest(6));

% Don't fix anything but the input pressure
dx = 1e-3;
LB = [u-dx, -pi/2, r_rest, L_rest*0.25, -Inf, 0];
UB = [u+dx, pi/2, r_rest*4, L_rest*2, Inf, Inf];

% objective function
Fnorm = @(x)FBstress_norm(x, u, params);

if nargin > 2 % if the user supplies an initial condition for solver, use it
    [IC,fval,exitFlag,output] = simulannealbnd(Fnorm,x0,LB,UB);
else
    [IC,fval,exitFlag,output] = simulannealbnd(Fnorm,params.x_rest,LB,UB);
end
    
    
end
