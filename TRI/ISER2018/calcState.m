function [xsol,exitflag,output] = calcState( p, xind0, params, fload)
%calcState - Solve for a the end effector state given input pressures
%   Detailed explanation goes here

% check for optional argument fload
if ~exist('fload','var')
     % fload input does not exist, so default it to 0
      fload = zeros(6,1);
end

num = params.num;   % number of actuators

% solve for x
% x0 = zeros(2,1);
fun = @(xind) xequation(xind, p, params);
options = optimoptions('lsqnonlin', 'Algorithm', 'Levenberg-Marquardt', 'FiniteDifferenceType', 'central', 'MaxFunctionEvaluations', 100, 'MaxIterations', 10);
% [xindsol,fval,exitflag,output] = fsolve(fun, xind0);
[xindsol,resnorm,residual,exitflag,output] = lsqnonlin(fun, xind0, [], [], options);

xsol = xind2xcoupled(xindsol)'; % get full state vector in terms of the independent states

end


function residual = xequation(xind, p, params)
%xequation - equation that must be solved to find x

x = xind2xcoupled(xind)'; % get full state vector in terms of the independent states

D = params.D;
q = x2q(x);
Jq = calcJq(q, params);

% calculate the elastomer force
felast = calcFelast( x, params );

residual = D*Jq' * p + felast;

end

% Plot the residual value over a rang of points to help visualize minimum
function showResidual(p, params)

x = linspace(-0.009, 0.001, 100);
y = linspace(-1.5, 1, 100);
[X,Y] = meshgrid(x,y);

for i = 1:length(x)
    for j = 1:length(y)
        residual = xequation([x(i), y(j)]', p, params);
        Z(i,j) = norm(residual, Inf);
    end
end

figure
surf(X,Y,Z);

end