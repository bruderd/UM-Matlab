function [ic, dicdx, dicdu] = instant_cost(x, u, params)
    
%     % Dubin's Car, desired position (5, 6)
%     errorx = 5 - x(1);
%     errory = 6 - x(2);
%     erroru = 1 - u(1);
%     ic = errorx^2 + errory^2 + 100*erroru^2;
%     dicdx = [-5*2 + 2*x(1), -6*2 + 2*x(2), 0];
%     dicdu = [100*(-1*2 + 2*u(1)), 0];
     
%     % No cost
%     ic = 0;
%     dicdx = zeros(1,params.n);
%     dicdu = zeros(1,params.m);    
    
    % Keep constant input
    ic = 100*(u'*u - params.u0'*params.u0)^2;
    dicdx = zeros(1,params.n);
    dicdu = 100*2*(u'*u - params.u0'*params.u0)*2*u';

end