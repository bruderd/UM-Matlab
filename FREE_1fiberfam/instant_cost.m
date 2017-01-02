function [ic, dicdx, dicdu] = instant_cost(x, u, params)
    
    % Desired twist, 1 fiber-family family FREE
    phi_desired = params.phi_desired;
    ic = 10*(phi_desired - x(5))^2 + 0.01*u^2;
    dicdx = 10*[0 0 0 0 -2*(phi_desired - x(5))];
    dicdu = 2*u*0.01;
    
%     % Desired pressure
%     P_desired = 22;
%     ic =  (P_desired - x(1))^2;
%     dicdx = [-2*(P_desired - x(1)) 0 0 0 0];
%     dicdu = 0;

end