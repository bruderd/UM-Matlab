function [fc, dfcdx, dfcdu] = final_cost(x, u, params)

    % 1-fiber family FREE
    phi_desired = params.phi_desired;
    fc = 10000*(phi_desired - x(5))^2;
    dfcdx = 10000*[0 0 0 0 -2*(phi_desired - x(5))];
    dfcdu = 0;
    
%     % Desired pressure
%     P_desired = 22;
%     fc = P_desired - x(1);
%     dfcdx = [-1 0 0 0 0];
%     dfcdu = 0;
    
end