function [fc, dfcdx, dfcdu] = final_cost(x, u, params)
    
    % Dubin's Car, desired position (5, 6, pi)
    errorx = 5 - x(1);
    errory = 6 - x(2);
    errortheta = 3*pi/2 - x(3);
    fc = errorx^2 + errory^2 + errortheta^2;
    dfcdx = [-5*2 + 2*x(1), -6*2 + 2*x(2), -(3*pi/2)*2 + 2*x(3)];
    dfcdu = [0 0];

%     % No cost
%     fc = 0;
%     dfcdx = zeros(1,params.n);
%     dfcdu = zeros(1,params.m);  

end