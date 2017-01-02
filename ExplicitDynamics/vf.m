function [f, dfdx, dfdu] = vf(x,u,params)

%     f = x + u;
%     dfdx = 1;
%     dfdu = 1;

%     % Shannon's LQR
%     A = params.A;
%     B = params.B;
%     C = params.C;
%     D = params.D;
%     
%     f = A*x' + B*u';
%     dfdx = A;
%     dfdu = [B(1) 0; 0 B(2)];
    
    % Dubin's Car
    v = u(1);
    w = u(2);
    theta = x(3);
    
    f = [v*cos(theta) ; v*sin(theta) ; w];
    dfdx = [0 0 -v*sin(theta);...
            0 0 v*cos(theta);...
            0 0 0];
    dfdu = [cos(theta) 0;...
            sin(theta) 0;...
            0 1];
    

end