function [f, dfdx, dfdu, dfdxdot] = vf(x, u, xdot, params)

    n = params.n;
    m = params.m;

%% Simple Practice System
%     f = x + u;
%     dfdx = 1;
%     dfdu = 1;
%     dfdxdot = 0;

%% Shannon's LQR
%     A = params.A;
%     B = params.B;
%     C = params.C;
%     D = params.D;
%     
%     f = A*x' + B*u';
%     dfdx = A;
%     dfdu = [B(1) 0; 0 B(2)];
    
%% Dubin's Car
%     v = u(1);
%     w = u(2);
%     theta = x(3);
%     
%     f = [v*cos(theta) ; v*sin(theta) ; w];
%     dfdx = [0 0 -v*sin(theta);...
%             0 0 v*cos(theta);...
%             0 0 0];
%     dfdu = [cos(theta) 0;...
%             sin(theta) 0;...
%             0 1];
%     dfdxdot = zeros(n,n);
    
%% simple non-linear implicit dynamics: xdot^2 = -x*u
%     f = -xdot - x*u;
%     dfdx = -u;
%     dfdu = -x;
%     df_ddxdt = -eye(n);
    
%% Weird made up test case for single input system
%     f = zeros(n,1);
%     dfdx = zeros(n,n);
%     dfdu = zeros(n,1);
%     dfdxdot = zeros(n,n);
%     
%     for i = 1:n
%         f(i,1) = -xdot(i)^2 + xdot(i) - x(i) + i*u;
%         
%         dfdx(i,i) = -1;
%         
%         dfdu(i,1) = i;
%         
%         dfdxdot(i,i) = -2*xdot(i) + 1;  
%     end
    
%% Crazy system in R4 with single input
%     f = [-(xdot(1) - u)^2 + cos(x(2)*x(3)) - 10*xdot(4)/x(1)];
%          (x(1)*xdot(2))^2 - sin(xdot(3)*x(4));...
%          xdot(1)*xdot(2)*xdot(3)*xdot(4) - x(4)*x(2) + x(3);...
%          -tan(xdot(3)^3 * x(1)) + (x(4)^2 + 2*x(2) - x(1))/(3*xdot(3))];
%     
%     dfdx = [10*xdot(4)/x(1)^2, -sin(x(2)*x(3)), -sin(x(2)*x(3)), 0];
%             2*x(1)*xdot(2)^2, 0, 0, -cos(xdot(3)*x(4));...
%             0, -x(4), 1, -x(2);...
%             -sec(xdot(3)^3*x(1))^2 - 1/(3*xdot(3)), 2/(3*xdot(3)), 0, 2*x(4)/(3*xdot(3))];
%     
%     dfdu = [2*(xdot(1) - u)];
%             0;...
%             0;...
%             0];
%     
%     dfdxdot = [-2*(xdot(1) - u), 0, 0, -10/x(1)];
%                0, 2*x(1)^2*xdot(2), -cos(xdot(3)*x(4)), 0;...
%                xdot(2)*xdot(3)*xdot(4), xdot(1)*xdot(3)*xdot(4), xdot(1)*xdot(2)*xdot(4), xdot(1)*xdot(2)*xdot(3);...
%                0, 0, -3*xdot(3)^2*x(1)/(cos(xdot(3)^3*x(1))) + -(x(4)^2 +2*x(2) - x(1))/(3*xdot(3)^2), 0];
%      
%% Crazy 1-D system with single input
%     f = [-10*xdot(1)*u/x(1)];
% 
%     dfdx = [10*xdot(1)*u/x(1)^2,];
% 
%     dfdu = [-10*xdot(1)/x(1)]; 
%     
%     dfdxdot = [-10*u/x(1)];
    
%% Less crazy 1-D system with single input
    f = [xdot(1)*xdot(1)*u/x(1) - x(1) + 1/x(1)];

    dfdx = [-xdot(1)*xdot(1)*u/x(1)^2 - 1 - 1/x(1)^2];

    dfdu = [xdot(1)*xdot(1)/x(1)]; 
    
    dfdxdot = [2*xdot(1)*u/x(1)];
    
    %% Same as last system but no states in denominator
%     f = [xdot(1)*xdot(1)*u - x(1)];
% 
%     dfdx = [-1];
% 
%     dfdu = [xdot(1)*xdot(1)]; 
%     
%     dfdxdot = [2*xdot(1)*u];
end