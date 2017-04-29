function [fc, dfcdx, dfcdu] = final_cost(x, u, params)

%     fc = -10*(sum(x)*sum(u));
%     dfcdx = -[10*sum(u) 10+sum(u)];
%     dfcdu = -[10*sum(x)];
    
%     fc = 10*(sum(x));
%     dfcdx = [10 10 10];
%     dfcdu = [0 0 0];
    
%     fc = x(1)^2 + x(2)^2 + x(3)^2;
%     dfcdx = [2*x(1) 2*x(2) 2*x(3)];   %must make dimensions consistent
%     dfcdu = [0 0 0];

%     % For Shannon's LQR example
%     fc = 0;
%     dfcdx = [0 0];
%     dfcdu = [0];
    
%     % Dubin's Car
%     x_ = x(1);
%     y_ = x(2);
%     fc = x_^2 + y_^2;
%     dfcdx = [2*x_ 2*y_ 0];
%     dfcdu = [0 0];

    fc = u'*u + x'*x;
    dfcdx = 2*x';
    dfcdu = 2*u';

end