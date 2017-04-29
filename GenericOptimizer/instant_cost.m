function [ic, dicdx, dicdu] = instant_cost(x, u, params)

%     ic = -sum(x)*sum(u);
%     dicdx = -[sum(u) sum(u)];   %must make dimensions consistent, not just 1-D. Although for this special case it should still work
%     dicdu = -[sum(x)];
    
%     ic = x(1)^2 + x(2)^2 + x(3)^2;
%     dicdx = [2*x(1) 2*x(2) 2*x(3)];   %must make dimensions consistent
%     dicdu = [0 0 0];
    
%     % for Shannon's LQR problem
%     R = [1 0; 0 1];
%     Q = [1];
%     ic = x*R*x' + u'*Q*u;
%     dicdx = 2*x*R;
%     dicdu = 2*u*Q;
    
    % Dubin's Car
    ic = u'*u + x'*x;
    dicdx = 2*x';
    dicdu = 2*u';
    
    

end