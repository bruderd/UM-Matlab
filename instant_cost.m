function [ic, dicdx, dicdu] = instant_cost(x, u)

    ic = -sum(x)*sum(u);
    dicdx = -[sum(u) sum(u)];   %must make dimensions consistent, not just 1-D. Although for this special case it should still work
    dicdu = -[sum(x)];
    
%     ic = sum(x) + sum(u);
%     dicdx = 1;  %must make dimensions consistent, not just 1-D. Although for this special case it should still work
%     dicdu = 1;

%     ic = x(end)*u(end);
%     dicdx = zeros(1, length(x));
%     dicdx(end) = u(end);
%     dicdu = zeros(1, length(u));
%     dicdu(end) = x(end);
    
%     ic = x(end) + u(end);
%     dicdx = zeros(1, length(x));
%     dicdx(end) = 1;
%     dicdu = zeros(1, length(u));
%     dicdu(end) = 1;

end