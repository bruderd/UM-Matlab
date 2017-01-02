function f = test_f(arg, params)

    n = params.n;
    m = params.m;
    
%% For testing dfdx
    x = arg;        
%     u = ones(m,1);  % u is fixed
%    dxdt = zeros(n,1);  % dxdt if fixed
    u = [1]';  % u is fixed
    dxdt = [1 5 7 3 7]';  % dxdt is fixed

% %% For testing dfdu
%     x = params.x0;    % x is fixed    
%     u = arg;  
%     dxdt = zeros(n,1);  % dxdt is fixed
     
% %% For testing dfdxdot
%     x = params.x0;    % x is fixed    
%     u = ones(m,1);  
%     dxdt = arg;  % dxdt is fixed

%% Set value of output    
    [f, ~, ~, ~] = vf(x, u, dxdt, params);

end