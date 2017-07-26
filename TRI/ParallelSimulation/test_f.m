function f = test_f(arg, params, choice)

    m = params.numFREEs;
    
if nargin == 2
    choice = 1; % if not specified, choose to check dfdx
end
    
    
%% Choose which gradient is to be tested
if choice == 1
    x = arg;        
    u = ones(m,1);  % u is fixed
    xdot = ones(6,1);  % xdot is fixed
elseif choice == 2
    x = ones(6,1);    % x is fixed    
    u = ones(m,1);    % u is fixed
    xdot = arg;  
elseif choice == 3
    x = ones(6,1);    % x is fixed      
    u = arg;          
    xdot = zeros(6,1);  % xdot is fixed    
end
                     
%% Set value of output    
    [f, ~, ~, ~] = vf(x, u, xdot, params);

end