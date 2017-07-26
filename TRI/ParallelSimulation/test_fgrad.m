function fgrad = test_fgrad(arg, params, choice)

    m = params.numFREEs;
    
if nargin == 2
    choice = 1; % if not specified, choose to check dfdx
end
      
%% Choose which gradient is to be tested
if choice == 1
    x = arg;        
    u = ones(m,1);  % u is fixed
    xdot = ones(6,1);  % xdot is fixed
    [~, fgrad, ~, ~] = vf(x, u, xdot, params);
elseif choice == 2
    x = ones(6,1);    % x is fixed    
    u = ones(m,1);    % u is fixed
    xdot = arg;  
    [~, ~, ~, fgrad] = vf(x, u, xdot, params);
elseif choice == 3
    x = ones(6,1);      % x is fixed    
    u = arg;            % u is fixed
    xdot = zeros(6,1);   
    [~, ~, fgrad, ~] = vf(x, u, xdot, params);
end
                     

end