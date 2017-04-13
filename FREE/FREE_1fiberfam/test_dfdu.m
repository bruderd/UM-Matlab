function dfdu = test_dfdu(arg, params, x_const, dxdt_const)

    n = params.n;
    m = params.m;
    
    if nargin == 2
        x_const = ones(n,1);
        dxdt_const = zeros(n,1);
    end
    
    % For testing dfdu
    x = x_const;    % x is fixed    
    u = arg;  
    dxdt = dxdt_const;  % dxdt if fixed

    [~, ~, dfdu, ~] = vf(x, u, dxdt, params);

end