function f = test_f_dfdxdot(arg, params, x_const, u_const)

    n = params.n;
    m = params.m;
    
    if nargin == 2
        x_const = ones(n,1);
        u_const = ones(m,1);
    end

%% For testing dfdxdot
    x = x_const;    % x is fixed    
    u = u_const;    % u is fixed
    dxdt = arg;  

%% Set value of output    
    [f, ~, ~, ~] = vf(x, u, dxdt, params);

end