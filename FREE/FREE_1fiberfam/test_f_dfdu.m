function f = test_f_dfdu(arg, params, x_const, dxdt_const)

    n = params.n;
    m = params.m;
    
    if nargin == 2
        x_const = ones(n,1);
        dxdt_const = zeros(n,1);
    end
    
%% For testing dfdu
    x = x_const;        
    u = arg;            % u is fixed
    dxdt = dxdt_const;  % dxdt is fixed
     

%% Set value of output    
    [f, ~, ~, ~] = vf(x, u, dxdt, params);

end