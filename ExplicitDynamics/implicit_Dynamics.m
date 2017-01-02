function [con, coneq, grad_con, grad_coneq] = implicit_Dynamics(s, params)

    % Extract values form 'params' struct
    N = params.N;
    n = params.n;
    m = params.m;
    dt = params.dt;
    x0 = params.x0;
    u0 = params.u0;
    
    [x, u] = state_decode(s, params);

    %% Calculate inequality constraints (con) and gradient of them (grad_con)
    con = [];
    grad_con = [];
    
    %% Calculate equality constraints (coneq) and gradient (grad_coneq)
    coneq = zeros(n, 2*N+1);
    grad_coneq = zeros(n, (2*N+1), (2*N+1)*(n+m));

    % IC equality constraint
    coneq(1:n , 1) = x(:,1) - x0;
    for j = 1:n
        grad_coneq(j, 1 , j) = 1;
    end
    % Dynamics as equality constraints
    for k = 1:N
        
        [f, dfdx, dfdu] = vf(x(:,k), u(:,k), params);
        
        coneq(1:n, k+1) =  - x(:,k+1) + x(:,k) + f*dt;
       
        grad_coneq( 1:n, k+1, n*(k-1) + 1 : n*(k-1) + n) = eye(n) + dfdx*dt;
        grad_coneq( 1:n, k+1, n*k + 1 : n*k + n) = -1*eye(n);
        grad_coneq( 1:n, k+1, (N+1)*n + m*(k-1) + 1 : (N+1)*n + m*(k-1) + m) = dfdu*dt;
      
    end

    %% Converts constraints into form that fmincon likes
    coneq = matrix2vector(coneq);
    grad_coneq = matrix3D_2_matrix2D(grad_coneq);
    
end