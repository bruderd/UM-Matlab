function [con, coneq, grad_con, grad_coneq] = dynamics_LQR(s, params)

    % Extract values form 'params' struct
    N = params.N;
    n = params.n;
    m = params.m;
    dt = params.dt;
    x0 = params.x0;
    u0 = params.u0;
    A = params.A;
    B = params.B;
    C = params.C;
    D = params.D;
    
    [x, u] = state_decode(s, params);

    % Calculate inequality constraints (con)
    con = [];
    
    % Calculate gradient of inequality constraints (grad_con)
    grad_con = [];
    
    % Calculate equality constraints (coneq) and gradient (grad_coneq) for
    % the system with dynamics given by xdot = -I*x + u (assumes u has same
    % dimension as x)
    coneq = zeros(N+1, n);
    grad_coneq = zeros((N+1), n, (N+1)*(n+m));

    % IC equality constraint
%     coneq(1 , 1:n) = x(1,:) - x0' + u(1,:) - u0';
    coneq(1 , 1:n) = x(1,:) - x0';
    for j = 1:n
        grad_coneq(1, j , j) = 1;
%         grad_coneq(1, j , (N+1)*n + j) = 1;
    end
    % Dynamics as equality constraints
    for k = 1:N
        
        [f, dfdx, dfdu] = vf(x(k,:), u(k,:), params);
        
        coneq(k+1 , 1:n) =  - x(k+1,:)' + x(k,:)' + ( f )*dt;
        
        %TODO: MUST FIX THIS CHUNK OF CODE, THEN YOU CAN TRY TO RUN    
        %probably some major issues here.....
        grad_coneq(k+1, 1:n , n*(k-1) + 1 : n*(k-1) + n) = eye(n) + dfdx*dt;
        grad_coneq(k+1, 1:n , n*k + 1 : n*k + n) = -1*eye(n);
        grad_coneq(k+1, 1:n , (N+1)*n + m*(k-1) + 1 : (N+1)*n + m*(k-1) + n) = dfdu*dt;
        
%         for j = 1:n
%             grad_coneq(k+1, j , n*(k-1) + j) = 1 + dt;
%             grad_coneq(k+1, j , n*k + j) = -1;
%             grad_coneq(k+1, j , (N+1)*n + m*(k-1) + j) = dt;
%         end
      
    end

    coneq = matrix2vector(coneq);
    grad_coneq = matrix3D_2_matrix2D(grad_coneq);

    
end