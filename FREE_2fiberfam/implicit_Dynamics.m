function [con, coneq, grad_con, grad_coneq] = implicit_Dynamics(s, params)

    % Extract values form 'params' struct
    N = params.N;
    n = params.n;
    m = params.m;
    dt = params.dt;
    x0 = params.x0;
    u0 = params.u0;
    Pmax = params.Pmax;
    
    [x, u, dxdt] = state_decode(s, params);

    %% Calculate inequality constraints (con) and gradient of them (grad_con)
    con = zeros(n + 1, 2*(N+1));
    grad_con = zeros(n + 1, 2*N+1 + 1, (2*N+1)*n + (N+1)*m);
    
    for p = 1:N+1
        % Upper bounds on states
        con(1, p) = x(1,p) - Pmax;           %pressure limitation Pmax
        con(2, p) = x(2,p)*sign(x(2,p)) - 1.5;          %fiber angle cannot exceed 90 degrees in magnitude  
        con(3, p) = x(3,p)*sign(x(3,p)) - 1.5;          %fiber angle cannot exceed 90 degrees in magnitude  
        con(4, p) = x(4,p) - x0(4)*5;      %radius can triple, but no more
        con(5, p) = x(5,p) - (x0(5)*1.5);   %length cannot increase more than 50%
 
        % Lower bounds on states
        con(1, (N+1)+p) = -x(1,p) + 0;
        con(2, (N+1)+p) = -x(2,p)*sign(x(2,p)) + 0;            %fiber angle cannot cross over 0 degrees 
        con(3, (N+1)+p) = -x(3,p)*sign(x(3,p)) + 0;            %fiber angle cannot cross over 0 degrees   
        con(4, (N+1)+p) = -x(4,p) + x0(4)*0.5;           %radius cannot decrease in real life but ALLOW THIS FOR NOW (1/6/2017)
        con(5, (N+1)+p) = -x(5,p) + (x0(5)*0.5);     %length cannot decrease more than 50%
             
    end
    
    for q = 1:N+1
        % Gradient of upper bounds on states
        grad_con(1, q, n*(q-1)+1) = 1;
        grad_con(2, q, n*(q-1)+2) = 1*sign(x(2,q));
        grad_con(3, q, n*(q-1)+3) = 1*sign(x(3,q));
        grad_con(4, q, n*(q-1)+4) = 1;
        grad_con(5, q, n*(q-1)+5) = 1;

        % Gradient of lower bounds on states
        grad_con(1, (N+1)+q, n*(q-1)+1) = -1;
        grad_con(2, (N+1)+q, n*(q-1)+2) = -1*sign(x(2,q));
        grad_con(3, (N+1)+q, n*(q-1)+3) = -1*sign(x(3,q));
        grad_con(4, (N+1)+q, n*(q-1)+4) = -1;
        grad_con(5, (N+1)+q, n*(q-1)+5) = -1;
        
    end
    
    % Extra twist constraint, relaxed with fudge-factor
    fudge = 10;
    for j = 1:N+1
        % Constraint to ensure that the twist of the end-point is consistent for both fibers (upper and lower bounds, respectively)
        con(6, j) = ( x(5,j)/x(4,j) )*(tan(x(2,j)) - tan(x(3,j))) - ( x0(5)/x0(4) )*(tan(x0(2)) - tan(x0(3))) - fudge;
        con(6, (N+1)+j) = -(( x(5,j)/x(4,j) )*(tan(x(2,j)) - tan(x(3,j))) - ( x0(5)/x0(4) )*(tan(x0(2)) - tan(x0(3)))) - fudge;
    
        % Gradient of extra twist constraint
        grad_con(6, j, n*(j-1)+1) = 0;
        grad_con(6, j, n*(j-1)+2) = ( x(5,j)/x(4,j) )*sec(x(2,j))^2;
        grad_con(6, j, n*(j-1)+3) = -( x(5,j)/x(4,j) )*sec(x(3,j))^2;
        grad_con(6, j, n*(j-1)+4) = -( x(5,j)/x(4,j)^2 )*(tan(x(2,j)) - tan(x(3,j)));
        grad_con(6, j, n*(j-1)+5) = ( 1/x(4,j) )*(tan(x(2,j)) - tan(x(3,j)));   
        
        grad_con(6, (N+1)+j, n*(j-1)+1) = 0;
        grad_con(6, (N+1)+j, n*(j-1)+2) = -( x(5,j)/x(4,j) )*sec(x(2,j))^2;
        grad_con(6, (N+1)+j, n*(j-1)+3) = ( x(5,j)/x(4,j) )*sec(x(3,j))^2;
        grad_con(6, (N+1)+j, n*(j-1)+4) = ( x(5,j)/x(4,j)^2 )*(tan(x(2,j)) - tan(x(3,j)));
        grad_con(6, (N+1)+j, n*(j-1)+5) = -( 1/x(4,j) )*(tan(x(2,j)) - tan(x(3,j)));   
    end
    
    
    %% Calculate equality constraints (coneq) and gradient (grad_coneq)
    coneq = zeros(n + 0, 2*N+1);
    grad_coneq = zeros(n + 0, (2*N+1), (2*N+1)*n + (N+1)*m);

    % IC equality constraint
    coneq(1:n , 1) = x(:,1) - x0;
    for j = 1:n
        grad_coneq(j, 1 , j) = 1;
    end
    % Dynamics as equality constraints
    for k = 1:N
        
        [f, dfdx, dfdu, df_ddxdt] = vf2(x(:,k), u(:,k), dxdt(:,k), params);
        
        [h, dh_ddxdt, dhdx_k, dhdx_kplus1, dhdu] = vh(dxdt(:,k), x(:,k), x(:,k+1), u(:,k), params);
        
        coneq(1:n, k+1) =  f;
        coneq(1:n, (N+1) + k) = h;
       
        grad_coneq(1:n, k+1, n*(k-1) + 1 : n*(k-1) + n) = dfdx;
     %   grad_coneq( 1:n, k+1, n*k + 1 : n*k + n) = -1*eye(n);
        grad_coneq(1:n, k+1, (N+1)*n + m*(k-1) + 1 : (N+1)*n + m*(k-1) + m) = dfdu;
        grad_coneq(1:n, k+1, (N+1)*n + (N+1)*m + (k-1)*n + 1 : (N+1)*n + (N+1)*m + (k-1)*n + n) = df_ddxdt;
        
        grad_coneq(1:n, (N+1) + k, n*(k-1) + 1 : n*(k-1) + n) = dhdx_k;
        grad_coneq(1:n, (N+1) + k, n*k + 1 : n*k + n) = dhdx_kplus1;
        grad_coneq(1:n, (N+1) + k, (N+1)*n + m*(k-1) + 1 : (N+1)*n + m*(k-1) + m) = dhdu;
        grad_coneq(1:n, (N+1) + k, (N+1)*n + (N+1)*m + (k-1)*n + 1 : (N+1)*n + (N+1)*m + (k-1)*n + n) = dh_ddxdt;
      
    end
    
%     for j = 1:N
%         % Constraint to ensure that the twist of the end-point is consistent for both fibers
%         coneq(6, j) = ( x(5,j)/x(4,j) )*(tan(x(2,j)) - tan(x(3,j))) - ( x0(5)/x0(4) )*(tan(x0(2)) - tan(x0(3)));
%         coneq(6, (N+1)+j) = 0;
%     
%         % Gradient of extra twist constraint
%         grad_coneq(6, j, n*(j-1)+1) = 0;
%         grad_coneq(6, j, n*(j-1)+2) = ( x(5,j)/x(4,j) )*sec(x(2,j))^2;
%         grad_coneq(6, j, n*(j-1)+3) = -( x(5,j)/x(4,j) )*sec(x(3,j))^2;
%         grad_coneq(6, j, n*(j-1)+4) = -( x(5,j)/x(4,j)^2 )*(tan(x(2,j)) - tan(x(3,j)));
%         grad_coneq(6, j, n*(j-1)+5) = ( 1/x(4,j) )*(tan(x(2,j)) - tan(x(3,j)));   
%     end

    %% Converts constraints into form that fmincon likes
    coneq = matrix2vector(coneq);
    grad_coneq = matrix3D_2_matrix2D(grad_coneq);
    
    con = matrix2vector(con);
    grad_con = matrix3D_2_matrix2D(grad_con);
    
end