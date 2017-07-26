function [c, grad_c] = obj(s, params)

    % Extract values form 'params' struct
    N = params.N;
    n = params.n;
    m = params.m;
    dt = params.dt;
    
    [x, u, dxdt] = state_decode(s, params);

    % Calculate running cost
    rc = 0;
    drcdx = zeros( (N + 1) , n );
    drcdu = zeros( (N + 1) , m);
    for k = 1:N+1
        [ic , dicdx, dicdu] = instant_cost(x(:,k), u(:,k), params);
        rc = rc + ic * dt;
        
        drcdx( k, : ) = dicdx * dt;
        drcdu( k, : ) = dicdu * dt;
    end
    
    % Calculate final cost
    [fc, dfcdx, dfcdu] = final_cost(x(:,end), u(:,end), params);
    
    % Set value for cost and gradient of cost
    c = rc + fc;
    
    %possible fix to problem...
    drcdx(end,:) = drcdx(end,:) + dfcdx;
    drcdu(end,:) = drcdu(end,:) + dfcdu;
    dcdx = drcdx;
    dcdu = drcdu;
    
%     dcdx = [drcdx ; dfcdx];
%     dcdu = [drcdu ; dfcdu];
    grad_c = state_encode(dcdx', dcdu', zeros(1, N*n), params);
    
    
    
end