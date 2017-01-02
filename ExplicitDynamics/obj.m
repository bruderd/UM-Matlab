function [c, grad_c] = obj(s, params)

    % Extract values form 'params' struct
    N = params.N;
    n = params.n;
    m = params.m;
    dt = params.dt;
    
    [x, u] = state_decode(s, params);

    % Calculate running cost
    rc = 0;
    drcdx = zeros( (N + 1) , n );
    drcdu = zeros( (N + 1) , m);
    for k = 1:N+1
        [ic , dicdx, dicdu] = instant_cost(x(:,k), u(:,k));
        rc = rc + ic * dt;
        
% This is bad code that was commented out. state_encode takes care of
% indexing shit so this stuff was redundant and messing everything up!
%         drcdx( k, ( ( k - 1 ) * n + 1 ):( k*n ) ) = dicdx * dt;
%         drcdu( k, ( ( k - 1 ) * m + 1 ):( k*m ) ) = dicdu * dt;
        drcdx( k, : ) = dicdx * dt;
        drcdu( k, : ) = dicdu * dt;
    end
    
    % Calculate final cost
    [fc, dfcdx, dfcdu] = final_cost_none(x(:,end), u(:,end));
    
    % Set value for cost and gradient of cost
    c = rc + fc;
    
    %possible fix to problem...
    drcdx(end,:) = drcdx(end,:) + dfcdx;
    drcdu(end,:) = drcdu(end,:) + dfcdu;
    dcdx = drcdx;
    dcdu = drcdu;
    
%     dcdx = [drcdx ; dfcdx];
%     dcdu = [drcdu ; dfcdu];
    grad_c = state_encode(dcdx', dcdu', params);
    
    
    
end