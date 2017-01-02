function s = state_encode(x, u, params)
% converts states and inputs ICs into decision variable s
    x0 = params.x0;
    u0 = params.u0;
    N = params.N;
    m = params.m;
    n = params.n;

    s = zeros((N+1)*n + (N+1)*m,1);
    for k = 0 : N
        s(k*n+1:k*n+n) = x(:,k+1)';
    end
        
    for j = 0 : N
        s(k*n+n+j*m+1 : k*n+n+j*m+m) = u(:,j+1)';
    end
       
%    s = s';
end