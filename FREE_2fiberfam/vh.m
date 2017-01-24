function [h, dh_ddxdt, dhdx_k, dhdx_kplus1, dhdu] = vh(dxdt, x_k,x_kplus1,u,params)
% h is the constraint that imposes that the decision variable dxdt is the
% derivative of x

    n = params.n;
    m = params.m;

    dt = params.dt;

    h = -dxdt*dt + x_kplus1 - x_k;
    
    dh_ddxdt = -eye(n)*dt;
    
    dhdx_k = -eye(n);
    
    dhdx_kplus1 = eye(n);
    
    dhdu = zeros(n,m);
    
end