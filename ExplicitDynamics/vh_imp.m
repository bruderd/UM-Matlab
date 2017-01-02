function [h, dhdx, dhdu] = vh_imp(x_k,x_kplus1,u,dxdt,params)

    dt = params.dt;

    h = -dxdt*dt + x_kplus1 - x_k;
    

end