function [f, dfdx, dfdu, df_ddxdt] = vf2(x, u, dxdt, params)

    n = params.n;
    m = params.m;
    
    [c1, c2, c3] = deal(params.Felast_consts(1), params.Felast_consts(2), params.Felast_consts(3));
    [c4, c5, c6] = deal(params.Melast_consts(1), params.Melast_consts(2), params.Melast_consts(3));
    
    [P0, gama0, betta0, r0, L0] = deal(params.x_rest(1), params.x_rest(2), params.x_rest(3), params.x_rest(4), params.x_rest(5));
    [P, gama, betta, r, L] = deal(x(1), x(2), x(3), x(4), x(5));
    [dP, dgama, dbetta, dr, dL] = deal(dxdt(1), dxdt(2), dxdt(3), dxdt(4), dxdt(5));
    
    
%% 2 fibers! F_elast = [c1 c2 c3]*[L^2 L 1]?; M_elast = [c4 c5 c6]*[theta^2 theta 1];  (1/17/2017)     
% f = [-dP + 0.5*(u - P);...  % The constant in front of (u-P) is arbitrary 
%      c2*dL + 2*c1*L*dL + pi*r^2*dP + 2*pi*P*r*dr - 2*pi*r^2*(cot(betta)^2 + cot(gama)^2)*dP + 2*pi*P*r^2*(2*cot(betta)*dbetta*(cot(betta)^2 + 1) + 2*cot(gama)*dgama*(cot(gama)^2 + 1)) - 4*pi*P*r*(cot(betta)^2 + cot(gama)^2)*dr;...
%      2*pi*r^3*dP*(cot(betta) + cot(gama)) - 2*pi*P*r^3*((cot(betta)^2 + 1)*dbetta + (cot(gama)^2 + 1)*dgama) + (c5*tan(gama)*dL)/r + (2*c4*tan(gama)^2*L*dL)/r^2 + 6*pi*P*r^2*dr*(cot(betta) + cot(gama)) + (c5*L*(tan(gama)^2 + 1)*dgama)/r - (2*c4*tan(gama)^2*L^2*dr)/r^3 - (c5*tan(gama)*L*dr)/r^2 + (2*c4*tan(gama)*L^2*(tan(gama)^2 + 1)*dgama)/r^2;...
%      sin(gama)*dgama + (cos(gama0)*dL)/L0;...
%      sin(betta)*dbetta + (cos(betta0)*dL)/L0];
%  
% dfdx = [0];     % not needed for ode15i 
% dfdu = [0];     % not needed for ode15i
% df_ddxdt = [0]; % not needed for ode15i


%% 2 fiber system: F_elast = c1*(L0-L); M_elast = c4*(-1)*phi     (1/18/2017)
% f = [-dP + 0.5*(u - P);...  % The constant in front of (u-P) is arbitrary 
%      pi*r^2*dP - c1*dL + 2*pi*P*r*dr - 2*pi*r^2*(cot(betta)^2 + cot(gama)^2)*dP + 2*pi*P*r^2*(2*cot(betta)*dbetta*(cot(betta)^2 + 1) + 2*cot(gama)*dgama*(cot(gama)^2 + 1)) - 4*pi*P*r*(cot(betta)^2 + cot(gama)^2)*dr;...
%      c4*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*pi*P*r^3*((cot(betta)^2 + 1)*dbetta + (cot(gama)^2 + 1)*dgama) + 2*pi*r^3*dP*(cot(betta) + cot(gama)) + 6*pi*P*r^2*dr*(cot(betta) + cot(gama));...
%      sin(gama)*dgama + (cos(gama0)*dL)/L0;...
%      sin(betta)*dbetta + (cos(betta0)*dL)/L0];
%  
% dfdx = [0];     % not needed for ode15i 
% dfdu = [0];     % not needed for ode15i
% df_ddxdt = [0]; % not needed for ode15i
 

%% 2 fiber system: F_elast = c1*(L0-L); M_elast = c4*(-1)*phi. Added 1/2 in front of gama-betta summation terms    (1/18/2017)
% f = [-dP + 0.5*(u - P);...  % The constant in front of (u-P) is arbitrary 
%      pi*r^2*dP - c1*dL + 2*pi*P*r*dr - pi*r^2*(cot(betta)^2 + cot(gama)^2)*dP + pi*P*r^2*(2*cot(betta)*dbetta*(cot(betta)^2 + 1) + 2*cot(gama)*dgama*(cot(gama)^2 + 1)) - 2*pi*P*r*(cot(betta)^2 + cot(gama)^2)*dr;...
%      c4*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - pi*P*r^3*((cot(betta)^2 + 1)*dbetta + (cot(gama)^2 + 1)*dgama) + pi*r^3*dP*(cot(betta) + cot(gama)) + 3*pi*P*r^2*dr*(cot(betta) + cot(gama));...
%      sin(gama)*dgama + (cos(gama0)*dL)/L0;...
%      sin(betta)*dbetta + (cos(betta0)*dL)/L0];
%      
% dfdx = [0];     % not needed for ode15i 
% dfdu = [0];     % not needed for ode15i
% df_ddxdt = [0]; % not needed for ode15i


%% 2 fiber system: F_elast = c1*(L0-L), M_elast = c4*(-1)*phi. Added 1/2 in front of gama-betta summation terms. Added additional constraint (phi_gama = phi_betta)   (1/23/2017)
f = [-dP + 0.5*(u - P);...  % The constant in front of (u-P) is arbitrary 
    pi*r^2*dP - c1*dL + 2*pi*P*r*dr - pi*r^2*(cot(betta)^2 + cot(gama)^2)*dP + pi*P*r^2*(2*cot(betta)*dbetta*(cot(betta)^2 + 1) + 2*cot(gama)*dgama*(cot(gama)^2 + 1)) - 2*pi*P*r*(cot(betta)^2 + cot(gama)^2)*dr;...
    c4*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - pi*P*r^3*((cot(betta)^2 + 1)*dbetta + (cot(gama)^2 + 1)*dgama) + pi*r^3*dP*(cot(betta) + cot(gama)) + 3*pi*P*r^2*dr*(cot(betta) + cot(gama));...
    sin(gama)*dgama + (cos(gama0)*dL)/L0;...
    sin(betta)*dbetta + (cos(betta0)*dL)/L0 + (tan(betta)*dL)/r - (tan(gama)*dL)/r + (L*(tan(betta)^2 + 1)*dbetta)/r - (L*(tan(gama)^2 + 1)*dgama)/r - (tan(betta)*L*dr)/r^2 + (tan(gama)*L*dr)/r^2];
%     (tan(betta)*dL)/r - (tan(gama)*dL)/r + (L*(tan(betta)^2 + 1)*dbetta)/r - (L*(tan(gama)^2 + 1)*dgama)/r - (tan(betta)*L*dr)/r^2 + (tan(gama)*L*dr)/r^2];

dfdx = [0];     % not needed for ode15i 
dfdu = [0];     % not needed for ode15i
df_ddxdt = [0]; % not needed for ode15i

end