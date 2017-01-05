function [f, dfdx, dfdu, df_ddxdt] = vf(x, u, dxdt, params)

    n = params.n;
    m = params.m;
    
    [c1, c2, c3] = deal(params.Felast_consts(1), params.Felast_consts(2), params.Felast_consts(3));
    [c4, c5, c6] = deal(params.Melast_consts(1), params.Melast_consts(2), params.Melast_consts(3));
    
    [P0, gama0, r0, L0] = deal(params.x_rest(1), params.x_rest(2), params.x_rest(3), params.x_rest(4));
    [P, gama, r, L] = deal(x(1), x(2), x(3), x(4));
    [dP, dgama, dr, dL] = deal(dxdt(1), dxdt(2), dxdt(3), dxdt(4));

    
    f = [-dP + 0.5*(u - P);...  % The constant in front of (u-P) is arbitrary
         c2*dL + 2*c1*L*dL - pi*r^2*dP - 2*pi*P*r*dr + 2*pi*cot(gama)^2*r^2*dP + 4*pi*cot(gama)^2*P*r*dr - 4*pi*cot(gama)*P*r^2*dgama*(cot(gama)^2 + 1);...
%          c5*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*c4*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*pi*cot(gama)*r^3*dP + 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) - 6*pi*cot(gama)*P*r^2*dr;...  
         c5*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*c4*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) + 2*pi*cot(gama)*r^3*dP - 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) + 6*pi*cot(gama)*P*r^2*dr;...        % flipped the sign of M_elast in force balance equations
         sin(gama)*dgama + (cos(gama0)*dL)/L0];
 
         
 % -phi + (-tan(gama)*L/r + tan(gama0)*L0/r0)      this is the phi constraint (only here as a reminder for when I make it a constraint later
     
         
    dfdx = [-1*0.5,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                     0;...
            4*pi*dr*r*cot(gama)^2 - 2*pi*dr*r - 4*pi*dgama*r^2*cot(gama)*(cot(gama)^2 + 1),                                                                                                                                                                                                                                                                                                              4*pi*P*dgama*r^2*(cot(gama)^2 + 1)^2 - 4*pi*dP*r^2*cot(gama)*(cot(gama)^2 + 1) - 8*pi*P*dr*r*cot(gama)*(cot(gama)^2 + 1) + 8*pi*P*dgama*r^2*cot(gama)^2*(cot(gama)^2 + 1),                                                                                                                                                                                                                                                                                        4*pi*P*dr*cot(gama)^2 - 2*pi*dP*r - 2*pi*P*dr + 4*pi*dP*r*cot(gama)^2 - 8*pi*P*dgama*r*cot(gama)*(cot(gama)^2 + 1),                                                                                                                                                                                                                                               2*c1*dL;...
%             2*pi*dgama*r^3*(cot(gama)^2 + 1) - 6*pi*dr*r^2*cot(gama), c5*((dL*(tan(gama)^2 + 1))/r - (L*dr*(tan(gama)^2 + 1))/r^2 + (2*L*dgama*tan(gama)*(tan(gama)^2 + 1))/r) + 2*c4*((L*tan(gama))/r - (L0*tan(gama0))/r0)*((dL*(tan(gama)^2 + 1))/r - (L*dr*(tan(gama)^2 + 1))/r^2 + (2*L*dgama*tan(gama)*(tan(gama)^2 + 1))/r) + 2*pi*dP*r^3*(cot(gama)^2 + 1) + (2*L*c4*(tan(gama)^2 + 1)*((dL*tan(gama))/r + (L*dgama*(tan(gama)^2 + 1))/r - (L*dr*tan(gama))/r^2))/r + 6*pi*P*dr*r^2*(cot(gama)^2 + 1) - 4*pi*P*dgama*r^3*cot(gama)*(cot(gama)^2 + 1), 6*pi*P*dgama*r^2*(cot(gama)^2 + 1) - 2*c4*((L*tan(gama))/r - (L0*tan(gama0))/r0)*((dL*tan(gama))/r^2 + (L*dgama*(tan(gama)^2 + 1))/r^2 - (2*L*dr*tan(gama))/r^3) - 6*pi*dP*r^2*cot(gama) - 12*pi*P*dr*r*cot(gama) - c5*((dL*tan(gama))/r^2 + (L*dgama*(tan(gama)^2 + 1))/r^2 - (2*L*dr*tan(gama))/r^3) - (2*L*c4*tan(gama)*((dL*tan(gama))/r + (L*dgama*(tan(gama)^2 + 1))/r - (L*dr*tan(gama))/r^2))/r^2, c5*((dgama*(tan(gama)^2 + 1))/r - (dr*tan(gama))/r^2) + 2*c4*((dgama*(tan(gama)^2 + 1))/r - (dr*tan(gama))/r^2)*((L*tan(gama))/r - (L0*tan(gama0))/r0) + (2*c4*tan(gama)*((dL*tan(gama))/r + (L*dgama*(tan(gama)^2 + 1))/r - (L*dr*tan(gama))/r^2))/r;...
            - 2*dgama*pi*(cot(gama)^2 + 1)*r^3 + 6*dr*pi*cot(gama)*r^2, c5*((dL*(tan(gama)^2 + 1))/r - (L*dr*(tan(gama)^2 + 1))/r^2 + (2*L*dgama*tan(gama)*(tan(gama)^2 + 1))/r) + 2*c4*((L*tan(gama))/r - (L0*tan(gama0))/r0)*((dL*(tan(gama)^2 + 1))/r - (L*dr*(tan(gama)^2 + 1))/r^2 + (2*L*dgama*tan(gama)*(tan(gama)^2 + 1))/r) - 2*pi*dP*r^3*(cot(gama)^2 + 1) + (2*L*c4*(tan(gama)^2 + 1)*((dL*tan(gama))/r + (L*dgama*(tan(gama)^2 + 1))/r - (L*dr*tan(gama))/r^2))/r - 6*pi*P*dr*r^2*(cot(gama)^2 + 1) + 4*pi*P*dgama*r^3*cot(gama)*(cot(gama)^2 + 1), 6*pi*dP*r^2*cot(gama) - 2*c4*((L*tan(gama))/r - (L0*tan(gama0))/r0)*((dL*tan(gama))/r^2 + (L*dgama*(tan(gama)^2 + 1))/r^2 - (2*L*dr*tan(gama))/r^3) - c5*((dL*tan(gama))/r^2 + (L*dgama*(tan(gama)^2 + 1))/r^2 - (2*L*dr*tan(gama))/r^3) + 12*pi*P*dr*r*cot(gama) - 6*pi*P*dgama*r^2*(cot(gama)^2 + 1) - (2*L*c4*tan(gama)*((dL*tan(gama))/r + (L*dgama*(tan(gama)^2 + 1))/r - (L*dr*tan(gama))/r^2))/r^2, c5*((dgama*(tan(gama)^2 + 1))/r - (dr*tan(gama))/r^2) + 2*c4*((dgama*(tan(gama)^2 + 1))/r - (dr*tan(gama))/r^2)*((L*tan(gama))/r - (L0*tan(gama0))/r0) + (2*c4*tan(gama)*((dL*tan(gama))/r + (L*dgama*(tan(gama)^2 + 1))/r - (L*dr*tan(gama))/r^2))/r;...
            0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                        dgama*cos(gama),                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                     0];
   
    dfdu = [1*0.5;...
            0;...
            0;...
            0];
        
    df_ddxdt = [-1,                                                                                                                               0,                                                                                                           0,                                                                            0;...
                2*pi*r^2*cot(gama)^2 - pi*r^2,                                                                                         -4*pi*P*r^2*cot(gama)*(cot(gama)^2 + 1),                                                                             4*pi*P*r*cot(gama)^2 - 2*pi*P*r,                                                                  c2 + 2*L*c1;...
%                 -2*pi*r^3*cot(gama), 2*pi*P*r^3*(cot(gama)^2 + 1) + (L*c5*(tan(gama)^2 + 1))/r + (2*L*c4*((L*tan(gama))/r - (L0*tan(gama0))/r0)*(tan(gama)^2 + 1))/r, - 6*pi*P*r^2*cot(gama) - (L*c5*tan(gama))/r^2 - (2*L*c4*tan(gama)*((L*tan(gama))/r - (L0*tan(gama0))/r0))/r^2, (c5*tan(gama))/r + (2*c4*tan(gama)*((L*tan(gama))/r - (L0*tan(gama0))/r0))/r;...
                2*pi*r^3*cot(gama), (L*c5*(tan(gama)^2 + 1))/r - 2*pi*P*r^3*(cot(gama)^2 + 1) + (2*L*c4*((L*tan(gama))/r - (L0*tan(gama0))/r0)*(tan(gama)^2 + 1))/r, 6*pi*P*r^2*cot(gama) - (L*c5*tan(gama))/r^2 - (2*L*c4*tan(gama)*((L*tan(gama))/r - (L0*tan(gama0))/r0))/r^2, (c5*tan(gama))/r + (2*c4*tan(gama)*((L*tan(gama))/r - (L0*tan(gama0))/r0))/r;...
                0,                                                                                                                       sin(gama),                                                                                                           0,                                                                cos(gama0)/L0]; 
        
            
%     % Possible fixing of a sign error in derivatives 1/2/2017 --BAD
%     f = [1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 -1] * f;
%     dfdx = [1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 -1] * dfdx;
%     df_ddxdt = [1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 -1] * df_ddxdt;
    
    
    
%     % Just left-over from previous code. Left here as an example
%     % simple non-linear implicit dynamics: xdot^2 = -x*u
%     f = -dxdt - x*u;
%     dfdx = -u;
%     dfdu = -x;
%     df_ddxdt = -eye(n);
%      

end