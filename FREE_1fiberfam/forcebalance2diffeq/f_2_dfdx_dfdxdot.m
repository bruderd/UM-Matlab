%% Calculate dfdx, dfdxdot from f (only the second & third row of f)

clear

syms r0 L0 gama0 P gama r L dP dgama dr dL c1 c2 c3 c4 c5 c6 x

x = [P gama r L];
xdot = [dP dgama dr dL];

%f3 = [c5*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) + 2*c4*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) + 2*pi*cot(gama)*r^3*dP - 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) + 6*pi*cot(gama)*P*r^2*dr];
%f3 = [c5*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*c4*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*pi*cot(gama)*r^3*dP + 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) - 6*pi*cot(gama)*P*r^2*dr];
%f3 = [c5*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*c4*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) + 2*pi*cot(gama)*r^3*dP - 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) + 6*pi*cot(gama)*P*r^2*dr];

% f2 = [c2*dL + 2*c1*L*dL + pi*r^2*dP + 2*pi*P*r*dr - 2*pi*cot(gama)^2*r^2*dP - 4*pi*cot(gama)^2*P*r*dr + 4*pi*cot(gama)*P*r^2*dgama*(cot(gama)^2 + 1)];
% f3 = [c5*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*c4*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*pi*cot(gama)*r^3*dP + 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) - 6*pi*cot(gama)*P*r^2*dr];

%% 1/6/2017: Made it so F_elast is a function of (L-L0)
% f2 = c2*dL - 2*c1*dL*(L0 - L) + pi*r^2*dP + 2*pi*P*r*dr - 2*pi*cot(gama)^2*r^2*dP - 4*pi*cot(gama)^2*P*r*dr + 4*pi*cot(gama)*P*r^2*dgama*(cot(gama)^2 + 1);
% f3 = c5*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*c4*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*pi*cot(gama)*r^3*dP + 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) - 6*pi*cot(gama)*P*r^2*dr;                                                                                                                                                                                                                                                                                                                                                                                                                 
% f4 = sin(gama)*dgama + (cos(gama0)*dL)/L0;
%  
%% 1/7/2017: Got rid of a minus sign in front of phi that I had added earlier
f2 = c2*dL + 2*c1*L*dL + pi*r^2*dP + 2*pi*P*r*dr - 2*pi*cot(gama)^2*r^2*dP - 4*pi*cot(gama)^2*P*r*dr + 4*pi*cot(gama)*P*r^2*dgama*(cot(gama)^2 + 1);
f3 = 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) - 2*c4*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*pi*cot(gama)*r^3*dP - c5*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 6*pi*cot(gama)*P*r^2*dr;
f4 = sin(gama)*dgama + (cos(gama0)*dL)/L0;


%% Calculate gradients
df2dx = jacobian(f2, x);
df2dxdot = jacobian(f2, xdot);

df3dx = jacobian(f3, x);
df3dxdot = jacobian(f3, xdot);

df4dx = jacobian(f4, x);
df4dxdot = jacobian(f4, xdot);
