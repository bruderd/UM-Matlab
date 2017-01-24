%% Calculate dfdx, dfdxdot from f (only the second & third row of f)

clear

syms r0 L0 gama0 betta0 P gama betta r L dP dgama dbetta dr dL c1 c2 c3 c4 c5 c6 x u

x = [P gama betta r L];
xdot = [dP dgama dbetta dr dL];

%% 2 families of fibers (1/24/2017)
f = [-dP + 0.5*(u - P);...  % The constant in front of (u-P) is arbitrary 
     pi*r^2*dP - c1*dL + 2*pi*P*r*dr - pi*r^2*(cot(betta)^2 + cot(gama)^2)*dP + pi*P*r^2*(2*cot(betta)*dbetta*(cot(betta)^2 + 1) + 2*cot(gama)*dgama*(cot(gama)^2 + 1)) - 2*pi*P*r*(cot(betta)^2 + cot(gama)^2)*dr;...
     c4*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - pi*P*r^3*((cot(betta)^2 + 1)*dbetta + (cot(gama)^2 + 1)*dgama) + pi*r^3*dP*(cot(betta) + cot(gama)) + 3*pi*P*r^2*dr*(cot(betta) + cot(gama));...
     sin(gama)*dgama + (cos(gama0)*dL)/L0;...
     sin(betta)*dbetta + (cos(betta0)*dL)/L0];

f1 = f(1);
f2 = f(2);
f3 = f(3);
f4 = f(4);
f5 = f(5);


%% Calculate gradients
df1dx = jacobian(f1, x);
df1dxdot = jacobian(f1, xdot);

df2dx = jacobian(f2, x);
df2dxdot = jacobian(f2, xdot);

df3dx = jacobian(f3, x);
df3dxdot = jacobian(f3, xdot);

df4dx = jacobian(f4, x);
df4dxdot = jacobian(f4, xdot);

df5dx = jacobian(f5, x);
df5dxdot = jacobian(f5, xdot);

% Put these together into dfdx, dfdxdot
dfdx = [df1dx; df2dx; df3dx; df4dx; df5dx];
dfdxdot = [df1dxdot; df2dxdot; df3dxdot; df4dxdot; df5dxdot];