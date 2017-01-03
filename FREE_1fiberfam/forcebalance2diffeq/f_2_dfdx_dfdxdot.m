%% Calculate dfdx, dfdxdot from f (only the second third row of f)

clear

syms r0 L0 gama0 P gama r L dP dgama dr dL c1 c2 c3 c4 c5 c6 x

x = [P gama r L];
xdot = [dP dgama dr dL];

%f3 = [c5*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) + 2*c4*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) + 2*pi*cot(gama)*r^3*dP - 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) + 6*pi*cot(gama)*P*r^2*dr];
f3 = [c5*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*c4*((L0*tan(gama0))/r0 - (tan(gama)*L)/r)*((tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2) - 2*pi*cot(gama)*r^3*dP + 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) - 6*pi*cot(gama)*P*r^2*dr];


df3dx = jacobian(f3, x);
df3dxdot = jacobian(f3, xdot);
