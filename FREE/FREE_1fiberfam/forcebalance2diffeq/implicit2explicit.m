%% implicit_2_explicit

syms P gama r L r0 L0 gama0 dP dgama dr dL c1 c2 c3 c4 c5 c6

assume(r0, 'real')
assumeAlso(L0, 'real')
assumeAlso(gama0, 'real')
assumeAlso(dP, 'real')
assumeAlso(dgama, 'real')
assumeAlso(dr, 'real')
assumeAlso(dL, 'real')
assumeAlso(c1, 'real')
assumeAlso(c2, 'real')
assumeAlso(c3, 'real')
assumeAlso(c4, 'real')
assumeAlso(c5, 'real')
assumeAlso(c6, 'real')     

assumeAlso(P, 'real')
assumeAlso(gama, 'real')
assumeAlso(r, 'real')
assumeAlso(L, 'real')

%% Version where F_elast = (L-L0), M_elast = theta (1/12/2017)
% f = [-dP + 0.5*(u - P);...  % The constant in front of (u-P) is arbitrary
%      dL + pi*r^2*dP + 2*pi*P*r*dr - 2*pi*cot(gama)^2*r^2*dP - 4*pi*cot(gama)^2*P*r*dr + 4*pi*cot(gama)*P*r^2*dgama*(cot(gama)^2 + 1);...
%      (tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2 - 2*pi*cot(gama)*r^3*dP + 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) - 6*pi*cot(gama)*P*r^2*dr;...
%      sin(gama)*dgama + (cos(gama0)*dL)/L0];

sol = solve([0 == dL + pi*r^2*dP + 2*pi*P*r*dr - 2*pi*cot(gama)^2*r^2*dP - 4*pi*cot(gama)^2*P*r*dr + 4*pi*cot(gama)*P*r^2*dgama*(cot(gama)^2 + 1),...
            0 == (tan(gama)*dL)/r + (L*(tan(gama)^2 + 1)*dgama)/r - (tan(gama)*L*dr)/r^2 - 2*pi*cot(gama)*r^3*dP + 2*pi*P*r^3*dgama*(cot(gama)^2 + 1) - 6*pi*cot(gama)*P*r^2*dr,...
            0 == sin(gama)*dgama + (cos(gama0)*dL)/L0],...
            [dgama, dr, dL]);

